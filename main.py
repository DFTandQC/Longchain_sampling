#!/usr/bin/env python3
"""
Multi-Molecule Cluster Sampling Pipeline - Main Entry Point

Run the full sampling workflow: generate multi-frame clusters with sampling engine,
then split output frames into individual seed XYZ files.

Supports both single-molecule (legacy) and multi-molecule cluster configurations.

Usage (Legacy - single molecule):
    python main.py --nseeds 100
    python main.py --nseeds 50 --N 3 --dmin 8.0 --dmax 12.0

Usage (Multi-molecule via JSON config):
    python main.py --config cluster_config.json --nseeds 100

Usage (Multi-molecule via command line):
    python main.py --molecules '[{"name":"PT","file":"opt-PT.xyz","count":2},{"name":"H2SO4","file":"h2so4.xyz","count":1}]' --nseeds 100

Usage (List available monomers):
    python main.py --list-monomers

Usage (Select molecules and override counts):
    python main.py --use "PT,H2SO4,NO2" --nseeds 100
    python main.py --use "PT,H2SO4,NO,NO2" --counts "PT=2,H2SO4=1,NO=2,NO2=1" --nseeds 200

Usage (Parallel sampling - run 5 independent jobs):
    python main.py --nseeds 100 --parallel-jobs 5
    python main.py --config config.json --nseeds 100 --parallel-jobs 3
"""
from pathlib import Path
import subprocess
import sys
import glob
import argparse
import shutil

from lib.config import ClusterConfig, get_argument_parser, load_molecules_from_args, list_available_monomers
from lib.extractor import split_multiframe
import importlib
import concurrent.futures
import os


def run_sampling_engine(args_list, sampling_workers=1):
    """Invoke the sampling engine (lib/sampling.py) as a module."""
    # Try to run the sampling engine in-process to avoid Python startup overhead.
    try:
        mod = importlib.import_module("lib.sampling")
        print(f"[MAIN] Running sampling engine in-process (workers={sampling_workers})")
        # Temporarily set sampling_workers as attribute for run_with_args to access
        import sys as sys_module
        sys_module.argv = ["lib.sampling"] + args_list
        
        # Parse and inject sampling_workers into args
        from lib.config import get_argument_parser
        ap = get_argument_parser()
        args = ap.parse_args(args_list)
        args.sampling_workers = sampling_workers
        
        # Call run_with_args with parsed args object + sampling_workers
        mod.run_with_args_parsed(args)
        return 0
    except Exception as e:
        print(f"[MAIN] In-process run failed, falling back to subprocess: {e}")
        cmd = [sys.executable, "-m", "lib.sampling"] + args_list
        if sampling_workers > 1:
            cmd.extend(["--sampling-workers", str(sampling_workers)])
        print("[MAIN] Running sampling engine as subprocess:", " ".join(cmd))
        res = subprocess.run(cmd)
        return res.returncode


def run_single_parallel_job(job_id, base_args, output_dir, base_seed):
    """Run a single sampling job for parallel mode"""
    try:
        # Create output directory
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Build command with modified output directory
        cmd = [sys.executable, "main.py"]
        
        # Add all arguments except --parallel-jobs, --out, and --seed
        skip_next = False
        for i, arg in enumerate(base_args):
            if skip_next:
                skip_next = False
                continue
            if arg == "--parallel-jobs":
                skip_next = True
                continue
            if arg == "--out":
                skip_next = True
                continue
            if arg == "--seed":
                skip_next = True
                continue
            cmd.append(arg)
        
        # Add unique seed for this parallel job
        unique_seed = base_seed + job_id * 100000
        cmd.extend(["--seed", str(unique_seed)])
        
        # Add output directory
        cmd.extend(["--out", output_dir])
        
        print(f"\n[PARALLEL Job {job_id}] Starting in: {output_dir}")
        print(f"[PARALLEL Job {job_id}] Command: {' '.join(cmd)}")
        
        # Run subprocess
        result = subprocess.run(cmd, capture_output=False, text=True)
        
        if result.returncode == 0:
            print(f"[PARALLEL Job {job_id}] ✓ Completed successfully")
            return True
        else:
            print(f"[PARALLEL Job {job_id}] ✗ Failed with return code {result.returncode}")
            return False
            
    except Exception as e:
        print(f"[PARALLEL Job {job_id}] ✗ Exception: {e}")
        return False


def merge_parallel_seeds_xyz(base_out, num_jobs):
    """
    Merge seeds_xyz from all parallel jobs into a unified folder structure.
    Preserves original file structure with job-specific subfolders.
    """
    base_out = Path(base_out)
    merged_dir = base_out / "seeds_xyz_combined"
    
    print("\n" + "="*70)
    print("MERGING PARALLEL JOBS - seeds_xyz")
    print("="*70)
    
    merged_dir.mkdir(parents=True, exist_ok=True)
    
    total_files = 0
    
    for job_id in range(1, num_jobs + 1):
        job_folder = base_out / f"sampling_{job_id:02d}"
        source_dir = job_folder / "seeds_xyz"
        
        if not source_dir.exists():
            print(f"[MERGE] Job {job_id}: No seeds_xyz found (skipping)")
            continue
        
        # Create job-specific subfolder in merged directory
        dest_dir = merged_dir / f"sampling_{job_id:02d}"
        dest_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy all files from this job's seeds_xyz
        xyz_files = list(source_dir.glob("*.xyz"))
        
        for src_file in xyz_files:
            dest_file = dest_dir / src_file.name
            try:
                shutil.copy2(src_file, dest_file)
                total_files += 1
            except Exception as e:
                print(f"[MERGE] Error copying {src_file.name}: {e}")
        
        print(f"[MERGE] Job {job_id}: Copied {len(xyz_files)} files to {dest_dir.name}/")
    
    print(f"\n[MERGE] ✓ Total files merged: {total_files}")
    print(f"[MERGE] Combined output: {merged_dir.resolve()}")
    print("="*70)


def main():
    ap = argparse.ArgumentParser(
        description="Generate multi-molecule cluster seeds: sample geometries then extract into individual files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples (Legacy - single molecule):
  python main.py --nseeds 100              # Generate 100 seeds (default N=2 PT)
  python main.py --nseeds 50 --N 3         # Generate 50 trimers
  python main.py --nseeds 20 --dmin 8 --dmax 12  # Custom distance range

Examples (Multi-molecule):
  python main.py --config config.json --nseeds 100
  python main.py --use "PT,H2SO4,NO2" --nseeds 100
  python main.py --use "PT,H2SO4,NO,NO2" --counts "PT=2,H2SO4=1,NO=2,NO2=1" --nseeds 200
  python main.py --list-monomers

Examples (Parallel sampling):
  python main.py --nseeds 100 --parallel-jobs 5          # 5 parallel jobs (default)
  python main.py --config config.json --nseeds 100 --parallel-jobs 3
  python main.py --use "PT,H2SO4" --nseeds 50 --parallel-jobs 10
        """
    )
    
    # Reuse config defaults
    cfg_def = ClusterConfig()
    ap.add_argument("--list-monomers", action="store_true",
                    help="List all available monomers in monomer/ directory and exit")
    ap.add_argument("--config", type=str, default="config.json",
                    help="JSON config file for multi-molecule clusters (default: config.json with presets)")
    ap.add_argument("--molecules", type=str, default=None,
                    help="Molecule spec as JSON string (overrides --mono and --N)")
    ap.add_argument("--use", type=str, default=None,
                    help="Specify which molecules to sample (comma-separated, e.g., \"PT,H2SO4,NO2\")")
    ap.add_argument("--counts", type=str, default=None,
                    help="Override molecule counts (e.g., \"PT=2,H2SO4=1,NO=3\")")
    ap.add_argument("--mono", default=cfg_def.mono_file,
                    help="[Legacy] Monomer xyz file (default: %(default)s)")
    ap.add_argument("--out", default=cfg_def.out_dir,
                    help="Output directory (default: %(default)s)")
    ap.add_argument("--nseeds", type=int, default=cfg_def.nseeds,
                    help="Number of cluster seeds to generate (default: %(default)s)")
    ap.add_argument("--N", type=int, default=cfg_def.N,
                    help="[Legacy] Monomers per cluster (default: %(default)s)")
    ap.add_argument("--clash_cut", type=float, default=cfg_def.clash_cut,
                    help="Clash detection threshold Angstrom (default: %(default)s)")
    ap.add_argument("--max_trials_add", type=int, default=cfg_def.max_trials_add,
                    help="Max placement attempts per monomer (default: %(default)s)")
    ap.add_argument("--dmin", type=float, default=cfg_def.dmin,
                    help="Min COM distance Angstrom (default: %(default)s)")
    ap.add_argument("--dmax", type=float, default=cfg_def.dmax,
                    help="Max COM distance Angstrom (default: %(default)s)")
    ap.add_argument("--lateral", type=float, default=cfg_def.lateral,
                    help="Max lateral offset Angstrom (default: %(default)s)")
    ap.add_argument("--jitter_deg", type=float, default=cfg_def.jitter_deg,
                    help="Angle jitter degrees (default: %(default)s)")
    ap.add_argument("--seed", type=int, default=cfg_def.random_seed,
                    help="Random seed for reproducibility (default: %(default)s)")
    ap.add_argument("--workers", type=int, default=1,
                    help="Worker threads for extraction (default: %(default)s)")
    ap.add_argument("--sampling-workers", type=int, default=1,
                    help="Worker processes for parallel sampling (default: %(default)s)")
    
    # Parallel jobs argument
    ap.add_argument("--parallel-jobs", type=int, default=1,
                    help="Number of parallel sampling jobs (default: %(default)s, set >1 to run multiple independent sampling runs)")
    
    # Filtering arguments
    ap.add_argument("--enable_filter", action="store_true", default=cfg_def.enable_filter,
                    help="Enable post-generation filtering of seeds")
    ap.add_argument("--contact_cut", type=float, default=cfg_def.contact_cut,
                    help="Contact cutoff distance in Angstrom (default: %(default)s)")
    ap.add_argument("--min_contacts", type=int, default=cfg_def.min_contacts,
                    help="Minimum inter-molecular contacts required (default: %(default)s)")
    ap.add_argument("--max_rg", type=float, default=cfg_def.max_rg,
                    help="Maximum radius of gyration in Angstrom (default: %(default)s)")
    ap.add_argument("--rmsd_dedup", type=float, default=cfg_def.rmsd_dedup,
                    help="RMSD threshold for deduplication (optional)")
    ap.add_argument("--keep_best", type=int, default=cfg_def.keep_best,
                    help="Keep best N seeds by contact count (optional)")
    ap.add_argument("--max_attempts", type=int, default=cfg_def.max_attempts,
                    help="Max generation attempts for nseeds valid seeds (default: %(default)s)")
    ap.add_argument("--core_dist_max", type=float, default=cfg_def.core_dist_max,
                    help="Max PT core-to-core distance in Angstrom (default: %(default)s)")
    ap.add_argument("--head_core_max", type=float, default=cfg_def.head_core_max,
                    help="Max head-to-PT-core distance in Angstrom (default: %(default)s)")
    ap.add_argument("--pt_k", type=int, default=cfg_def.pt_k,
                    help="Number of O atoms for PT core definition (default: %(default)s)")
    
    args = ap.parse_args()

    # Handle --list-monomers
    if args.list_monomers:
        list_available_monomers()
        return

    # Handle parallel jobs mode
    if args.parallel_jobs > 1:
        print("\n" + "="*70)
        print(f"PARALLEL SAMPLING MODE: Running {args.parallel_jobs} independent jobs")
        print("="*70)
        
        base_out = Path(args.out)
        base_out.mkdir(parents=True, exist_ok=True)
        
        # Prepare job arguments
        base_args = sys.argv[1:]  # Get original command line arguments
        
        jobs = []
        for job_id in range(1, args.parallel_jobs + 1):
            folder_name = f"sampling_{job_id:02d}"
            output_dir = str(base_out / folder_name)
            jobs.append((job_id, base_args, output_dir))
        
        # Run jobs sequentially (or in parallel with ThreadPoolExecutor if desired)
        results = []
        print(f"Starting {args.parallel_jobs} sampling jobs with independent random seeds...\n")
        print(f"Base seed: {args.seed}, Job seeds: {args.seed + 100000}, {args.seed + 200000}, ...\n")
        
        # Use ProcessPoolExecutor or run sequentially
        # Sequential approach for stability (each job is a full Python process):
        for job_id, base_args_item, output_dir in jobs:
            success = run_single_parallel_job(job_id, base_args_item, output_dir, args.seed)
            results.append((job_id, success, output_dir))
        
        # Print summary
        print("\n" + "="*70)
        print("PARALLEL JOBS SUMMARY")
        print("="*70)
        successful = sum(1 for _, s, _ in results if s)
        failed = len(results) - successful
        print(f"Successful: {successful}/{args.parallel_jobs}")
        print(f"Failed: {failed}/{args.parallel_jobs}")
        print("\nOutput locations:")
        for job_id, success, output_dir in results:
            status = "✓" if success else "✗"
            print(f"  {status} Job {job_id}: {output_dir}")
        print("="*70)
        
        # Merge seeds_xyz from all jobs if all were successful
        if failed == 0:
            merge_parallel_seeds_xyz(args.out, args.parallel_jobs)
        
        return 0 if failed == 0 else 1

    # Load and display molecules
    print("[MAIN] Loading molecule configuration...")
    try:
        molecules = load_molecules_from_args(args)
        if not molecules:
            print("[ERROR] Failed to load molecules")
            return 1
    except Exception as e:
        print(f"[ERROR] {e}")
        return 1
    
    # Display final molecule composition
    print("\n[MOLECULE COMPOSITION]")
    active_mols = [m for m in molecules if m.enabled]
    disabled_mols = [m for m in molecules if not m.enabled]
    
    if active_mols:
        print("  Enabled:")
        for mol in active_mols:
            print(f"    {mol.name:15s} x {mol.count}  ({mol.file})")
    
    if disabled_mols:
        print("  Disabled:")
        for mol in disabled_mols:
            print(f"    {mol.name:15s} (not sampling)")
    
    # Build sampling engine argument list
    sampling_args = [
        "--out", args.out,
        "--nseeds", str(args.nseeds),
        "--clash_cut", str(args.clash_cut),
        "--max_trials_add", str(args.max_trials_add),
        "--dmin", str(args.dmin),
        "--dmax", str(args.dmax),
        "--lateral", str(args.lateral),
        "--jitter_deg", str(args.jitter_deg),
        "--seed", str(args.seed),
    ]
    
    # Add filter arguments if specified
    if args.enable_filter:
        sampling_args.append("--enable_filter")
    sampling_args.extend([
        "--contact_cut", str(args.contact_cut),
        "--min_contacts", str(args.min_contacts),
        "--max_rg", str(args.max_rg),
        "--max_attempts", str(args.max_attempts),
        "--core_dist_max", str(args.core_dist_max),
        "--head_core_max", str(args.head_core_max),
        "--pt_k", str(args.pt_k),
    ])
    if args.rmsd_dedup is not None:
        sampling_args.extend(["--rmsd_dedup", str(args.rmsd_dedup)])
    if args.keep_best is not None:
        sampling_args.extend(["--keep_best", str(args.keep_best)])
    
    # Build molecules list (only enabled ones)
    import json
    active_molecules = [
        {
            "name": m.name,
            "file": m.file,
            "count": m.count
        }
        for m in active_mols
    ]
    
    if active_molecules:
        sampling_args.extend(["--molecules", json.dumps(active_molecules)])
        print(f"\n[MAIN] Sampling {len(active_mols)} molecule type(s)")
    else:
        print("[ERROR] No molecules enabled for sampling")
        return 1

    # Step 1: Run sampling engine
    print("\n" + "="*70)
    print("STEP 1: Sampling Cluster Geometries")
    print("="*70)
    rc = run_sampling_engine(sampling_args, sampling_workers=args.sampling_workers)
    if rc != 0:
        print(f"[ERROR] Sampling engine exited with code {rc}")
        sys.exit(rc)

    # Step 2: Extract multi-frame XYZ into per-seed files
    print("\n" + "="*70)
    print("STEP 2: Extracting Individual Seed Files")
    print("="*70)
    out_dir = Path(args.out)
    ptn_files = sorted(out_dir.glob("seeds_*.xyz")) + sorted(out_dir.glob("PTN_seeds_N*.xyz"))
    
    if not ptn_files:
        print(f"[WARN] No multi-frame seed files found in {out_dir}")
        print(f"Individual seed files should already be in {out_dir / 'seeds'}")
        return 0

    total = 0
    # Use threads for I/O-bound extraction when workers > 1
    if args.workers and args.workers > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(split_multiframe, p, out_dir / "seeds_xyz"): p for p in ptn_files}
            for fut in concurrent.futures.as_completed(futures):
                p = futures[fut]
                try:
                    n, outdir = fut.result()
                except Exception as e:
                    print(f"[ERROR] splitting {p.name}: {e}")
                    continue
                print(f"[EXTRACT] {p.name} -> {n} files in {outdir}")
                total += n
    else:
        for ptn_file in ptn_files:
            n, outdir = split_multiframe(ptn_file, outdir=out_dir / "seeds_xyz")
            print(f"[EXTRACT] {ptn_file.name} -> {n} files in {outdir}")
            total += n

    print("\n" + "="*70)
    if total > 0:
        print(f"SUCCESS: Extracted {total} additional seed files")
        print(f"Output directory: {(out_dir / 'seeds_xyz').resolve()}")
    else:
        print(f"SUCCESS: Seed files generated in {(out_dir / 'seeds').resolve()}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())

