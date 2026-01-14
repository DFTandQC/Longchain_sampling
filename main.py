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
"""
from pathlib import Path
import subprocess
import sys
import glob
import argparse

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
        return

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

if __name__ == "__main__":
    main()

