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
"""
from pathlib import Path
import subprocess
import sys
import glob
import argparse

from lib.config import ClusterConfig, get_argument_parser
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
  python main.py --config cluster_config.json --nseeds 100
  python main.py --molecules '[{"name":"PT","file":"opt-PT.xyz","count":2},{"name":"H2SO4","file":"h2so4.xyz","count":1}]' --nseeds 100
        """
    )
    
    # Reuse config defaults
    cfg_def = ClusterConfig()
    ap.add_argument("--config", type=str, default="config.json",
                    help="JSON config file for multi-molecule clusters (default: config.json with presets)")
    ap.add_argument("--molecules", type=str, default=None,
                    help="Molecule spec as JSON string (overrides --mono and --N)")
    ap.add_argument("--mono", default=cfg_def.mono_file,
                    help="[Legacy] Monomer xyz file (default: %(default)s)")
    ap.add_argument("--out", default=cfg_def.out_dir,
                    help="Output directory (default: %(default)s)")
    ap.add_argument("--nseeds", type=int, default=cfg_def.nseeds,
                    help="Number of cluster seeds to generate (default: %(default)s)")
    ap.add_argument("--N", type=int, default=cfg_def.N,
                    help="[Legacy] Monomers per cluster (default: %(default)s)")
    ap.add_argument("--clash_cut", type=float, default=cfg_def.clash_cut,
                    help="Clash detection threshold Å (default: %(default)s)")
    ap.add_argument("--max_trials_add", type=int, default=cfg_def.max_trials_add,
                    help="Max placement attempts per monomer (default: %(default)s)")
    ap.add_argument("--dmin", type=float, default=cfg_def.dmin,
                    help="Min COM distance Å (default: %(default)s)")
    ap.add_argument("--dmax", type=float, default=cfg_def.dmax,
                    help="Max COM distance Å (default: %(default)s)")
    ap.add_argument("--lateral", type=float, default=cfg_def.lateral,
                    help="Max lateral offset Å (default: %(default)s)")
    ap.add_argument("--jitter_deg", type=float, default=cfg_def.jitter_deg,
                    help="Angle jitter degrees (default: %(default)s)")
    ap.add_argument("--seed", type=int, default=cfg_def.random_seed,
                    help="Random seed for reproducibility (default: %(default)s)")
    ap.add_argument("--workers", type=int, default=1,
                    help="Worker threads for extraction (default: %(default)s)")
    ap.add_argument("--sampling-workers", type=int, default=1,
                    help="Worker processes for parallel sampling (default: %(default)s)")
    args = ap.parse_args()

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
    
    # Handle multi-molecule vs legacy mode
    if args.config:
        sampling_args.extend(["--config", args.config])
        print(f"[MAIN] Using multi-molecule config: {args.config}")
    elif args.molecules:
        sampling_args.extend(["--molecules", args.molecules])
        print(f"[MAIN] Using molecules from command line")
    else:
        # Legacy mode
        sampling_args.extend([
            "--mono", args.mono,
            "--N", str(args.N)
        ])
        print(f"[MAIN] Using legacy single-molecule mode: {args.mono} (N={args.N})")

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
                print(f"[EXTRACT] {p.name} → {n} files in {outdir}")
                total += n
    else:
        for ptn_file in ptn_files:
            n, outdir = split_multiframe(ptn_file, outdir=out_dir / "seeds_xyz")
            print(f"[EXTRACT] {ptn_file.name} → {n} files in {outdir}")
            total += n

    print("\n" + "="*70)
    if total > 0:
        print(f"SUCCESS: Extracted {total} additional seed files")
        print(f"Output directory: {(out_dir / 'seeds_xyz').resolve()}")
    else:
        print(f"SUCCESS: Seed files generated in {(out_dir / 'seeds').resolve()}")

if __name__ == "__main__":
    main()
