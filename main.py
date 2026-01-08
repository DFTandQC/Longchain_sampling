#!/usr/bin/env python3
"""
PT Cluster Sampling Pipeline - Main Entry Point

Run the full sampling workflow: generate multi-frame clusters with sampling engine,
then split output frames into individual seed XYZ files.

Usage:
    python main.py --nseeds 100
    python main.py --nseeds 50 --N 3 --dmin 8.0 --dmax 12.0
"""
from pathlib import Path
import subprocess
import sys
import glob
import argparse

from lib.config import ClusterConfig, get_argument_parser
from lib.extractor import split_multiframe


def run_sampling_engine(args_list):
    """Invoke the sampling engine (lib/sampling.py) as a module."""
    cmd = [sys.executable, "-m", "lib.sampling"] + args_list
    print("[MAIN] Running sampling engine:", " ".join(cmd))
    res = subprocess.run(cmd)
    return res.returncode


def main():
    ap = argparse.ArgumentParser(
        description="Generate PT cluster seeds: sample geometries then extract into individual files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --nseeds 100              # Generate 100 seeds (default N=2)
  python main.py --nseeds 50 --N 3         # Generate 50 trimers
  python main.py --nseeds 20 --dmin 8 --dmax 12  # Custom distance range
        """
    )
    
    # Reuse config defaults
    cfg_def = ClusterConfig()
    ap.add_argument("--mono", default=cfg_def.mono_file,
                    help="Monomer xyz file (default: %(default)s)")
    ap.add_argument("--out", default=cfg_def.out_dir,
                    help="Output directory (default: %(default)s)")
    ap.add_argument("--nseeds", type=int, default=cfg_def.nseeds,
                    help="Number of cluster seeds to generate (default: %(default)s)")
    ap.add_argument("--N", type=int, default=cfg_def.N,
                    help="Monomers per cluster (default: %(default)s)")
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
    args = ap.parse_args()

    # Build sampling engine argument list
    sampling_args = [
        "--mono", args.mono,
        "--out", args.out,
        "--nseeds", str(args.nseeds),
        "--N", str(args.N),
        "--clash_cut", str(args.clash_cut),
        "--max_trials_add", str(args.max_trials_add),
        "--dmin", str(args.dmin),
        "--dmax", str(args.dmax),
        "--lateral", str(args.lateral),
        "--jitter_deg", str(args.jitter_deg),
        "--seed", str(args.seed),
    ]

    # Step 1: Run sampling engine
    print("\n" + "="*70)
    print("STEP 1: Sampling Cluster Geometries")
    print("="*70)
    rc = run_sampling_engine(sampling_args)
    if rc != 0:
        print(f"[ERROR] Sampling engine exited with code {rc}")
        sys.exit(rc)

    # Step 2: Extract multi-frame XYZ into per-seed files
    print("\n" + "="*70)
    print("STEP 2: Extracting Individual Seed Files")
    print("="*70)
    out_dir = Path(args.out)
    ptn_files = sorted(out_dir.glob("PTN_seeds_N*.xyz"))
    
    if not ptn_files:
        print(f"[ERROR] No PTN_seeds_N*.xyz files found in {out_dir}")
        sys.exit(1)

    total = 0
    for ptn_file in ptn_files:
        n, outdir = split_multiframe(ptn_file, outdir=out_dir / "seeds_xyz")
        print(f"[EXTRACT] {ptn_file.name} → {n} files in {outdir}")
        total += n

    print("\n" + "="*70)
    print(f"SUCCESS: Generated {total} seed files")
    print(f"Output directory: {(out_dir / 'seeds_xyz').resolve()}")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
