from pathlib import Path
import re


def split_multiframe(inpath, outdir=None):
    """Split multi-frame XYZ file into per-seed files.
    
    Args:
        inpath: Path to multi-frame XYZ file
        outdir: Output directory (default: seeds_xyz next to inpath)
    
    Returns:
        (num_files, output_directory)
    """
    inp = Path(inpath)
    if outdir is None:
        outdir = inp.parent / "seeds_xyz"
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    lines = inp.read_text().splitlines()
    i = 0
    k = 0

    while i < len(lines):
        n = int(lines[i].strip())
        comment = lines[i+1].strip()

        block = lines[i:i+2+n]
        k += 1

        # Prefer to extract seed name from comment if available
        m = re.search(r"(seed_\d+)", comment)
        name = m.group(1) if m else f"seed_{k:04d}"

        (outdir / f"{name}.xyz").write_text("\n".join(block) + "\n")
        i += 2 + n

    return k, outdir


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Split multi-frame XYZ into per-seed XYZ files.")
    ap.add_argument("infile", help="Input multi-frame XYZ (e.g. out/PTN_seeds_N2.xyz)")
    ap.add_argument("--out", help="Output directory (default: seeds_xyz next to infile)")
    args = ap.parse_args()
    n, outdir = split_multiframe(args.infile, args.out)
    print(f"Wrote {n} structures to {outdir}")
