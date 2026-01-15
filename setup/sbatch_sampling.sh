#!/bin/bash
#SBATCH --job-name=xtb_opt
#SBATCH --account=hvehkama
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4000
#SBATCH --array=1-200%40
#SBATCH --output=logs.out
#SBATCH --error=logs.err

set -euo pipefail

# ====== Directories ======
SEED_DIR="seeds_xyz"
OUT_DIR="opt_xtb"
mkdir -p logs "$OUT_DIR"

# ====== Modules ======
module purge
module load xtb

# ====== Thread settings ======
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ====== Gather file list and pick one by array index ======
mapfile -t files < <(ls -1 "${SEED_DIR}"/*.xyz | sort)

idx=$((SLURM_ARRAY_TASK_ID-1))
if [[ $idx -lt 0 || $idx -ge ${#files[@]} ]]; then
  echo "Array index out of range: SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID, nfiles=${#files[@]}"
  exit 2
fi

f="${files[$idx]}"
base=$(basename "$f" .xyz)

workdir="${OUT_DIR}/${base}"
mkdir -p "$workdir"

# ====== Run in scratch (faster; avoids concurrent writes to home) ======
rundir="${TMPDIR:-$workdir}/${base}"
mkdir -p "$rundir"
cp "$f" "$rundir/seed.xyz"
cd "$rundir"

echo "Input: $f"
echo "Run dir: $(pwd)"
echo "Threads: $SLURM_CPUS_PER_TASK"

# ====== xTB quick optimization ======
# You may change `--opt tight` to `--opt vtight` for a stricter (and slower) optimization
srun xtb seed.xyz --gfn 1 --opt tight > xtb.out

# ====== Copy results back to permanent directory ======
cp -f xtb.out "$workdir/"
cp -f seed.xyz "$workdir/${base}_seed.xyz"

if [[ -f xtbopt.xyz ]]; then
  cp -f xtbopt.xyz "$workdir/${base}_xtbopt.xyz"
fi

# Some xTB builds also produce xtbopt.log / xtbrestart; copy them if present
for extra in xtbopt.log xtbrestart; do
  [[ -f "$extra" ]] && cp -f "$extra" "$workdir/"
done

echo "Done: $base"