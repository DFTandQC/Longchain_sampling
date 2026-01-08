# PT Cluster Sampling Pipeline

**A comprehensive multi-monomer cluster generation system using factorial design**. Generates diverse initial structures for subsequent molecular dynamics or quantum chemistry calculations via random sampling and Head-Insert placement strategy.

---

## Project Structure

```
sampling/
├── main.py              ← Main program (user entry point)
├── lib/                 ← Core library folder
│   ├── __init__.py
│   ├── sampling.py      ← Sampling engine (core algorithm)
│   ├── config.py        ← Parameter definitions and CLI
│   └── extractor.py     ← Multi-frame XYZ splitter
├── opt-PT-B97-3c.xyz    ← Monomer structure (input)
├── out/                 ← Output directory
│   ├── seeds/           ← xcontrol constraint files
│   ├── PTN_seeds_N2.xyz ← Multi-frame cluster file
│   └── seeds_xyz/       ← Individual XYZ files (final output)
└── README.md
```

---

## Theory and Principles

### Overall Workflow

```
┌─────────────────────────────────────────────────────────────┐
│ 1. Parameter Setup (CLI or code)                            │
│    - Cluster size N (number of monomers)                   │
│    - Sampling parameters (distance, rotation, clash)       │
└─────────────────────┬───────────────────────────────────────┘
                      ↓
┌─────────────────────────────────────────────────────────────┐
│ 2. Sampling Engine (lib/sampling.py)                       │
│    - Read monomer structure                                │
│    - Build random clusters                                │
│    - Output multi-frame XYZ + xcontrol constraints        │
└─────────────────────┬───────────────────────────────────────┘
                      ↓
┌─────────────────────────────────────────────────────────────┐
│ 3. XYZ Extractor (lib/extractor.py)                       │
│    - Split multi-frame XYZ into single files              │
│    - Naming: seed_0001.xyz, seed_0002.xyz, ...            │
└─────────────────────┬───────────────────────────────────────┘
                      ↓
┌─────────────────────────────────────────────────────────────┐
│ Final Output: out/seeds_xyz/seed_XXXX.xyz                │
│ Ready for: MD, QM, MD optimization, etc.                 │
└─────────────────────────────────────────────────────────────┘
```

### Core Algorithms

#### 1. **Head-Insert Placement Strategy** (lib/sampling.py: `place_head_insert`)

Each monomer has a "head" region (near oxygen atoms). New molecules are positioned as follows:

1. **Define Head Vector**:
   - Head point = centroid of O atoms closest to cluster COM
   - Head vector = unit vector from molecule COM to head point

2. **Orientation**:
   - Select random anchor molecule in existing cluster
   - Rotate new molecule so its head vector points opposite to anchor's head
   - Add small angle jitter for increased diversity

3. **Translation**:
   - COM distance: random within `[dmin, dmax]` range
   - Lateral offset: random perpendicular to head vector (max `lateral`)

#### 2. **Clash Detection** (lib/sampling.py: `clashes_with_cluster`)

Before placement, verify: minimum distance between any atom in new molecule and any atom in cluster is >= `clash_cut` (default 1.2 Å)

#### 3. **Constraint Generation** (lib/sampling.py: `build_OO_pairs_for_cluster`)

For each growth edge (anchor → new):
- Select `kO` nearest O atoms from each molecule's head region
- Randomly pick `pairs_per_edge` O-O pairs between two molecules
- Generate distance constraints for each pair (target distance `OO_target` Å)
- Output in xTB-compatible xcontrol format

#### 4. **Output Format**

**Single XYZ File**(`seed_XXXX.xyz`):
```
178                           # Number of atoms
seed_0001_N2 (N=2)           # Comment line
C    -0.47702853 -0.76711012 0.29100357
C     0.24519277  0.08823437 1.32656798
...
```

**xcontrol Constraints** (`seed_XXXX.xcontrol`):
```
$constrain
  force constant=0.500
  distance: 45, 67, 3.400
  distance: 48, 70, 3.400
...
$end
```

---

## Usage

### Basic Usage

#### Quick Test (5 structures)

```bash
cd e:\Research\Scripts\sampling
python main.py --nseeds 5
```

**Output**: `out/seeds_xyz/seed_0001.xyz` ... `seed_0005.xyz`

#### Generate 100 samples

```bash
python main.py --nseeds 100
```

#### Change cluster size to 3 monomers

```bash
python main.py --nseeds 50 --N 3
```

#### Custom distance and sampling parameters

```bash
python main.py --nseeds 100 \
  --N 3 \
  --dmin 8.0 \
  --dmax 12.0 \
  --clash_cut 1.15 \
  --lateral 2.5
```

---

## Parameter Guide

### Cluster Structure Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--N` | 2 | Number of monomers per cluster |
| `--nseeds` | 200 | Total number of clusters to generate |

### Distance and Layout Parameters

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `--dmin` | 10.0 | Å | Minimum COM distance between new and anchor molecules |
| `--dmax` | 15.0 | Å | Maximum COM distance between new and anchor molecules |
| `--lateral` | 3.0 | Å | Maximum lateral offset value |
| `--jitter_deg` | 25.0 | ° | Maximum additional rotation angle |

### Clash Detection and Placement Parameters

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `--clash_cut` | 1.20 | Å | Clash threshold; reject if atom pair distance < this value |
| `--max_trials_add` | 2000 | - | Maximum attempts to place individual molecule |

### Constraint Generation Parameters

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `--kO` | 4 | - | Number of O atoms selected to define head region |
| `--pairs_per_edge` | 8 | - | Number of O-O distance constraint pairs per growth edge |
| `--OO_target` | 3.40 | Å | Target O-O distance |

### Other Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--mono` | `opt-PT-B97-3c.xyz` | Input monomer file |
| `--out` | `out` | Output directory |
| `--seed` | 7 | Random seed (for reproducibility) |

---

## Output Structure

After execution, the output directory `out/` contains:

```
out/
├── seeds/                    # xcontrol constraint file directory
│   ├── seed_0001_N2.xcontrol
│   ├── seed_0002_N2.xcontrol
│   └── ...
├── PTN_seeds_N2.xyz          # Multi-frame XYZ (all N=2 clusters)
└── seeds_xyz/                # Individual XYZ files (final product)
    ├── seed_0001.xyz
    ├── seed_0002.xyz
    └── ...
```

### Verify Output

Check the number of generated seeds:

```bash
ls out/seeds_xyz/ | wc -l
```

View contents of a single seed:

```bash
head -20 out/seeds_xyz/seed_0001.xyz
```

---

## Advanced Usage

### Integration with MD/QM Software

#### Gromacs MD Preparation

```bash
# Prepare Gromacs topology for each seed
for xyz in out/seeds_xyz/*.xyz; do
    # Use your conversion tool (e.g., Open Babel)
    obabel "$xyz" -O "${xyz%.xyz}.gro"
done
```

#### xTB Geometry Optimization

```bash
# Run xTB optimization using generated xcontrol
for xyz in out/seeds_xyz/*.xyz; do
    seed=$(basename "$xyz" .xyz)
    xtb "$xyz" -opt -c tight --input "out/seeds/${seed}.xcontrol"
done
```

### Parameter Sweep

Sample multiple parameter combinations:

```bash
for dmin in 8 10 12; do
    for N in 2 3 4; do
        python main.py \
            --N $N \
            --dmin $dmin \
            --dmax $((dmin+5)) \
            --nseeds 20 \
            --out "out_N${N}_d${dmin}"
    done
done
```

---

## Troubleshooting

### Issue: Fewer generated seeds than requested nseeds

**Cause**: Sampling failed. Some configurations cannot satisfy clash constraints.

**Solutions**:
- Increase `--dmin` / `--dmax` (larger inter-molecular distances)
- Decrease `--clash_cut` (more lenient clash detection)
- Increase `--max_trials_add` (more placement attempts)
- Decrease `--N` (simpler clusters)

### Issue: Slow sampling

**Cause**: Strict clash detection or attempting to place too many molecules.

**Solutions**:
- Increase `--clash_cut`
- Decrease `--max_trials_add`
- Decrease `--N` or `--nseeds`

### Issue: Import Error `ModuleNotFoundError: No module named 'lib'`

**Cause**: Missing `__init__.py` in `lib` directory (already included).

**Verify**:
```bash
python -c "from lib.config import ClusterConfig; print('OK')"
```

---

## Dependencies

- **Python 3.7+**
- **NumPy** (numerical computations)

Install dependencies:

```bash
pip install numpy
```

---

## License and Citation

If using this tool in scientific publications, please cite:

```
PT Cluster Sampling Pipeline (2026)
Source: e:\Research\Scripts\sampling\
```

---

## Support and Troubleshooting

Having issues? Check:
1. Correct Python version installed (3.7+)
2. NumPy is installed
3. Input monomer file exists
4. Output directory has write permissions
