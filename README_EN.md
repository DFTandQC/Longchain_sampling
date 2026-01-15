# Multi-Molecule Cluster Sampling System

> Multi-molecule mixed cluster initial structure generation engine

## Quick Start

```bash
# List available molecules
python main.py --list-monomers

# Select molecules and generate
python main.py --use "PT,H2SO4,NO2" --nseeds 100

# Custom counts
python main.py --use "PT,H2SO4,NO,NO2" \
  --counts "PT=2,H2SO4=1,NO=2,NO2=1" --nseeds 200
```

## Core Features

- ✅ **Head-Insert Placement Strategy** - Incrementally build multi-molecule clusters
- ✅ **PT Ester-core Alignment** - Three-branch strategy (pure PT/mixed/single PT)
- ✅ **Post-processing Filtering** - Rg/contacts/dedup/PT constraints
- ✅ **Smart Monomer Discovery** - Automatic file matching
- ✅ **Flexible Molecule Selection** - CLI selection and count override
- ✅ **Parallel Sampling** - Multi-process acceleration

## Documentation

📚 **Full Documentation**: [docs/FULL_DOCUMENTATION_EN.md](docs/FULL_DOCUMENTATION_EN.md) (English)
- System architecture and data flow
- Core algorithm details
- PT special strategy
- Post-processing filtering mechanisms
- Parameter tuning guide
- Usage examples
- Troubleshooting

📚 **Chinese Documentation**: [docs/FULL_DOCUMENTATION.md](docs/FULL_DOCUMENTATION.md)
- Complete Chinese user guide

## Project Structure

```
sampling/
├── setup_environment.py         # Auto-configuration script 🆕
├── requirements.txt             # Dependency list 🆕
├── main.py                      # Entry point
├── config.json                  # Main config (6 molecules)
├── lib/
│   ├── config.py               # Config loading
│   ├── sampling.py             # Sampling algorithm
│   └── extractor.py            # File extraction
├── monomer/                     # Monomer library (6 types)
│   ├── opt-PT-B97-3c.xyz
│   ├── opt-cisSA-B97-3c.xyz   # H2SO4
│   ├── opt-transSA-B97-3c.xyz
│   ├── opt-NO-B97-3c.xyz
│   ├── opt-NO2-B97-3c.xyz
│   └── opt-SO2-B97-3c.xyz
├── test/                        # Test files
│   ├── test_pt_core.py         # PT core tests
│   ├── test_filtering.py       # Filtering tests
│   ├── tests.py                # Main test suite
│   └── smoke_test.py           # Quick verification
├── out/                         # Output (auto-created)
│   ├── seeds/                  # Individual seeds
│   ├── seeds_xyz/              # Extracted
│   └── seeds_PT_combined.xyz  # Multi-frame combined
└── docs/
    ├── FULL_DOCUMENTATION_EN.md # Full documentation
    └── FULL_DOCUMENTATION.md    # Chinese version
```

## Common Commands

### Basic Usage

```bash
# List available molecules
python main.py --list-monomers

# Use default configuration
python main.py --nseeds 100

# Select specific molecules
python main.py --use "PT,H2SO4" --nseeds 50
```

### Enable Filtering

```bash
python main.py --use "PT,H2SO4,NO,NO2" \
  --counts "PT=2,H2SO4=1,NO=2,NO2=1" \
  --enable_filter \
  --min_contacts 20 \
  --max_rg 25.0 \
  --nseeds 200
```

### Fine-tune Parameters

```bash
# Tight clustering
python main.py --use "PT,H2SO4" \
  --dmin 7.0 --dmax 10.0 --lateral 1.5 \
  --nseeds 100

# Loose clustering
python main.py --use "PT,H2SO4" \
  --dmin 12.0 --dmax 18.0 --lateral 4.0 \
  --nseeds 100
```

### Parallel Acceleration

#### Parallel within Single Sampling Task

```bash
# Use 4 processes to accelerate sampling
python main.py --use "PT,H2SO4,NO2" \
  --sampling-workers 4 --nseeds 500
```

#### Multiple Independent Sampling Tasks (Recommended)

```bash
# Run 5 independent sampling tasks, each generating 100 seeds
python main.py --nseeds 100 --parallel-jobs 5

# Using config file, run 3 independent tasks
python main.py --config config.json --nseeds 100 --parallel-jobs 3

# Run 10 tasks, sampling specific molecule combinations
python main.py --use "PT,H2SO4,NO2" --nseeds 50 --parallel-jobs 10
```

**Output Structure**:
```
out/
├── sampling_01/  # Task 1 output
├── sampling_02/  # Task 2 output
├── sampling_03/  # Task 3 output
└── ...
```

**Advantages**:
- ✅ Completely independent sampling processes
- ✅ Auto-creates separated output folders
- ✅ Great for generating more diverse initial structures
- ✅ Outputs can be directly merged

## Key Parameters

| Parameter | Recommended (Pure PT) | Recommended (Mixed) | Description |
|-----------|----------------------|---------------------|-------------|
| `dmin` | 10.0 | 7.0-8.0 | COM distance minimum (Å) |
| `dmax` | 15.0 | 10.0-11.0 | COM distance maximum (Å) |
| `lateral` | 3.0 | 1.5-2.0 | Lateral displacement (Å) |
| `core_dist_max` | — | 6.0-8.0 | PT core-core distance (Å) |
| `head_core_max` | — | 5.0-6.0 | Non-PT head to PT core (Å) |

## Output Files

```
out/
├── seeds/
│   ├── seed_10001_PT2_H2SO41_NO21.xyz      # Geometry
│   ├── seed_10001_PT2_H2SO41_NO21.xcontrol # Constraints
│   └── ...
├── seeds_PT_combined.xyz                    # Multi-frame combined
└── seeds_xyz/                               # Extracted (optional)
```

## Adding New Molecules

1. Optimize molecular structure (xTB/ORCA/Gaussian)
2. Place at `monomer/opt-MOLECULE.xyz`
3. Update `config.json` or use `--use`
4. Test: `python main.py --use "MOLECULE,PT" --nseeds 10`

Supports molecules with O/N/S atoms or pure hydrocarbons.

## Troubleshooting

### Generating 0 Seeds

```bash
# Relax constraints
python main.py --use "PT,NO2" \
  --core_dist_max 10.0 --dmin 7.0 --dmax 11.0 --nseeds 100
```

### Molecules Too Dispersed

```bash
# Tighten distances
python main.py --use "PT,H2SO4" \
  --dmin 6.0 --dmax 9.0 --lateral 1.5 \
  --enable_filter --min_contacts 25 --nseeds 100
```

### High Filter Rejection Rate

```bash
# Relax filter conditions
python main.py --use "PT,H2SO4" \
  --enable_filter --min_contacts 15 --max_rg 30.0 --nseeds 100
```

## Quick Configuration (New Environment)

On first run on a new computer:

```bash
# Auto-configure environment (recommended)
python setup_environment.py
```

Auto-executes:
- ✓ Check Python version (3.8+)
- ✓ Install dependencies (NumPy)
- ✓ Verify project structure
- ✓ Run test verification
- ✓ Show quick start guide

Or install manually:

```bash
pip install -r requirements.txt
```

**Dependencies**:
- Python 3.8+
- NumPy 1.20+

## Version

**v2.0** (2026-01-15)
- PT ester-core placement strategy
- Post-processing filtering system
- Monomer discovery and molecule selection
- Complete documentation

---

📚 See [docs/FULL_DOCUMENTATION_EN.md](docs/FULL_DOCUMENTATION_EN.md) for complete user guide
