"""
Configuration module for multi-molecule cluster generation.
Centralizes all parameter definitions and defaults.
Supports mixed clusters with different molecular types.
"""

from dataclasses import dataclass, field
from typing import Dict, List
import argparse
import json


@dataclass
class MoleculeSpec:
    """Specification for a molecular component in the cluster."""
    name: str              # Name of molecule type (e.g., "PT", "H2SO4")
    file: str              # Path to monomer XYZ file
    count: int             # Number of this molecule in each cluster


@dataclass
class ClusterConfig:
    """Parameters for multi-molecule cluster generation."""
    
    # Input/Output
    out_dir: str = "out"
    
    # Molecule composition (new multi-molecule support)
    # Default: single PT molecule with count=2 (backward compatibility)
    molecules: List[MoleculeSpec] = field(default_factory=lambda: [
        MoleculeSpec(name="PT", file="monomer/opt-PT-B97-3c.xyz", count=2)
    ])
    
    # Legacy single-molecule support (for backward compatibility)
    mono_file: str = "monomer/opt-PT-B97-3c.xyz"
    N: int = 2  # Number of monomers per cluster (deprecated, use molecules instead)
    
    # Clash detection
    clash_cut: float = 1.20  # Å
    max_trials_add: int = 2000
    
    # COM distance and positioning
    # For single-component clusters: looser spacing (10-15 Angstrom) is acceptable
    # For mixed clusters: tighter spacing (7-11 Angstrom) prevents natural dispersion
    # Mixed clusters need ~20% tighter packing to maintain structural cohesion
    dmin: float = 10.0  # Angstrom (use 7.0 for mixed clusters via --dmin)
    dmax: float = 15.0  # Angstrom (use 11.0 for mixed clusters via --dmax)
    lateral: float = 3.0  # Angstrom (use 1.5 for mixed clusters via --lateral; tighter lateral displacement)
    jitter_deg: float = 25.0  # degrees
    
    # Restraint generation
    kO: int = 4  # Number of nearest heteroatoms to define head region
    pairs_per_edge: int = 8  # Restraint pairs per attachment edge
    OO_target: float = 3.40  # Angstrom, target distance between constraint atoms
    k_fc: float = 0.50  # Force constant for restraints
    
    # Seed generation
    nseeds: int = 200
    random_seed: int = 7
    
    # Post-generation filtering
    enable_filter: bool = False
    contact_cut: float = 3.8  # Angstrom
    min_contacts: int = 20
    max_rg: float = 25.0  # Angstrom
    rmsd_dedup: float = None  # Optional: if None, skip dedup
    keep_best: int = None  # Optional: limit final number saved
    max_attempts: int = 10000  # Max attempts to generate nseeds valid seeds
    core_dist_max: float = 8.0  # Angstrom, max PT core-to-core distance (if PT exists)
    head_core_max: float = 6.0  # Angstrom, max non-PT head-to-PT-core distance (if PT exists)
    pt_k: int = 8  # Number of O atoms for PT core definition
    
    def to_dict(self):
        """Convert to dictionary for easy access."""
        return self.__dict__
    
    @classmethod
    def from_args(cls, args):
        """Create config from argparse Namespace."""
        return cls(
            mono_file=args.mono,
            out_dir=args.out,
            N=args.N,
            clash_cut=args.clash_cut,
            max_trials_add=args.max_trials_add,
            dmin=args.dmin,
            dmax=args.dmax,
            lateral=args.lateral,
            jitter_deg=args.jitter_deg,
            kO=args.kO,
            pairs_per_edge=args.pairs_per_edge,
            OO_target=args.OO_target,
            k_fc=args.k_fc,
            nseeds=args.nseeds,
            random_seed=args.seed,
        )


def get_argument_parser():
    """Create and return argument parser with all cluster parameters.

    Single, clean implementation (no duplicates) that supports legacy single-molecule
    mode and multi-molecule JSON configuration. Guidance in help strings suggests
    tighter `--dmin/--dmax/--lateral` values for compact mixed-cluster sampling.
    """
    cfg_def = ClusterConfig()
    ap = argparse.ArgumentParser(
        description="Generate multi-molecule cluster seeds with Head-insert biased growth."
    )

    # I/O
    ap.add_argument("--mono", default=cfg_def.mono_file,
                    help="[Legacy] Monomer xyz file (use --config for multi-molecule)")
    ap.add_argument("--out", default=cfg_def.out_dir,
                    help="Output directory")
    ap.add_argument("--nseeds", type=int, default=cfg_def.nseeds,
                    help="Number of cluster seeds to generate")

    # Multi-molecule configuration
    ap.add_argument("--config", type=str, default=None,
                    help="JSON config file specifying molecule composition (overrides --mono and --N)")
    ap.add_argument("--molecules", type=str, default=None,
                    help='Molecule spec as JSON string e.g. "[{\"name\":\"PT\",\"file\":\"opt-PT.xyz\",\"count\":2}]"')

    # Cluster size (legacy for single-molecule mode)
    ap.add_argument("--N", type=int, default=cfg_def.N,
                    help="[Legacy] Number of monomers per cluster (use --molecules for multi-molecule)")

    # Clash and placement
    ap.add_argument("--clash_cut", type=float, default=cfg_def.clash_cut,
                    help="Reject if any inter-molecular atom distance < clash_cut (Angstrom)")
    ap.add_argument("--max_trials_add", type=int, default=cfg_def.max_trials_add,
                    help="Max trials to place each new monomer")
    ap.add_argument("--dmin", type=float, default=cfg_def.dmin,
                    help="Min COM distance Angstrom (default 10.0; consider 7.0 for compact mixed clusters)")
    ap.add_argument("--dmax", type=float, default=cfg_def.dmax,
                    help="Max COM distance Angstrom (default 15.0; consider 11.0 for compact mixed clusters)")
    ap.add_argument("--lateral", type=float, default=cfg_def.lateral,
                    help="Max lateral offset Angstrom (default 3.0; consider 1.5 for compact mixed clusters)")
    ap.add_argument("--jitter_deg", type=float, default=cfg_def.jitter_deg,
                    help="Max angle jitter (deg)")

    # Restraints
    ap.add_argument("--kO", type=int, default=cfg_def.kO,
                    help="Number of nearest heteroatoms to define head-region for restraints (per monomer)")
    ap.add_argument("--pairs_per_edge", type=int, default=cfg_def.pairs_per_edge,
                    help="Number of restraint pairs per attachment edge (heteroatom-based)")
    ap.add_argument("--OO_target", type=float, default=cfg_def.OO_target,
                    help="Target distance (Angstrom) for restraints (atom-pair based)")
    ap.add_argument("--k_fc", type=float, default=cfg_def.k_fc,
                    help="Force constant for distance restraints")

    # Random seed
    ap.add_argument("--seed", type=int, default=cfg_def.random_seed,
                    help="Random seed for reproducibility")

    # Post-generation filtering
    ap.add_argument("--enable_filter", action="store_true", default=cfg_def.enable_filter,
                    help="Enable post-generation filtering of seeds")
    ap.add_argument("--contact_cut", type=float, default=cfg_def.contact_cut,
                    help="Contact cutoff distance in Angstrom (default 3.8)")
    ap.add_argument("--min_contacts", type=int, default=cfg_def.min_contacts,
                    help="Minimum number of inter-molecular contacts required (default 20)")
    ap.add_argument("--max_rg", type=float, default=cfg_def.max_rg,
                    help="Maximum allowed radius of gyration in Angstrom (default 25.0)")
    ap.add_argument("--rmsd_dedup", type=float, default=cfg_def.rmsd_dedup,
                    help="RMSD threshold for duplicate detection; if not set, dedup disabled (optional)")
    ap.add_argument("--keep_best", type=int, default=cfg_def.keep_best,
                    help="Limit output to best N seeds by contact count (optional)")
    ap.add_argument("--max_attempts", type=int, default=cfg_def.max_attempts,
                    help="Maximum generation attempts to produce nseeds valid seeds (default 10000)")
    ap.add_argument("--core_dist_max", type=float, default=cfg_def.core_dist_max,
                    help="Max PT core-to-core distance for mixed PT systems in Angstrom (default 8.0)")
    ap.add_argument("--head_core_max", type=float, default=cfg_def.head_core_max,
                    help="Max non-PT head-to-PT-core distance for mixed PT systems in Angstrom (default 6.0)")
    ap.add_argument("--pt_k", type=int, default=cfg_def.pt_k,
                    help="Number of O atoms for PT core definition (default 8)")

    return ap



def load_molecules_from_args(args) -> List[MoleculeSpec]:
    """Parse molecule configuration from command-line arguments."""
    # Priority: --config > --molecules > legacy (--mono + --N)
    
    if hasattr(args, 'config') and args.config:
        # Load from JSON config file
        with open(args.config, 'r') as f:
            cfg = json.load(f)
        molecules = [
            MoleculeSpec(name=m['name'], file=m['file'], count=m['count'])
            for m in cfg.get('molecules', [])
        ]
        return molecules if molecules else None
    
    if hasattr(args, 'molecules') and args.molecules:
        # Parse from JSON string
        mol_list = json.loads(args.molecules)
        molecules = [
            MoleculeSpec(name=m['name'], file=m['file'], count=m['count'])
            for m in mol_list
        ]
        return molecules if molecules else None
    
    # Legacy mode: single molecule with count N
    mono = getattr(args, 'mono', 'opt-PT-B97-3c.xyz')
    N = getattr(args, 'N', 2)
    return [MoleculeSpec(name="PT", file=mono, count=N)]


def create_config_from_args(args) -> 'ClusterConfig':
    """Create config from argparse Namespace, supporting both legacy and multi-molecule modes."""
    molecules = load_molecules_from_args(args)
    
    cfg = ClusterConfig(
        out_dir=args.out,
        molecules=molecules,
        mono_file=getattr(args, 'mono', 'opt-PT-B97-3c.xyz'),
        N=getattr(args, 'N', 2),
        clash_cut=args.clash_cut,
        max_trials_add=args.max_trials_add,
        dmin=args.dmin,
        dmax=args.dmax,
        lateral=args.lateral,
        jitter_deg=args.jitter_deg,
        kO=args.kO,
        pairs_per_edge=args.pairs_per_edge,
        OO_target=args.OO_target,
        k_fc=args.k_fc,
        nseeds=args.nseeds,
        random_seed=args.seed,
        enable_filter=args.enable_filter,
        contact_cut=args.contact_cut,
        min_contacts=args.min_contacts,
        max_rg=args.max_rg,
        rmsd_dedup=args.rmsd_dedup,
        keep_best=args.keep_best,
        max_attempts=args.max_attempts,
        core_dist_max=args.core_dist_max,
        head_core_max=args.head_core_max,
        pt_k=args.pt_k,
    )
    return cfg