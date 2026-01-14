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
    # For single-component clusters: looser spacing (10–15 Å) is acceptable
    # For mixed clusters: tighter spacing (7–11 Å) prevents natural dispersion
    # Mixed clusters need ~20% tighter packing to maintain structural cohesion
    dmin: float = 10.0  # Å (use 7.0 for mixed clusters via --dmin)
    dmax: float = 15.0  # Å (use 11.0 for mixed clusters via --dmax)
    lateral: float = 3.0  # Å (use 1.5 for mixed clusters via --lateral; tighter lateral displacement)
    jitter_deg: float = 25.0  # degrees
    
    # Restraint generation
    kO: int = 4  # Number of nearest heteroatoms to define head region
    pairs_per_edge: int = 8  # Restraint pairs per attachment edge
    OO_target: float = 3.40  # Å, target distance between constraint atoms
    k_fc: float = 0.50  # Force constant for restraints
    
    # Seed generation
    nseeds: int = 200
    random_seed: int = 7
    
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
                    help="Reject if any inter-molecular atom distance < clash_cut (Å)")
    ap.add_argument("--max_trials_add", type=int, default=cfg_def.max_trials_add,
                    help="Max trials to place each new monomer")
    ap.add_argument("--dmin", type=float, default=cfg_def.dmin,
                    help="Min COM distance Å (default 10.0; consider 7.0 for compact mixed clusters)")
    ap.add_argument("--dmax", type=float, default=cfg_def.dmax,
                    help="Max COM distance Å (default 15.0; consider 11.0 for compact mixed clusters)")
    ap.add_argument("--lateral", type=float, default=cfg_def.lateral,
                    help="Max lateral offset Å (default 3.0; consider 1.5 for compact mixed clusters)")
    ap.add_argument("--jitter_deg", type=float, default=cfg_def.jitter_deg,
                    help="Max angle jitter (deg)")

    # Restraints
    ap.add_argument("--kO", type=int, default=cfg_def.kO,
                    help="Number of nearest heteroatoms to define head-region for restraints (per monomer)")
    ap.add_argument("--pairs_per_edge", type=int, default=cfg_def.pairs_per_edge,
                    help="Number of restraint pairs per attachment edge (heteroatom-based)")
    ap.add_argument("--OO_target", type=float, default=cfg_def.OO_target,
                    help="Target distance (Å) for restraints (atom-pair based)")
    ap.add_argument("--k_fc", type=float, default=cfg_def.k_fc,
                    help="Force constant for distance restraints")

    # Random seed
    ap.add_argument("--seed", type=int, default=cfg_def.random_seed,
                    help="Random seed for reproducibility")

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
    )
    return cfg