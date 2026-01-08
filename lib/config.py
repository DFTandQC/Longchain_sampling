"""
Configuration module for PT cluster generation.
Centralizes all parameter definitions and defaults.
"""

from dataclasses import dataclass
import argparse


@dataclass
class ClusterConfig:
    """Parameters for PT cluster generation."""
    
    # Input/Output
    mono_file: str = "opt-PT-B97-3c.xyz"
    out_dir: str = "out"
    
    # Cluster size
    N: int = 2 # Number of monomers per cluster
    
    # Clash detection
    clash_cut: float = 1.20  # Å
    max_trials_add: int = 2000
    
    # COM distance and positioning
    dmin: float = 10.0  # Å
    dmax: float = 15.0  # Å
    lateral: float = 3.0  # Å (max lateral offset)
    jitter_deg: float = 25.0  # degrees
    
    # Restraint generation
    kO: int = 4  # Number of nearest O atoms to define head region
    pairs_per_edge: int = 8  # O-O restraints per attachment edge
    OO_target: float = 3.40  # Å, target O-O distance
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
    """Create and return argument parser with all cluster parameters."""
    cfg_def = ClusterConfig()
    ap = argparse.ArgumentParser(
        description="Generate PT_N cluster seeds with random N and Head-insert biased growth."
    )
    
    # I/O
    ap.add_argument("--mono", default=cfg_def.mono_file, 
                    help="Monomer xyz file")
    ap.add_argument("--out", default=cfg_def.out_dir, 
                    help="Output directory")
    ap.add_argument("--nseeds", type=int, default=cfg_def.nseeds, 
                    help="Number of cluster seeds to generate")
    
    # Cluster size
    ap.add_argument("--N", type=int, default=cfg_def.N, 
                    help="Number of monomers per cluster")
    
    # Clash and placement
    ap.add_argument("--clash_cut", type=float, default=cfg_def.clash_cut, 
                    help="Reject if any inter-molecular atom distance < clash_cut (Å)")
    ap.add_argument("--max_trials_add", type=int, default=cfg_def.max_trials_add, 
                    help="Max trials to place each new monomer")
    ap.add_argument("--dmin", type=float, default=cfg_def.dmin, 
                    help="Min COM distance when attaching new monomer (Å)")
    ap.add_argument("--dmax", type=float, default=cfg_def.dmax, 
                    help="Max COM distance when attaching new monomer (Å)")
    ap.add_argument("--lateral", type=float, default=cfg_def.lateral, 
                    help="Max lateral offset (Å)")
    ap.add_argument("--jitter_deg", type=float, default=cfg_def.jitter_deg, 
                    help="Max angle jitter (deg)")
    
    # Restraints
    ap.add_argument("--kO", type=int, default=cfg_def.kO, 
                    help="Number of nearest O atoms to define head-region for restraints (per monomer)")
    ap.add_argument("--pairs_per_edge", type=int, default=cfg_def.pairs_per_edge, 
                    help="Number of O-O restraints per attachment edge")
    ap.add_argument("--OO_target", type=float, default=cfg_def.OO_target, 
                    help="Target O-O distance (Å) for restraints")
    ap.add_argument("--k_fc", type=float, default=cfg_def.k_fc, 
                    help="Force constant for distance restraints")
    
    # Random seed
    ap.add_argument("--seed", type=int, default=cfg_def.random_seed, 
                    help="Random seed for reproducibility")
    
    return ap
