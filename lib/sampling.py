#!/usr/bin/env python3
from pathlib import Path
import numpy as np

from .config import ClusterConfig, MoleculeSpec, get_argument_parser, create_config_from_args


# ---------------- I/O ----------------
def read_xyz(path: str):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].strip())
    comment = lines[1].rstrip("\n")
    elems, xyz = [], []
    for ln in lines[2:2+n]:
        s = ln.split()
        elems.append(s[0])
        xyz.append([float(s[1]), float(s[2]), float(s[3])])
    return np.array(elems, dtype=object), np.array(xyz, dtype=float), comment

def write_xyz(path: str, elems, X, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(elems)}\n")
        f.write(comment + "\n")
        for e, (x,y,z) in zip(elems, X):
            f.write(f"{e:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")

def write_xyz_multiframe(path: str, frames, elems, comments):
    with open(path, "w") as f:
        for X, cmt in zip(frames, comments):
            f.write(f"{len(elems)}\n{cmt}\n")
            for e, (x,y,z) in zip(elems, X):
                f.write(f"{e:2s} {x:16.8f} {y:16.8f} {z:16.8f}\n")

# ---------------- geometry ----------------
def com(X):
    return X.mean(axis=0)

def unit(v, eps=1e-12):
    n = np.linalg.norm(v)
    return v*0.0 if n < eps else v/n

def skew(v):
    x,y,z = v
    return np.array([[0,-z, y],
                     [z, 0,-x],
                     [-y,x, 0]], float)

def axis_angle(axis, theta):
    axis = unit(axis)
    K = skew(axis)
    return np.eye(3) + np.sin(theta)*K + (1-np.cos(theta))*(K@K)

def rodrigues_align(a, b):
    """Return rotation matrix R such that R@a ~ b (a,b are 3-vectors)."""
    a = unit(a); b = unit(b)
    v = np.cross(a, b)
    s = np.linalg.norm(v)
    c = float(np.dot(a, b))
    if s < 1e-12:
        if c > 0:
            return np.eye(3)
        # 180° rotation around any axis perpendicular to a
        tmp = np.array([1.0,0.0,0.0])
        if abs(np.dot(tmp, a)) > 0.9:
            tmp = np.array([0.0,1.0,0.0])
        axis = unit(np.cross(a, tmp))
        return axis_angle(axis, np.pi)
    vx = skew(v)
    return np.eye(3) + vx + vx@vx*((1-c)/(s**2))

def random_unit(rng):
    v = rng.normal(size=3)
    return unit(v)

def random_twist_about(axis, rng):
    return axis_angle(axis, rng.uniform(0, 2*np.pi))

# ---------------- Head definition (flexible, supports multiple element types) ----------------
def get_head_atoms_priority(elems):
    """
    Determine which atoms define the molecular 'head' region.
    
    Priority order:
    1. Oxygen atoms (O) - standard for acids, etc.
    2. Nitrogen atoms (N) - for amines, ammonium
    3. Sulfur atoms (S) - for thiols, sulfides
    4. All atoms - fallback if no heteroatoms found
    
    This allows robust head definition for diverse molecules:
    - H2SO4: uses O atoms
    - NH3: uses N atom
    - H2O: uses O atom
    - Thiols: uses S atom
    """
    atom_types = np.unique(elems)
    
    # Try priority order
    for target in ["O", "N", "S"]:
        if target in atom_types:
            return np.where(elems == target)[0]
    
    # Fallback: use all atoms
    return np.arange(len(elems))

def head_point_from_atoms(elems, X, atom_indices=None):
    """
    Calculate head point as centroid of specified atoms (or heteroatoms if not specified).
    
    Args:
        elems: array of element symbols
        X: coordinate array
        atom_indices: specific atom indices to use; if None, use get_head_atoms_priority()
    
    Returns:
        3D head point coordinates
    """
    if atom_indices is None:
        atom_indices = get_head_atoms_priority(elems)
    
    if len(atom_indices) == 0:
        raise ValueError("No suitable atoms found to define head point")
    
    return X[atom_indices].mean(axis=0)

def head_point_O_centroid(elems, X):
    """Legacy function: get head point from O atoms specifically."""
    O = X[elems == "O"]
    if len(O) == 0:
        raise ValueError("No oxygen atoms found; cannot define head as O-centroid.")
    return O.mean(axis=0)

def head_vector(elems, X, eps=1e-12):
    """
    Calculate head vector from center-of-mass to head point.
    
    Robust implementation with fallback:
    1. Try primary heteroatom-based head point (O/N/S centroid)
    2. If vector norm is too small, fall back to farthest heteroatom from COM
    3. If still too small, use farthest atom overall
    
    This prevents degenerate cases where:
    - Molecule is symmetric (heteroatom centroid near COM)
    - Molecule has no heteroatoms
    - Multiple heteroatoms are distributed uniformly
    """
    c = com(X)
    
    # Try primary head point
    h = head_point_from_atoms(elems, X)
    vec = h - c
    norm_vec = np.linalg.norm(vec)
    
    # If vector is too small (symmetric or degenerate), use farthest atom
    if norm_vec < eps:
        # Find farthest heteroatom (O/N/S) from COM
        hetero_mask = np.isin(elems, ["O", "N", "S"])
        hetero_indices = np.where(hetero_mask)[0]
        
        if len(hetero_indices) > 0:
            hetero_coords = X[hetero_indices]
            dists = np.linalg.norm(hetero_coords - c, axis=1)
            farthest_hetero_idx = hetero_indices[np.argmax(dists)]
            h = X[farthest_hetero_idx]
            vec = h - c
            norm_vec = np.linalg.norm(vec)
        
        # Still degenerate? Use farthest atom overall
        if norm_vec < eps:
            dists = np.linalg.norm(X - c, axis=1)
            farthest_idx = np.argmax(dists)
            h = X[farthest_idx]
            vec = h - c
            norm_vec = np.linalg.norm(vec)
    
    return unit(vec, eps=eps), c, h

def pick_nearest_atoms_to_head(elems, X, head_point, k=4):
    """
    Return 0-based indices of k atoms closest to head_point.
    
    Priority:
    1. O atoms (oxygen - most common in acids)
    2. N atoms (nitrogen - amines, ammonium)
    3. S atoms (sulfur - thiols, sulfides)
    4. All atoms (if fewer heteroatoms exist than k)
    
    This replaces O-specific logic and works for any molecule type.
    """
    # Try oxygen atoms first
    idxO = np.where(elems == "O")[0]
    
    if len(idxO) >= k:
        Ocoords = X[idxO]
        d = np.linalg.norm(Ocoords - head_point, axis=1)
        return idxO[np.argsort(d)[:k]]
    
    # If not enough oxygen atoms, include other heteroatoms
    hetero_mask = np.isin(elems, ["O", "N", "S"])
    idx_hetero = np.where(hetero_mask)[0]
    
    if len(idx_hetero) >= k:
        hetero_coords = X[idx_hetero]
        d = np.linalg.norm(hetero_coords - head_point, axis=1)
        return idx_hetero[np.argsort(d)[:k]]
    
    # Fallback: use all atoms
    d = np.linalg.norm(X - head_point, axis=1)
    return np.argsort(d)[:k]

# Backward compatibility alias
def pick_nearest_O_indices(elems, X, head_point, k=4):
    """
    [DEPRECATED] Renamed to pick_nearest_atoms_to_head for clarity.
    Return 0-based indices of k atoms (O atoms preferred) closest to head_point.
    """
    return pick_nearest_atoms_to_head(elems, X, head_point, k=k)

# ---------------- clash checks ----------------
def min_dist_between_sets(XA, XB):
    d = XA[:,None,:] - XB[None,:,:]
    return np.linalg.norm(d, axis=2).min()

def clashes_with_cluster(Xnew, cluster_xyz_list, clash_cut):
    for Xold in cluster_xyz_list:
        if min_dist_between_sets(Xold, Xnew) < clash_cut:
            return True
    return False

# ---------------- xcontrol (distance restraints) ----------------
def write_xcontrol_distance_constraints(path, pairs_1based, target=3.40, k_fc=0.50):
    """Write distance constraints in xTB xcontrol format."""
    lines = []
    lines.append("$constrain")
    lines.append(f"  force constant={k_fc:.3f}")
    for (i,j) in pairs_1based:
        lines.append(f"  distance: {i}, {j}, {target:.3f}")
    lines.append("$end")
    Path(path).write_text("\n".join(lines) + "\n")

# ---------------- placing a new molecule: Head-Insert onto an anchor molecule ----------------
def place_head_insert(mono_elems, mono_X, anchor_elems, anchor_X, rng,
                      d_range=(8.0, 11.0),
                      lateral_range=1.5,
                      angle_jitter_deg=25.0):
    """
    Build Xnew (Natom,3) from monomer by aligning head vector towards anchor's head vector
    and translating to a controlled COM distance + lateral offset.
    
    Args:
        mono_elems: element array of monomer
        mono_X: coordinate array of monomer
        anchor_elems: element array of anchor molecule
        anchor_X: coordinate array of anchor molecule
        rng: random number generator
        d_range, lateral_range, angle_jitter_deg: placement parameters
    """
    vA, cA, hA = head_vector(anchor_elems, anchor_X)     # anchor's head direction
    vB, cB, hB = head_vector(mono_elems, mono_X)         # monomer's head direction (reference)

    # center new monomer at its COM
    Y = mono_X - cB

    # rotate so head faces anchor head: R vB = -vA
    R0 = rodrigues_align(vB, -vA)
    R = random_twist_about(vA, rng) @ R0

    # additional small jitter rotation (diversity)
    jitter_axis = random_unit(rng)
    jitter = axis_angle(jitter_axis, np.deg2rad(rng.uniform(-angle_jitter_deg, angle_jitter_deg)))
    R = jitter @ R

    # choose COM distance
    d = rng.uniform(*d_range)

    # lateral perturbation perpendicular to vA
    tmp = random_unit(rng)
    lat = tmp - np.dot(tmp, vA)*vA
    lat = unit(lat) * rng.uniform(0.0, lateral_range)

    Xnew = (Y @ R.T) + (cA + d*vA + lat)
    return Xnew, {"d": d}

# ---------------- building one cluster ----------------
def build_cluster_randomN(mono_elems, mono_X, rng,
                          N_fixed,
                          clash_cut,
                          max_trials_add,
                          d_range,
                          lateral_range,
                          angle_jitter_deg):
    """
    Build a cluster with a fixed number of monomers.

    Returns:
      cluster_list: list of (Natom,3) monomer geometries (placed)
      N: number of monomers (same as N_fixed)
      attach_edges: list of (anchor_idx, new_idx) edges describing growth
    """
    N = int(N_fixed)
    cluster = []
    attach_edges = []

    # place first molecule at origin (COM at 0)
    X0 = mono_X - com(mono_X)
    cluster.append(X0)

    for k in range(1, N):  # adding molecule k+1
        ok = False
        for _ in range(max_trials_add):
            anchor_idx = int(rng.integers(0, len(cluster)))
            anchor_X = cluster[anchor_idx]

            Xnew, info = place_head_insert(
                mono_elems, mono_X, mono_elems, anchor_X, rng,
                d_range=d_range,
                lateral_range=lateral_range,
                angle_jitter_deg=angle_jitter_deg
            )
            if not clashes_with_cluster(Xnew, cluster, clash_cut):
                cluster.append(Xnew)
                attach_edges.append((anchor_idx, k))  # anchor -> new
                ok = True
                break
        if not ok:
            return None, None, None  # fail
    return cluster, N, attach_edges

# ---------------- selecting O atoms and generating O...O restraint pairs ----------------
def build_OO_pairs_for_cluster(mono_elems, mono_X, cluster_list, attach_edges,
                              kO=4, n_pairs_per_edge=8, rng=None):
    """
    For each growth edge (anchor -> new), select kO nearest atoms to head in each monomer,
    then pick n_pairs_per_edge cross pairs to restrain (global numbering in the assembled cluster).
    
    Uses generalized head-point-based atom selection (not just O atoms).
    This works for molecules without oxygen or with O atoms far from the binding site.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    # Compute head point and select kO nearest atoms to it (robust method)
    _, _, h = head_vector(mono_elems, mono_X)
    atom_sel = pick_nearest_atoms_to_head(mono_elems, mono_X, h, k=kO)

    n_atoms = len(mono_elems)
    pairs = []

    # all possible cross pairs between selected atoms
    all_local_pairs = [(int(i), int(j)) for i in atom_sel for j in atom_sel]

    for (a_idx, new_idx) in attach_edges:
        # global offset for each monomer copy
        offA = a_idx * n_atoms
        offB = new_idx * n_atoms

        m = min(n_pairs_per_edge, len(all_local_pairs))
        chosen = rng.choice(len(all_local_pairs), size=m, replace=False)
        for t in chosen:
            i_local, j_local = all_local_pairs[t]
            i_global = offA + i_local + 1       # 1-based
            j_global = offB + j_local + 1       # 1-based
            pairs.append((i_global, j_global))

    return pairs

# ---------------- build OO pairs for mixed molecule clusters ----------------
def build_OO_pairs_for_mixed_cluster(molecules_data, cluster_list, molecule_indices, 
                                    attach_edges, kO=4, n_pairs_per_edge=8, rng=None):
    """
    Build restraint pairs for clusters with mixed molecule types.
    For each edge, select kO nearest atoms (robust head-based) from each molecule type.
    
    This replaces O-centroid-only logic with flexible head-point selection,
    allowing restraints for:
    - Molecules without oxygen atoms
    - Molecules with internal oxygen (O centroid near COM)
    - Any heteroatom-bearing molecule
    """
    if rng is None:
        rng = np.random.default_rng(0)
    
    pairs = []
    atom_sel_by_mol = []
    
    # Pre-compute atom indices for each molecule type (near head point)
    for mol_spec, elems, X in molecules_data:
        _, _, h = head_vector(elems, X)  # Robust head calculation
        atom_sel = pick_nearest_atoms_to_head(elems, X, h, k=kO)
        atom_sel_by_mol.append((atom_sel, len(elems)))  # Store both atom indices and atom count
    
    # Compute cumulative atom offsets for each position in cluster_list
    atom_offsets = [0]
    for i in range(len(cluster_list)):
        mol_idx = molecule_indices[i]
        _, natoms = atom_sel_by_mol[mol_idx]
        atom_offsets.append(atom_offsets[-1] + natoms)
    
    # For each edge, create cross-molecule restraints
    for (anchor_idx, new_idx) in attach_edges:
        anchor_mol_idx = molecule_indices[anchor_idx]
        new_mol_idx = molecule_indices[new_idx]
        
        atom_sel_anchor, _ = atom_sel_by_mol[anchor_mol_idx]
        atom_sel_new, _ = atom_sel_by_mol[new_mol_idx]
        
        anchor_offset = atom_offsets[anchor_idx]
        new_offset = atom_offsets[new_idx]
        
        # Create all possible cross pairs
        all_local_pairs = [(int(i), int(j)) for i in atom_sel_anchor for j in atom_sel_new]
        
        if len(all_local_pairs) == 0:
            continue
        
        m = min(n_pairs_per_edge, len(all_local_pairs))
        chosen = rng.choice(len(all_local_pairs), size=m, replace=False)
        
        for t in chosen:
            i_local, j_local = all_local_pairs[t]
            i_global = anchor_offset + i_local + 1  # 1-based (xcontrol format)
            j_global = new_offset + j_local + 1      # 1-based
            pairs.append((i_global, j_global))
    
    return pairs


# --- supporting function to compute atom offsets during assembly ---
def compute_atom_offsets(cluster_list, molecules_data, molecule_indices):
    """Compute atom offsets for each molecule in the assembled cluster."""
    offsets = [0]
    for i in range(len(cluster_list)):
        mol_idx = molecule_indices[i]
        mol_spec, elems, _ = molecules_data[mol_idx]
        offsets.append(offsets[-1] + len(elems))
    return offsets


def compute_molecule_radii(molecules_data):
    """
    Compute a characteristic radius for each molecule type: the maximum distance
    of any atom from its center-of-mass. Returns list of radii in the same order
    as `molecules_data`.

    The radius is used to adapt `dmin/dmax` and `lateral` for mixed pairs so
    larger molecules are placed at larger COM separations and vice versa.
    """
    radii = []
    for mol_spec, elems, X in molecules_data:
        c = com(X)
        r = float(np.max(np.linalg.norm(X - c, axis=1)))
        radii.append(r)
    return radii



# ---------------- single seed generation (worker function) ----------------
def _check_only_pt_molecules(molecules_data):
    """
    Check if all molecules are PT (Platinum cluster) molecules.
    Returns True only if there's a single molecule type named 'PT'.
    """
    if len(molecules_data) != 1:
        return False
    
    mol_spec, elems, X = molecules_data[0]
    return mol_spec.name.upper() == "PT"


def _generate_single_seed(s, cfg, molecules_data, rng_seed_offset):
    """
    Generate a single cluster seed with mixed molecules.
    
    Two modes:
    1. Pure PT mode: If only PT molecules, place them head-to-head
    2. Mixed mode: PT atoms placed first, then other molecules as close as possible
    
    Args:
        s: seed number
        cfg: ClusterConfig
        molecules_data: list of (MoleculeSpec, elems, X) tuples
        rng_seed_offset: offset for RNG seed
    
    Returns:
        dict with seed info, or None if failed
    """
    rng = np.random.default_rng(rng_seed_offset + s)
    
    # Check placement mode
    only_pt = _check_only_pt_molecules(molecules_data)
    
    # Build cluster by sequentially placing molecules
    cluster_list = []
    attach_edges = []
    molecule_indices = []  # track which molecule type each entry in cluster belongs to
    
    # Find PT molecule index if it exists
    pt_idx = None
    for i, (spec, _, _) in enumerate(molecules_data):
        if spec.name.upper() == "PT":
            pt_idx = i
            break
    
    # Determine first molecule to place
    if pt_idx is not None and not only_pt:
        # Mixed mode: start with PT molecule
        first_idx = pt_idx
    else:
        # Pure PT or no PT: start with first molecule
        first_idx = 0
    
    # Start with first molecule to place
    mol_spec, elems, X = molecules_data[first_idx]
    X0 = X - com(X)
    cluster_list.append(X0)
    molecule_indices.append(first_idx)
    
    # Track how many of each type have been placed
    placed_by_type = [0] * len(molecules_data)
    placed_by_type[first_idx] = 1
    
    # Place remaining molecules according to spec
    total_needed = sum(spec.count for spec, _, _ in molecules_data)
    count_placed = 1
    
    if only_pt:
        # Pure PT mode: place PT molecules using standard head-insert rule
        while count_placed < total_needed:
            mol_idx = 0  # Only one molecule type
            mol_spec, elems, X = molecules_data[mol_idx]
            
            ok = False
            for _ in range(cfg.max_trials_add):
                anchor_idx = int(rng.integers(0, len(cluster_list)))
                anchor_X = cluster_list[anchor_idx]
                anchor_elems = elems
                
                Xnew, info = place_head_insert(
                    elems, X, anchor_elems, anchor_X, rng,
                    d_range=(cfg.dmin, cfg.dmax),
                    lateral_range=cfg.lateral,
                    angle_jitter_deg=cfg.jitter_deg
                )
                if not clashes_with_cluster(Xnew, cluster_list, cfg.clash_cut):
                    cluster_list.append(Xnew)
                    molecule_indices.append(mol_idx)
                    attach_edges.append((anchor_idx, len(cluster_list) - 1))
                    placed_by_type[mol_idx] += 1
                    ok = True
                    count_placed += 1
                    break
            
            if not ok:
                return None  # fail
    
    else:
        # Mixed mode: place PT molecules first, then other molecules close to PT
        # Continue placing PT molecules
        while placed_by_type[pt_idx] < molecules_data[pt_idx][0].count:
            mol_spec, elems, X = molecules_data[pt_idx]
            
            ok = False
            for _ in range(cfg.max_trials_add):
                # Try anchoring to existing PT molecules first
                pt_anchors = [idx for idx in range(len(cluster_list)) 
                             if molecule_indices[idx] == pt_idx]
                
                if pt_anchors:
                    anchor_idx = rng.choice(pt_anchors)
                else:
                    anchor_idx = int(rng.integers(0, len(cluster_list)))
                
                anchor_X = cluster_list[anchor_idx]
                anchor_mol_idx = molecule_indices[anchor_idx]
                anchor_elems = molecules_data[anchor_mol_idx][1]
                
                Xnew, info = place_head_insert(
                    elems, X, anchor_elems, anchor_X, rng,
                    d_range=(cfg.dmin, cfg.dmax),
                    lateral_range=cfg.lateral,
                    angle_jitter_deg=cfg.jitter_deg
                )
                if not clashes_with_cluster(Xnew, cluster_list, cfg.clash_cut):
                    cluster_list.append(Xnew)
                    molecule_indices.append(pt_idx)
                    attach_edges.append((anchor_idx, len(cluster_list) - 1))
                    placed_by_type[pt_idx] += 1
                    ok = True
                    count_placed += 1
                    break
            
            if not ok:
                return None  # fail
        
        # Then place non-PT molecules close to PT
        while count_placed < total_needed:
            # Find first incomplete non-PT molecule type
            mol_idx = None
            for i in range(len(molecules_data)):
                if i != pt_idx and placed_by_type[i] < molecules_data[i][0].count:
                    mol_idx = i
                    break
            
            if mol_idx is None:
                break  # All molecules placed
            
            mol_spec, elems, X = molecules_data[mol_idx]
            
            ok = False
            for _ in range(cfg.max_trials_add):
                # Prefer anchoring to PT molecules
                pt_anchors = [idx for idx in range(len(cluster_list)) 
                             if molecule_indices[idx] == pt_idx]
                
                if pt_anchors:
                    anchor_idx = rng.choice(pt_anchors)
                else:
                    anchor_idx = int(rng.integers(0, len(cluster_list)))
                
                anchor_X = cluster_list[anchor_idx]
                anchor_mol_idx = molecule_indices[anchor_idx]
                anchor_elems = molecules_data[anchor_mol_idx][1]
                
                # Use slightly tighter distances for non-PT molecules (70-90% of standard)
                # to keep them close to PT but avoid clashes
                tighter_dmin = float(cfg.dmin) * 0.75
                tighter_dmax = float(cfg.dmax) * 0.85
                tighter_lateral = float(cfg.lateral) * 0.75
                
                Xnew, info = place_head_insert(
                    elems, X, anchor_elems, anchor_X, rng,
                    d_range=(tighter_dmin, tighter_dmax),
                    lateral_range=tighter_lateral,
                    angle_jitter_deg=cfg.jitter_deg
                )
                if not clashes_with_cluster(Xnew, cluster_list, cfg.clash_cut):
                    cluster_list.append(Xnew)
                    molecule_indices.append(mol_idx)
                    attach_edges.append((anchor_idx, len(cluster_list) - 1))
                    placed_by_type[mol_idx] += 1
                    ok = True
                    count_placed += 1
                    break
            
            if not ok:
                return None  # fail
    
    # Assemble mixed cluster
    elems_all = []
    X_all = []
    
    for i, X in enumerate(cluster_list):
        mol_idx = molecule_indices[i]
        mol_spec, elems, _ = molecules_data[mol_idx]
        elems_all.append(elems)
        X_all.append(X)
    
    elems_all = np.concatenate(elems_all)
    X_all = np.vstack(X_all)
    
    # Build restraint pairs for mixed molecules
    pairs = build_OO_pairs_for_mixed_cluster(
        molecules_data, cluster_list, molecule_indices, attach_edges,
        kO=cfg.kO, n_pairs_per_edge=cfg.pairs_per_edge, rng=rng
    )
    
    # Create composition string (actual counts from this seed)
    comp_parts = []
    for spec, _, _ in molecules_data:
        count = placed_by_type[molecules_data.index((spec, molecules_data[[m[0] for m in molecules_data].index(spec)][1], molecules_data[[m[0] for m in molecules_data].index(spec)][2]))]
    
    # Simpler approach: just use the counts from cfg.molecules
    comp_str = "_".join([f"{spec.name}{spec.count}" for spec in cfg.molecules])
    
    return {
        "s": s,
        "composition": comp_str,
        "molecule_indices": molecule_indices,
        "X_all": X_all,
        "pairs": pairs,
        "elems_all": elems_all
    }




def run_with_args_parsed(args, sampling_workers=None):
    """Run sampling with parsed argparse args object. Supports both legacy single-molecule and new multi-molecule modes."""
    if sampling_workers is None:
        sampling_workers = getattr(args, 'sampling_workers', 1)
    
    import concurrent.futures
    import os
    
    # Create config from parsed arguments (handles both legacy and new modes)
    cfg = create_config_from_args(args)
    
    # Load molecule data for each spec
    molecules_data = []
    for mol_spec in cfg.molecules:
        try:
            elems, X, comment = read_xyz(mol_spec.file)
            molecules_data.append((mol_spec, elems, X))
            print(f"[LOAD] {mol_spec.name}: {mol_spec.count}x from {mol_spec.file} ({len(elems)} atoms)")
        except Exception as e:
            print(f"[ERROR] Failed to load {mol_spec.file}: {e}")
            return 1
    
    out_dir = Path(cfg.out_dir)
    seeds_dir = out_dir / "seeds"
    seeds_dir.mkdir(parents=True, exist_ok=True)
    
    results = []
    n_ok = 0
    
    # Use multiprocessing for sampling if workers > 1
    if sampling_workers and sampling_workers > 1:
        print(f"[SAMPLING] Using {sampling_workers} worker processes")
        with concurrent.futures.ProcessPoolExecutor(max_workers=sampling_workers) as ex:
            futures = {
                ex.submit(_generate_single_seed, s, cfg, molecules_data, cfg.random_seed * 1000): s
                for s in range(1, cfg.nseeds + 1)
            }
            for fut in concurrent.futures.as_completed(futures):
                s = futures[fut]
                try:
                    result = fut.result()
                    if result is not None:
                        results.append(result)
                        n_ok += 1
                except Exception as e:
                    print(f"[WARN] seed {s:04d} failed: {e}")
    else:
        # Single-threaded fallback
        rng = np.random.default_rng(cfg.random_seed)
        for s in range(1, cfg.nseeds + 1):
            result = _generate_single_seed(s, cfg, molecules_data, cfg.random_seed * 1000)
            if result is None:
                print(f"[WARN] seed {s:04d} failed (could not place all molecules without clashes).")
                continue
            results.append(result)
            n_ok += 1
    
    # Write results to disk
    frames = []
    comments = []
    for result in results:
        s = result["s"]
        composition = result["composition"]
        X_all = result["X_all"]
        pairs = result["pairs"]
        elems_all = result["elems_all"]
        
        seed_name = f"seed_{s:04d}_{composition}"
        xyz_path = seeds_dir / f"{seed_name}.xyz"
        xc_path = seeds_dir / f"{seed_name}.xcontrol"
        
        write_xyz(xyz_path, elems_all, X_all, comment=f"{seed_name} HeadInsert mixed")
        write_xcontrol_distance_constraints(xc_path, pairs, target=cfg.OO_target, k_fc=cfg.k_fc)
        
        frames.append(X_all)
        comments.append(f"{seed_name}")
    
    # Write multi-frame file for all seeds with same composition
    if frames:
        # All seeds have same composition in a run, so write one multi-frame file
        nat = frames[0].shape[0]
        elems_template = results[0]["elems_all"]
        write_xyz_multiframe(out_dir / f"seeds_{cfg.molecules[0].name}_combined.xyz", 
                           frames, elems_template, comments)
    
    print(f"Done. Generated {n_ok} cluster seeds in {seeds_dir.resolve()}")
    print(f"Composition: {' + '.join([f'{s.count}{s.name}' for s in cfg.molecules])}")
    return 0



def run_with_args(argv=None):
    import concurrent.futures
    import os
    
    ap = get_argument_parser()
    args = ap.parse_args(argv)
    
    # Check for sampling_workers argument (not in base parser, add dynamically if needed)
    sampling_workers = getattr(args, 'sampling_workers', 1)
    
    return run_with_args_parsed(args, sampling_workers=sampling_workers)


def main():
    return run_with_args()


# Support running as: python -m lib.sampling
if __name__ == "__main__":
    main()
