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


def get_heteroatom_head_region(elems, X, k=4):
    """Return indices (list) of up to k atoms nearest to the heteroatom head point.

    Prefers O/N/S atoms; falls back to all atoms if heteroatoms are absent.
    """
    try:
        head_pt = head_point_from_atoms(elems, X)
    except Exception:
        head_pt = com(X)
    return list(pick_nearest_atoms_to_head(elems, X, head_pt, k=k))

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

# -------- PT ester-core point computation --------
def compute_pt_ester_core(elems, X, k=8):
    """
    Compute PT ester-core point: centroid of k oxygen atoms closest to molecular COM.
    
    Args:
        elems: element array
        X: coordinate array
        k: number of O atoms to use for core (default 8)
    
    Returns:
        3D core point (centroid of selected O atoms), or fallback COM if insufficient O atoms
    """
    c_mol = com(X)
    o_indices = np.where(elems == "O")[0]
    
    if len(o_indices) < k:
        # If fewer than k oxygen atoms, use all of them
        if len(o_indices) > 0:
            return X[o_indices].mean(axis=0)
        else:
            # No oxygen atoms; fall back to COM
            return c_mol
    
    # Find k closest O atoms to COM
    o_coords = X[o_indices]
    dists = np.linalg.norm(o_coords - c_mol, axis=1)
    k_closest_indices = np.argsort(dists)[:k]
    
    return o_coords[k_closest_indices].mean(axis=0)

def place_pt_to_pt_core(mono_elems, mono_X, anchor_X, anchor_core, rng,
                        d_range=(8.0, 11.0),
                        lateral_range=1.5,
                        angle_jitter_deg=25.0,
                        pt_k=8):
    """
    Place a PT molecule using ester-core-to-core strategy.
    
    Args:
        mono_elems: element array of PT monomer
        mono_X: coordinate array of PT monomer
        anchor_X: coordinate array of anchor PT (already placed)
        anchor_core: ester-core point of anchor PT
        rng: random number generator
        d_range, lateral_range, angle_jitter_deg: placement parameters
        pt_k: number of O atoms for core definition
    
    Returns:
        Xnew: placed PT coordinates
        info: metadata dict
    """
    # Compute new PT's core point
    c_mono = com(mono_X)
    core_mono = compute_pt_ester_core(mono_elems, mono_X, k=pt_k)
    
    # Center at original COM
    Y = mono_X - c_mono
    
    # Random rotation
    axis = random_unit(rng)
    theta = rng.uniform(0, 2*np.pi)
    R = axis_angle(axis, theta)
    
    # Jitter rotation
    jitter_axis = random_unit(rng)
    jitter_theta = np.deg2rad(rng.uniform(-angle_jitter_deg, angle_jitter_deg))
    jitter_R = axis_angle(jitter_axis, jitter_theta)
    R = jitter_R @ R
    
    # Compute direction from anchor core to new position
    shell_dir = random_unit(rng)
    d = rng.uniform(*d_range)
    
    # Lateral offset
    tmp = random_unit(rng)
    lat = tmp - np.dot(tmp, shell_dir) * shell_dir
    lat = unit(lat) * rng.uniform(0.0, lateral_range)
    
    # Place rotated monomer so its core aligns with shell position
    new_core_pos = anchor_core + d * shell_dir + lat
    Xnew = (Y @ R.T) + new_core_pos
    
    return Xnew, {"d": d, "core_pos": new_core_pos}

def compute_global_pt_core_centroid(cluster_list, molecule_indices, molecules_data, pt_idx, pt_k=8):
    """
    Compute global PT core centroid (average of all placed PT core points).
    
    Args:
        cluster_list: list of placed molecule coordinates
        molecule_indices: list indicating which molecule type each cluster entry is
        molecules_data: list of (MoleculeSpec, elems, X) tuples
        pt_idx: index of PT in molecules_data
        pt_k: number of O atoms for core definition
    
    Returns:
        3D global core centroid, or None if no PT placed yet
    """
    pt_cores = []
    
    for i, X in enumerate(cluster_list):
        if molecule_indices[i] == pt_idx:
            mol_spec, elems, _ = molecules_data[pt_idx]
            core = compute_pt_ester_core(elems, X, k=pt_k)
            pt_cores.append(core)
    
    if len(pt_cores) == 0:
        return None
    
    return np.mean(pt_cores, axis=0)

# -------- filtering functions --------
def radius_of_gyration(X):
    """Compute radius of gyration for atomic coordinates."""
    center = np.mean(X, axis=0)
    return np.sqrt(np.mean(np.sum((X - center)**2, axis=1)))

def count_inter_molecular_contacts(X, elems, molecule_indices, contact_cut=3.8):
    """
    Count inter-molecular atom pairs within contact_cut distance.
    
    Args:
        X: (N, 3) atomic coordinates
        elems: element symbols
        molecule_indices: list mapping each atom to molecule index
        contact_cut: contact distance threshold in Angstrom
    
    Returns:
        Number of inter-molecular contact pairs
    """
    n_contacts = 0
    
    for i in range(len(X)):
        for j in range(i+1, len(X)):
            if molecule_indices[i] != molecule_indices[j]:
                dist = np.linalg.norm(X[i] - X[j])
                if dist < contact_cut:
                    n_contacts += 1
    
    return n_contacts

def compute_rmsd_fingerprint(X):
    """
    Compute a structural fingerprint (histogram of inter-atomic distances).
    Returns a hashable tuple for duplicate detection.
    """
    dists = []
    for i in range(len(X)):
        for j in range(i+1, len(X)):
            dists.append(np.linalg.norm(X[i] - X[j]))
    
    dists = np.array(dists)
    # Histogram with 0.5 Å bins from 0 to 30 Å
    hist, _ = np.histogram(dists, bins=np.arange(0, 30.5, 0.5))
    return tuple(hist)

def filter_accept(cluster_list, molecule_indices, molecules_data, cfg, seen_fingerprints=None):
    """
    Check if a generated cluster seed passes all quality filters.
    
    Args:
        cluster_list: list of coordinate arrays for each placed molecule
        molecule_indices: list indicating which molecule type each entry is
        molecules_data: list of (MoleculeSpec, elems, X_template) tuples
        cfg: ClusterConfig with filter parameters
        seen_fingerprints: set of already-seen fingerprints for dedup (updated in-place)
    
    Returns:
        (pass: bool, reject_reason: str or None)
    """
    if not cfg.enable_filter:
        return True, None
    
    if seen_fingerprints is None:
        seen_fingerprints = set()
    
    # Concatenate all atoms
    X_list = []
    elems_list = []
    mol_idx_expanded = []
    
    for i, X_mol in enumerate(cluster_list):
        X_list.append(X_mol)
        mol_spec, elems, _ = molecules_data[molecule_indices[i]]
        elems_list.extend(elems)
        mol_idx_expanded.extend([i] * len(X_mol))
    
    X_all = np.vstack(X_list)
    
    # Check 1: Radius of gyration
    rg = radius_of_gyration(X_all)
    if rg > cfg.max_rg:
        return False, f"rg_too_large(rg={rg:.2f})"
    
    # Check 2: Inter-molecular contacts
    n_contacts = count_inter_molecular_contacts(X_all, elems_list, mol_idx_expanded, cfg.contact_cut)
    if n_contacts < cfg.min_contacts:
        return False, f"insufficient_contacts(n={n_contacts})"
    
    # Check 3: Duplicate detection (RMSD fingerprint)
    if cfg.rmsd_dedup is not None:
        fp = compute_rmsd_fingerprint(X_all)
        if fp in seen_fingerprints:
            return False, "duplicate"
        seen_fingerprints.add(fp)
    
    # Check 4: PT-specific rules (if applicable)
    # Only check if we have mixed molecules or multiple PT molecules with PT core placement enabled
    if len(molecules_data) > 1 or (len(cluster_list) > 1 and any(m[0].name == "PT" for m in molecules_data)):
        pt_idx = None
        for idx, (mol_spec, _, _) in enumerate(molecules_data):
            if mol_spec.name == "PT":
                pt_idx = idx
                break
        
        if pt_idx is not None:
            # Compute all PT core positions
            pt_cores = []
            pt_head_positions = []
            
            for i, X_mol in enumerate(cluster_list):
                if molecule_indices[i] == pt_idx:
                    mol_spec, elems, _ = molecules_data[pt_idx]
                    pt_core = compute_pt_ester_core(elems, X_mol, k=cfg.pt_k)
                    pt_cores.append(pt_core)
                    
                    # Get head region (nearest kO heteroatoms)
                    head_indices = get_heteroatom_head_region(elems, X_mol, k=cfg.kO)
                    if head_indices:
                        head_coords = X_mol[head_indices]
                        head_centroid = np.mean(head_coords, axis=0)
                        pt_head_positions.append(head_centroid)
            
            # PT core-to-core distance check
            if len(pt_cores) > 1:
                for i in range(len(pt_cores)):
                    for j in range(i+1, len(pt_cores)):
                        dist = np.linalg.norm(pt_cores[i] - pt_cores[j])
                        if dist > cfg.core_dist_max:
                            return False, f"pt_core_dist_too_large(d={dist:.2f})"
            
            # Non-PT head-to-PT-core distance check (for mixed systems)
            if pt_idx is not None and len(molecules_data) > 1 and len(pt_cores) > 0:
                pt_core_global = np.mean(pt_cores, axis=0)
                
                for i, X_mol in enumerate(cluster_list):
                    if molecule_indices[i] != pt_idx:
                        # This is a non-PT molecule; check its head to PT core
                        mol_spec, elems, _ = molecules_data[molecule_indices[i]]
                        head_indices = get_heteroatom_head_region(elems, X_mol, k=cfg.kO)
                        if head_indices:
                            head_coords = X_mol[head_indices]
                            head_centroid = np.mean(head_coords, axis=0)
                            dist = np.linalg.norm(head_centroid - pt_core_global)
                            if dist > cfg.head_core_max:
                                return False, f"head_core_dist_too_large(d={dist:.2f})"
    
    return True, None

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
    
    Three modes:
    1. Pure PT (core-based): PT core-to-core placement for better ester-region clustering
    2. Mixed with PT (core-based): PT molecules placed core-to-core, then non-PT near global PT core
    3. No PT (standard): Original head-insert behavior
    
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
    
    if pt_idx is not None:
        # PT-aware placement: use core-to-core strategy for PT molecules
        # Continue placing PT molecules using ester-core strategy
        while placed_by_type[pt_idx] < molecules_data[pt_idx][0].count:
            mol_spec, elems, X = molecules_data[pt_idx]
            
            ok = False
            for _ in range(cfg.max_trials_add):
                # Choose anchor: prefer existing PT, else any molecule
                pt_anchors = [idx for idx in range(len(cluster_list)) 
                             if molecule_indices[idx] == pt_idx]
                
                if pt_anchors:
                    anchor_idx = rng.choice(pt_anchors)
                else:
                    anchor_idx = int(rng.integers(0, len(cluster_list)))
                
                anchor_X = cluster_list[anchor_idx]
                anchor_mol_idx = molecule_indices[anchor_idx]
                
                # Use PT core-to-core placement if both are PT
                if anchor_mol_idx == pt_idx:
                    anchor_core = compute_pt_ester_core(elems, anchor_X, k=getattr(cfg, 'pt_k', 8))
                    # Keep PT cores close; make dmin/dmax compatible with core_dist_max to avoid impossible ranges
                    if cfg.core_dist_max > 0:
                        dmax_pt = min(cfg.dmax, cfg.core_dist_max)
                        dmin_pt = min(cfg.dmin, cfg.core_dist_max * 0.9)
                        if dmin_pt > dmax_pt:
                            dmin_pt = max(0.5 * dmax_pt, 0.1)  # keep a sane lower bound
                    else:
                        dmin_pt = cfg.dmin
                        dmax_pt = cfg.dmax
                    Xnew, info = place_pt_to_pt_core(
                        elems, X, anchor_X, anchor_core, rng,
                        d_range=(dmin_pt, dmax_pt),
                        lateral_range=cfg.lateral,
                        angle_jitter_deg=cfg.jitter_deg,
                        pt_k=getattr(cfg, 'pt_k', 8)
                    )
                else:
                    # Anchor is non-PT, use standard head-insert
                    anchor_elems = molecules_data[anchor_mol_idx][1]
                    Xnew, info = place_head_insert(
                        elems, X, anchor_elems, anchor_X, rng,
                        d_range=(cfg.dmin, cfg.dmax),
                        lateral_range=cfg.lateral,
                        angle_jitter_deg=cfg.jitter_deg
                    )
                
                # PT sanity: new PT core must be close to all existing PT cores
                new_core = compute_pt_ester_core(elems, Xnew, k=getattr(cfg, 'pt_k', 8))
                too_far = False
                for idx_existing, X_existing in enumerate(cluster_list):
                    if molecule_indices[idx_existing] == pt_idx:
                        existing_core = compute_pt_ester_core(molecules_data[pt_idx][1], X_existing, k=getattr(cfg, 'pt_k', 8))
                        if np.linalg.norm(new_core - existing_core) > cfg.core_dist_max:
                            too_far = True
                            break
                if too_far:
                    continue
                
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
        
        # After all PT placed, place non-PT molecules near global PT core
        if not only_pt and pt_idx is not None:
            # Compute global PT core centroid
            global_pt_core = compute_global_pt_core_centroid(
                cluster_list, molecule_indices, molecules_data, pt_idx,
                pt_k=getattr(cfg, 'pt_k', 8)
            )
            
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
                    
                    # Use slightly tighter distances for non-PT molecules to keep close to PT
                    tighter_dmin = float(cfg.dmin) * 0.75
                    tighter_dmax = float(cfg.dmax) * 0.85
                    tighter_lateral = float(cfg.lateral) * 0.75
                    
                    Xnew, info = place_head_insert(
                        elems, X, anchor_elems, anchor_X, rng,
                        d_range=(tighter_dmin, tighter_dmax),
                        lateral_range=tighter_lateral,
                        angle_jitter_deg=cfg.jitter_deg
                    )
                    
                    # Ensure non-PT head is close to PT core centroid
                    pt_cores_now = []
                    for idx_existing, X_existing in enumerate(cluster_list):
                        if molecule_indices[idx_existing] == pt_idx:
                            existing_core = compute_pt_ester_core(molecules_data[pt_idx][1], X_existing, k=getattr(cfg, 'pt_k', 8))
                            pt_cores_now.append(existing_core)
                    if pt_cores_now:
                        pt_core_global = np.mean(pt_cores_now, axis=0)
                        head_idx = get_heteroatom_head_region(elems, Xnew, k=cfg.kO)
                        if head_idx:
                            head_centroid = Xnew[head_idx].mean(axis=0)
                        else:
                            head_centroid = com(Xnew)
                        if np.linalg.norm(head_centroid - pt_core_global) > cfg.head_core_max:
                            continue
                    
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
        # No PT: use standard head-insert for all remaining molecules
        while count_placed < total_needed:
            # Cycle through molecule types in order
            mol_idx = None
            for i in range(len(molecules_data)):
                if placed_by_type[i] < molecules_data[i][0].count:
                    mol_idx = i
                    break
            
            if mol_idx is None:
                break
            
            mol_spec, elems, X = molecules_data[mol_idx]
            
            ok = False
            for _ in range(cfg.max_trials_add):
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
    
    # Create composition string based on actually placed counts
    comp_str = "_".join([
        f"{molecules_data[i][0].name}{placed_by_type[i]}" for i in range(len(molecules_data))
    ])
    
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
    
    # Initialize filter deduplication set if needed
    seen_fingerprints = set() if cfg.enable_filter and cfg.rmsd_dedup is not None else None
    
    # Track filter statistics
    filter_stats = {
        'rejected_by_rg': 0,
        'rejected_by_contacts': 0,
        'rejected_by_duplicate': 0,
        'accepted': 0
    }
    
    # Retry loop: generate until we have nseeds valid seeds or max_attempts reached
    attempt = 0
    max_attempts = cfg.max_attempts if cfg.enable_filter else 1
    nseeds_target = cfg.nseeds
    
    while n_ok < nseeds_target and attempt < max_attempts:
        attempt += 1
        
        # Generate seeds for this batch
        batch_results = []
        
        # Use multiprocessing for sampling if workers > 1
        if sampling_workers and sampling_workers > 1:
            if attempt == 1:
                print(f"[SAMPLING] Using {sampling_workers} worker processes")
            
            # For this batch, try to generate nseeds_target - n_ok more seeds
            seeds_to_generate = min(nseeds_target - n_ok, 100)  # Batch up to 100 per round
            
            with concurrent.futures.ProcessPoolExecutor(max_workers=sampling_workers) as ex:
                futures = {
                    ex.submit(_generate_single_seed, s + attempt*10000, cfg, molecules_data, cfg.random_seed * 1000 + attempt): s + attempt*10000
                    for s in range(1, seeds_to_generate + 1)
                }
                for fut in concurrent.futures.as_completed(futures):
                    s = futures[fut]
                    try:
                        result = fut.result()
                        if result is not None:
                            batch_results.append(result)
                    except Exception as e:
                        pass
        else:
            # Single-threaded fallback
            rng = np.random.default_rng(cfg.random_seed + attempt)
            seeds_to_generate = min(nseeds_target - n_ok, 100)
            
            for s in range(1, seeds_to_generate + 1):
                result = _generate_single_seed(s + attempt*10000, cfg, molecules_data, cfg.random_seed * 1000 + attempt)
                if result is None:
                    continue
                batch_results.append(result)
        
        # Filter batch results
        for result in batch_results:
            if n_ok >= nseeds_target:
                break
            
            # Extract cluster info for filtering
            X_all = result["X_all"]
            elems_all = result["elems_all"]
            molecule_indices = result["molecule_indices"]
            
            # Perform filtering
            cluster_list = []
            current_idx = 0
            for mol_idx in molecule_indices:
                mol_spec, elems, _ = molecules_data[mol_idx]
                n_atoms = len(elems)
                cluster_list.append(X_all[current_idx:current_idx + n_atoms])
                current_idx += n_atoms
            
            passes_filter, reject_reason = filter_accept(
                cluster_list, molecule_indices, molecules_data, cfg, seen_fingerprints
            )
            
            if passes_filter:
                results.append(result)
                filter_stats['accepted'] += 1
                n_ok += 1
            else:
                # Track rejection reason
                if reject_reason:
                    if 'rg_too_large' in reject_reason:
                        filter_stats['rejected_by_rg'] += 1
                    elif 'insufficient_contacts' in reject_reason:
                        filter_stats['rejected_by_contacts'] += 1
                    elif 'duplicate' in reject_reason:
                        filter_stats['rejected_by_duplicate'] += 1
                    elif 'pt_' in reject_reason or 'core_dist' in reject_reason or 'head_core' in reject_reason:
                        filter_stats['rejected_by_rg'] += 1  # PT checks are geometry-based
        
        if cfg.enable_filter and attempt < max_attempts and n_ok < nseeds_target:
            print(f"[FILTERING] Attempt {attempt}: {n_ok}/{nseeds_target} valid seeds accepted", end="")
            if reject_reason:
                print(f" (last rejection: {reject_reason})")
            else:
                print()
    
    # Report filtering statistics
    if cfg.enable_filter:
        print(f"\n[FILTER STATS]")
        print(f"  Accepted: {filter_stats['accepted']}")
        print(f"  Rejected (Rg): {filter_stats['rejected_by_rg']}")
        print(f"  Rejected (Contacts): {filter_stats['rejected_by_contacts']}")
        print(f"  Rejected (Duplicate): {filter_stats['rejected_by_duplicate']}")
        if n_ok < nseeds_target:
            print(f"  Warning: Only {n_ok}/{nseeds_target} valid seeds after {attempt} attempts")
    
    # Write results to disk
    frames = []
    comments = []
    
    # Apply keep_best filter if specified
    if cfg.keep_best is not None and len(results) > cfg.keep_best:
        # Sort by number of contacts (descending)
        scored_results = []
        for result in results:
            X_all = result["X_all"]
            elems_all = result["elems_all"]
            molecule_indices = result["molecule_indices"]
            
            cluster_list = []
            current_idx = 0
            for mol_idx in molecule_indices:
                mol_spec, elems, _ = molecules_data[mol_idx]
                n_atoms = len(elems)
                cluster_list.append(X_all[current_idx:current_idx + n_atoms])
                current_idx += n_atoms
            
            n_contacts = count_inter_molecular_contacts(X_all, elems_all, molecule_indices, cfg.contact_cut)
            scored_results.append((n_contacts, result))
        
        scored_results.sort(reverse=True, key=lambda x: x[0])
        results = [r[1] for r in scored_results[:cfg.keep_best]]
    
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
