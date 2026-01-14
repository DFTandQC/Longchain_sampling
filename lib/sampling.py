#!/usr/bin/env python3
from pathlib import Path
import numpy as np

from .config import ClusterConfig, get_argument_parser

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

# ---------------- PT head definition (no topology needed) ----------------
def head_point_O_centroid(elems, X):
    O = X[elems == "O"]
    if len(O) == 0:
        raise ValueError("No oxygen atoms found; cannot define head as O-centroid.")
    return O.mean(axis=0)

def head_vector(elems, X):
    c = com(X)
    h = head_point_O_centroid(elems, X)
    return unit(h - c), c, h

def pick_nearest_O_indices(elems, X, head_point, k=4):
    """Return 0-based indices of k oxygen atoms closest to head_point."""
    idxO = np.where(elems == "O")[0]
    Ocoords = X[idxO]
    d = np.linalg.norm(Ocoords - head_point, axis=1)
    sel = idxO[np.argsort(d)[:k]]
    return sel

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
def place_head_insert(mono_elems, mono_X, anchor_X, rng,
                      d_range=(8.0, 11.0),
                      lateral_range=1.5,
                      angle_jitter_deg=25.0):
    """
    Build Xnew (Natom,3) from monomer by aligning head vector towards anchor's head vector
    and translating to a controlled COM distance + lateral offset.
    """
    vA, cA, hA = head_vector(mono_elems, anchor_X)       # anchor's head direction
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
                mono_elems, mono_X, anchor_X, rng,
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
    For each growth edge (anchor -> new), select kO nearest O atoms to head in each monomer,
    then pick n_pairs_per_edge cross pairs to restrain (global numbering in the assembled cluster).
    """
    if rng is None:
        rng = np.random.default_rng(0)

    # indices (0-based) of O atoms closest to head within a monomer reference (same for all copies)
    h = head_point_O_centroid(mono_elems, mono_X)
    Osel = pick_nearest_O_indices(mono_elems, mono_X, h, k=kO)  # 0-based within monomer

    n_atoms = len(mono_elems)
    pairs = []

    # all possible cross pairs between selected O atoms
    all_local_pairs = [(int(i), int(j)) for i in Osel for j in Osel]

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

# ---------------- assemble cluster into one xyz ----------------
def assemble_cluster_xyz(mono_elems, cluster_list):
    elems_all = []
    X_all = []
    for X in cluster_list:
        elems_all.append(mono_elems)
        X_all.append(X)
    elems_all = np.concatenate(elems_all)
    X_all = np.vstack(X_all)
    return elems_all, X_all

# ---------------- main ----------------
def run_with_args(argv=None):
    ap = get_argument_parser()
    args = ap.parse_args(argv)

    # Create config from parsed arguments
    cfg = ClusterConfig.from_args(args)

    rng = np.random.default_rng(cfg.random_seed)

    mono_elems, mono_X, _ = read_xyz(cfg.mono_file)

    out_dir = Path(cfg.out_dir)
    seeds_dir = out_dir / "seeds"
    seeds_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    comments = []

    n_ok = 0
    for s in range(1, cfg.nseeds+1):
        cluster_list, N, edges = build_cluster_randomN(
            mono_elems, mono_X, rng,
            cfg.N,
            cfg.clash_cut,
            cfg.max_trials_add,
            d_range=(cfg.dmin, cfg.dmax),
            lateral_range=cfg.lateral,
            angle_jitter_deg=cfg.jitter_deg
        )
        if cluster_list is None:
            print(f"[WARN] seed {s:04d} failed (could not place all monomers without clashes).")
            continue

        # assemble xyz
        elems_all, X_all = assemble_cluster_xyz(mono_elems, cluster_list)

        # build O-O restraint pairs along growth edges
        pairs = build_OO_pairs_for_cluster(
            mono_elems, mono_X, cluster_list, edges,
            kO=cfg.kO, n_pairs_per_edge=cfg.pairs_per_edge,
            rng=rng
        )

        seed_name = f"seed_{s:04d}_N{N}"
        xyz_path = seeds_dir / f"{seed_name}.xyz"
        xc_path  = seeds_dir / f"{seed_name}.xcontrol"

        write_xyz(xyz_path, elems_all, X_all, comment=f"{seed_name} HeadInsert randomN")
        write_xcontrol_distance_constraints(xc_path, pairs, target=cfg.OO_target, k_fc=cfg.k_fc)

        frames.append(X_all)
        comments.append(f"{seed_name} (N={N})")
        n_ok += 1

    # multi-frame combined xyz (note: variable N means variable atom count; cannot mix in one xyz file safely)
    # So we write one multi-frame file PER N.
    byN = {}
    for X, cmt in zip(frames, comments):
        nat = X.shape[0]
        N = nat // len(mono_elems)
        byN.setdefault(N, {"frames": [], "comments": []})
        byN[N]["frames"].append(X)
        byN[N]["comments"].append(cmt)

    for N, pack in byN.items():
        elemsN = np.concatenate([mono_elems]*N)
        write_xyz_multiframe(out_dir / f"PTN_seeds_N{N}.xyz", pack["frames"], elemsN, pack["comments"])

    print(f"Done. Generated {n_ok} cluster seeds in {seeds_dir.resolve()}")
    print("Also wrote multi-frame xyz grouped by N: out/PTN_seeds_N{N}.xyz")


def main():
    return run_with_args()


# Support running as: python -m lib.sampling
if __name__ == "__main__":
    main()
