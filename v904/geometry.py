"""Cube geometry: the V9.03 steric gap, hard-core overlap and cluster bonds.

Three different "distance" notions are kept apart on purpose:

* `support_gap`   V9.03's steric gap, measured along the centre line with the
                  rounded-cube support function.  Unchanged, so the steric
                  term is the V9.03 term.
* `cores_overlap` exact separating-axis test on the SHARP cubes the vdW
                  integral runs over.  The integral diverges if they overlap,
                  so overlap is a hard-core constraint (weight 0).
* `connected`     the cluster definition that makes the target normalisable
                  (see `model.ClusterEnsemble`).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import lsq_linear


def support_gap(r_vec, Q1, Q2, edge_m, round_m):
    """Rounded-cube gap along the centre line, vectorised; V9.03 formula.

    r_vec (P, 3) = p2 - p1.  support(Q, d) = (L/2 - r) |Q^T d|_1 + r.
    """
    r_vec = np.asarray(r_vec, float).reshape(-1, 3)
    Q1 = np.asarray(Q1, float).reshape(-1, 3, 3)
    Q2 = np.asarray(Q2, float).reshape(-1, 3, 3)
    dist = np.linalg.norm(r_vec, axis=1)
    safe = np.where(dist > 1e-15, dist, 1.0)
    rhat = r_vec / safe[:, None]
    l1_1 = np.abs(np.einsum("pab,pa->pb", Q1, rhat)).sum(axis=1)
    l1_2 = np.abs(np.einsum("pab,pa->pb", Q2, rhat)).sum(axis=1)
    core = 0.5 * edge_m - round_m
    gap = dist - (core * l1_1 + round_m) - (core * l1_2 + round_m)
    return np.where(dist > 1e-15, gap, -edge_m)


def cores_overlap(r_vec, Q1, Q2, edge_m, tol_m=0.0):
    """Separating-axis test for pairs of sharp cubes; True where they overlap.

    Fifteen candidate axes: three faces of each cube and nine edge crosses.
    Two convex polyhedra are disjoint iff one of these separates them.
    """
    r_vec = np.asarray(r_vec, float).reshape(-1, 3)
    Q1 = np.asarray(Q1, float).reshape(-1, 3, 3)
    Q2 = np.asarray(Q2, float).reshape(-1, 3, 3)
    half = 0.5 * edge_m
    a = np.swapaxes(Q1, 1, 2)                      # (P, 3, 3) rows = axes
    b = np.swapaxes(Q2, 1, 2)
    cross = np.cross(a[:, :, None, :], b[:, None, :, :]).reshape(-1, 9, 3)
    axes = np.concatenate([a, b, cross], axis=1)   # (P, 15, 3)
    norm = np.linalg.norm(axes, axis=2)
    valid = norm > 1e-9                            # parallel edges give 0
    axes = axes / np.where(valid, norm, 1.0)[..., None]
    ra = half * np.abs(np.einsum("pik,pjk->pij", axes, a)).sum(axis=2)
    rb = half * np.abs(np.einsum("pik,pjk->pij", axes, b)).sum(axis=2)
    sep = np.abs(np.einsum("pik,pk->pi", axes, r_vec)) - ra - rb
    separated = np.any(valid & (sep > tol_m), axis=1)
    return ~separated


def core_distance(p1, Q1, p2, Q2, edge_m):
    """Exact minimum distance between two sharp cubes (0 if they overlap).

    A bounded linear least-squares problem: minimise |p1 + Q1 u - p2 - Q2 v|
    over u, v in [-L/2, L/2]^3.  Diagnostics only.
    """
    half = 0.5 * edge_m
    M = np.hstack([np.asarray(Q1), -np.asarray(Q2)])
    res = lsq_linear(M, np.asarray(p2) - np.asarray(p1), bounds=(-half, half),
                     method="bvls", tol=1e-14)
    return float(np.linalg.norm(M @ res.x - (np.asarray(p2) - np.asarray(p1))))


def connected(pos, bond_m, anchor=0):
    """True if the bond graph (centre distance < bond_m) is one component."""
    pos = np.asarray(pos)
    d2 = np.einsum("ijk,ijk->ij", pos[:, None] - pos[None], pos[:, None] - pos[None])
    adj = d2 < bond_m * bond_m
    seen = np.zeros(len(pos), bool)
    seen[anchor] = True
    frontier = seen.copy()
    while frontier.any():
        new = adj[frontier].any(axis=0) & ~seen
        seen |= new
        frontier = new
    return bool(seen.all())
