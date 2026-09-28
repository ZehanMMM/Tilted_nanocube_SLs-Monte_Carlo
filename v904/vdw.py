"""Sharp-cube Hamaker energy, the V10 integral taken to its voxel limit.

The integral is the one V6.20/V9/V10 all evaluate,

    U = -(A / pi^2) * int_{cube 1} int_{cube 2} dV1 dV2 / |r1 - r2|^6,

over two sharp cubes of edge L, each at its own position and orientation.
V9.03 approximated it with a 4^3 midpoint voxel sum and a d^2 >= 1e-19 m^2
floor; V10's latest version (`rotational_entropy_model.central_vdw`) refines
the same voxel sum without the floor until it stops changing.  At V9's
1.8 nm gap the 4^3 sum is 2.9x too weak, so the refinement matters.

Three evaluators of the SAME integral live here:

* `uniform_voxel_energy`         n^3 midpoint voxels per cube, no floor.  The
                                 literal V10 method for arbitrary orientation;
                                 O(n^6), for validation only.
* `parallel_voxel_energy`        V10's displacement-multiplicity form of the
                                 same sum, exact for co-oriented cubes.
* `HamakerEvaluator`             the n -> infinity limit, fast enough for MC.

The fast evaluator rests on two exact steps.  With R = r2 - r1,

    r^-6 = -(1/3) div_2 (R R^-6)      and      (n2.R) R^-6 = div_1 (n2 R^-4 / 4),

so Gauss's theorem applied once per cube turns the 6-D volume integral into

    int int dV1 dV2 / r^6 = -(1/12) sum_{faces i, j} (n_i . n_j)
                                     int_{F_i} int_{F_j} dS_i dS_j / R^4 ,

and the inner integral of R^-4 over a square has a closed form
(`square_inverse_quartic`).  Only a 2-D Gauss rule over the outer face is
left, and the integrand is smooth whenever the cubes do not touch.  Nothing is
fitted or tabulated; the tests show it is the limit of the voxel sums.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

HAMAKER_J = 2.0e-20
EDGE_M = 16.0e-9


# --------------------------------------------------------------------------
# reference evaluators (validation only)
# --------------------------------------------------------------------------
def voxel_offsets(n: int, edge_m: float = EDGE_M) -> np.ndarray:
    """Body-frame centres of the n^3 midpoint voxels, shape (n^3, 3)."""
    step = edge_m / n
    lin = np.linspace(-edge_m / 2 + step / 2, edge_m / 2 - step / 2, n)
    grid = np.stack(np.meshgrid(lin, lin, lin, indexing="ij"), axis=-1)
    return grid.reshape(-1, 3)


def uniform_voxel_energy(p1, Q1, p2, Q2, n: int, edge_m: float = EDGE_M,
                         hamaker_J: float = HAMAKER_J, d2_floor_m2: float = 0.0,
                         chunk: int = 2048) -> float:
    """n^3 x n^3 midpoint voxel sum for arbitrary orientations.

    d2_floor_m2 = 1e-19 and n = 4 reproduce V9.03 exactly; the V10 refinement
    uses d2_floor_m2 = 0.  Memory is bounded by `chunk` voxels at a time.
    """
    base = voxel_offsets(n, edge_m)
    v1 = base @ np.asarray(Q1).T + np.asarray(p1)
    v2 = base @ np.asarray(Q2).T + np.asarray(p2)
    volume = (edge_m / n) ** 3
    total = 0.0
    for start in range(0, len(v1), chunk):
        block = v1[start:start + chunk]
        d = v2[None, :, :] - block[:, None, :]
        d2 = np.einsum("ijk,ijk->ij", d, d)
        total += float(np.sum(1.0 / np.maximum(d2, d2_floor_m2) ** 3))
    return -(hamaker_J / np.pi ** 2) * volume ** 2 * total


def parallel_voxel_energy(center_vector_m, n: int, edge_m: float = EDGE_M,
                          hamaker_J: float = HAMAKER_J) -> float:
    """V10 `central_vdw`: the same voxel sum for co-oriented cubes.

    For co-oriented cubes every voxel pair is labelled by one integer
    displacement k, occurring prod(n - |k|) times, so n^6 terms collapse to
    (2n - 1)^3.  `center_vector_m` is in the shared body frame.
    """
    idx = np.arange(1 - n, n)
    disp = np.stack(np.meshgrid(idx, idx, idx, indexing="ij"), -1).reshape(-1, 3)
    weight = np.prod(n - np.abs(disp), axis=1).astype(float) / n ** 6
    dr = np.asarray(center_vector_m)[None, :] - disp * (edge_m / n)
    return float(-(hamaker_J / np.pi ** 2) * edge_m ** 6
                 * np.sum(weight / np.einsum("ij,ij->i", dr, dr) ** 3))


def richardson(values, ns, power: float = 2.0) -> float:
    """Extrapolate the last two midpoint sums to h -> 0 assuming error ~ h^p."""
    (v1, v2), (n1, n2) = values[-2:], ns[-2:]
    r = (n2 / n1) ** power
    return float((r * v2 - v1) / (r - 1.0))


# --------------------------------------------------------------------------
# closed-form inner integral
# --------------------------------------------------------------------------
def _primitive(x, y, c2):
    """F with d2F/dxdy = 1 / (x^2 + y^2 + c^2)^2."""
    sx = np.sqrt(x * x + c2)
    sy = np.sqrt(y * y + c2)
    return (x / sx * np.arctan(y / sx) + y / sy * np.arctan(x / sy)) / (2.0 * c2)


def square_inverse_quartic(u0, v0, z2, half):
    """int_{-h}^{h} int_{-h}^{h} du dv / ((u-u0)^2 + (v-v0)^2 + z^2)^2.

    Four corners of `_primitive`, with the square roots shared between
    corners that have the same x or y (4 instead of 8).
    """
    xp, xm = half - u0, -half - u0
    yp, ym = half - v0, -half - v0
    sxp, sxm = np.sqrt(xp * xp + z2), np.sqrt(xm * xm + z2)
    syp, sym = np.sqrt(yp * yp + z2), np.sqrt(ym * ym + z2)
    ax_p, ax_m = xp / sxp, xm / sxm
    ay_p, ay_m = yp / syp, ym / sym
    total = (ax_p * (np.arctan(yp / sxp) - np.arctan(ym / sxp))
             - ax_m * (np.arctan(yp / sxm) - np.arctan(ym / sxm))
             + ay_p * (np.arctan(xp / syp) - np.arctan(xm / syp))
             - ay_m * (np.arctan(xp / sym) - np.arctan(xm / sym)))
    return total / (2.0 * z2)


def face_rule(q: int, half: float):
    """q x q Gauss-Legendre nodes (u, v) and weights on [-half, half]^2."""
    x, w = np.polynomial.legendre.leggauss(q)
    x, w = x * half, w * half
    u, v = np.meshgrid(x, x, indexing="ij")
    return u.ravel(), v.ravel(), np.outer(w, w).ravel()


def cube_faces(pos, Q, half):
    """Face centres, normals and in-plane axes for a batch of cubes.

    pos (P, 3), Q (P, 3, 3) -> four arrays of shape (P, 6, 3).
    Face 2k is +axis k, face 2k+1 is -axis k.
    """
    pos = np.asarray(pos, float).reshape(-1, 3)
    Q = np.asarray(Q, float).reshape(-1, 3, 3)
    cols = np.swapaxes(Q, 1, 2)                     # cols[p, k] = Q[p][:, k]
    sign = np.array([1.0, -1.0])
    normal = (sign[None, None, :, None] * cols[:, :, None, :]).reshape(-1, 6, 3)
    t1 = np.repeat(cols[:, [1, 2, 0]], 2, axis=1)
    t2 = np.repeat(cols[:, [2, 0, 1]], 2, axis=1)
    centre = pos[:, None, :] + half * normal
    return centre, normal, t1, t2


# --------------------------------------------------------------------------
# fast converged evaluator
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class HamakerSettings:
    """Quadrature controls.  Defaults are checked in tests/test_vdw.py.

    Face pairs are sorted by a lower bound on their separation (units of L):

    * below near_L   adaptive: a panel of the outer face is split into four
                     until it is small compared with its distance to the
                     inner face (radius <= eta_geometric x distance), or
                     moderately small (<= eta_value x distance) AND its
                     q_near-point value agrees with its four children's to
                     `tol_kBT` (scaled by the panel's area share).  Sharp
                     edges and corners of two cubes can come within a few
                     tenths of a nm while the rounded steric gap still reads
                     1.8 nm, so no fixed rule is safe here, and a value-only
                     test is not guaranteed either: a spike narrower than the
                     node spacing is invisible to parent and children alike.
                     The geometric condition cannot be fooled that way.
                     (A production chain did reach a -458 kBT pair with sharp
                     cores 1.5e-5 nm apart; the value-blind `octree_energy`
                     reproduced it, so that was the true sharp-cube integral
                     -- a model defect cured by the D0 core minimum in
                     model.ClusterEnsemble -- not a quadrature error.)
    * below mid_L    fixed q_mid-point Gauss per face edge.
    * beyond         fixed q_far-point Gauss per face edge.

    Cube pairs whose centres are more than far_pair_L apart use a volume
    Gauss rule of `volume_order`^3 nodes per cube instead, where the surface
    form would lose digits to cancellation between opposite faces.
    """
    q_near: int = 6
    q_mid: int = 8
    q_far: int = 4
    near_L: float = 0.35
    mid_L: float = 1.0
    far_pair_L: float = 2.6
    volume_order: int = 4
    z_floor_L: float = 1e-6
    tol_kBT: float = 1e-5
    temperature_K: float = 298.15
    max_level: int = 14
    min_level: int = 0
    # 'hybrid' (default): a panel of radius rho at distance >= d from the other
    # face is accepted if rho <= eta_geometric d, or if rho <= eta_value d AND
    # its value agrees with its four children's.  'children' / 'embedded' are
    # value-only tests kept for comparison: 'embedded' measurably fails on the
    # contact benchmark (4e-4 kBT), and neither is guaranteed against spikes.
    scheme: str = "hybrid"
    eta_geometric: float = 0.5
    eta_value: float = 4.0


class HamakerEvaluator:
    """Converged sharp-cube Hamaker pair energies, batched over pairs."""

    def __init__(self, edge_m: float = EDGE_M, hamaker_J: float = HAMAKER_J,
                 settings: HamakerSettings = HamakerSettings()):
        self.edge_m = float(edge_m)
        self.half = 0.5 * self.edge_m
        self.hamaker_J = float(hamaker_J)
        self.settings = settings
        self._rules = {q: face_rule(q, self.half)
                       for q in {settings.q_near, settings.q_mid, settings.q_far}}
        x, w = np.polynomial.legendre.leggauss(settings.volume_order)
        x, w = x * self.half, w * self.half
        g = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1).reshape(-1, 3)
        self._vol_nodes = g
        self._vol_weights = np.einsum("i,j,k->ijk", w, w, w).ravel()
        self._z_floor2 = (settings.z_floor_L * self.edge_m) ** 2
        self._prefactor = self.hamaker_J / (12.0 * np.pi ** 2)
        # two unit-panel Gauss rules on [-1, 1]^2 (q_near and q_near - 2),
        # evaluated together; their difference is the panel's error estimate
        units = []
        for q in (settings.q_near, settings.q_near - 2):
            x, w = np.polynomial.legendre.leggauss(q)
            uu, vv = np.meshgrid(x, x, indexing="ij")
            units.append((uu.ravel(), vv.ravel(), np.outer(w, w).ravel()))
        self._n_hi = len(units[0][0])
        self._unit_u = np.concatenate([units[0][0], units[1][0]])
        self._unit_v = np.concatenate([units[0][1], units[1][1]])
        self._w_hi, self._w_lo = units[0][2], units[1][2]
        self._child = np.array([[-0.5, -0.5], [-0.5, 0.5], [0.5, -0.5], [0.5, 0.5]])
        # panels that reached max_level without passing their test; a correct
        # run keeps this at 0 (the core-distance constraint guarantees it)
        self.max_level_hits = 0

    def _inner(self, nodes, c2, a2, b2, n2):
        """Closed-form int_{F_j} R^-4 dS at every node; nodes (F, n, 3)."""
        d = nodes - c2[:, None, :]
        u0 = np.einsum("fnk,fk->fn", d, a2)
        v0 = np.einsum("fnk,fk->fn", d, b2)
        z = np.einsum("fnk,fk->fn", d, n2)
        return square_inverse_quartic(u0, v0, np.maximum(z * z, self._z_floor2), self.half)

    def _panel_values(self, fp, uc, vc, hw, face):
        """High- and low-order Gauss estimates over panels (centre uc, vc; half-width hw)."""
        c1, a1, b1, c2, a2, b2, n2 = face
        u = uc[:, None] + hw[:, None] * self._unit_u[None, :]
        v = vc[:, None] + hw[:, None] * self._unit_v[None, :]
        nodes = (c1[fp][:, None, :] + u[..., None] * a1[fp][:, None, :]
                 + v[..., None] * b1[fp][:, None, :])
        inner = self._inner(nodes, c2[fp], a2[fp], b2[fp], n2[fp])
        area = hw * hw
        return (inner[:, :self._n_hi] @ self._w_hi) * area,             (inner[:, self._n_hi:] @ self._w_lo) * area

    def _adaptive(self, face, cos):
        """Adaptive integral of the inner closed form over each near face pair.

        scheme 'children': a panel's q_near value is compared with the sum of
        its four children's; if they agree to the tolerance times the panel's
        area share the children's sum is kept, otherwise each child becomes a
        panel (its value is reused, not recomputed).  This is the classic
        robust test: a feature smaller than the panel shows up as a mismatch.

        scheme 'embedded': q_near against q_near - 2 on the same panel.
        Cheaper, but both rules can miss the same sharp feature, so it is
        kept only for the benchmark that shows it.

        The tolerance is on energy (tol_kBT) and is converted to the bare
        integral of each face pair through the prefactor and n_i . n_j.
        """
        s = self.settings
        F = len(cos)
        tol_J = s.tol_kBT * 1.380649e-23 * s.temperature_K
        tol_I = tol_J / (self._prefactor * np.maximum(np.abs(cos), 1e-12))
        result = np.zeros(F)
        fp = np.arange(F)
        uc = np.zeros(F)
        vc = np.zeros(F)
        hw = np.full(F, self.half)

        def split(fp, uc, vc, hw):
            uc = (uc[:, None] + self._child[None, :, 0] * hw[:, None]).ravel()
            vc = (vc[:, None] + self._child[None, :, 1] * hw[:, None]).ravel()
            return np.repeat(fp, 4), uc, vc, np.repeat(0.5 * hw, 4)

        for _ in range(s.min_level):
            fp, uc, vc, hw = split(fp, uc, vc, hw)
        if s.scheme == "embedded":
            for level in range(s.min_level, s.max_level + 1):
                hi, lo = self._panel_values(fp, uc, vc, hw, face)
                share = (hw / self.half) ** 2
                done = (np.abs(hi - lo) <= tol_I[fp] * share) | (level == s.max_level)
                np.add.at(result, fp[done], hi[done])
                if np.all(done):
                    break
                fp, uc, vc, hw = split(fp[~done], uc[~done], vc[~done], hw[~done])
            return result
        c1, a1, b1, c2, a2, b2, n2 = face

        def distance_lb(fp, uc, vc, hw):
            """Lower bound on the distance from each panel to its inner face."""
            centre = c1[fp] + uc[:, None] * a1[fp] + vc[:, None] * b1[fp]
            d = centre - c2[fp]
            du = np.maximum(np.abs(np.einsum("fk,fk->f", d, a2[fp])) - self.half, 0.0)
            dv = np.maximum(np.abs(np.einsum("fk,fk->f", d, b2[fp])) - self.half, 0.0)
            dz = np.einsum("fk,fk->f", d, n2[fp])
            return np.sqrt(du * du + dv * dv + dz * dz) - np.sqrt(2.0) * hw

        parent, _ = self._panel_values(fp, uc, vc, hw, face)
        for level in range(s.min_level, s.max_level + 1):
            cf, cu, cv, ch = split(fp, uc, vc, hw)
            kids = self._panel_values(cf, cu, cv, ch, face)[0].reshape(-1, 4)
            refined = kids.sum(axis=1)
            share = (hw / self.half) ** 2
            agree = np.abs(refined - parent) <= tol_I[fp] * share
            if s.scheme == "hybrid":
                rho = np.sqrt(2.0) * hw
                d_lb = distance_lb(fp, uc, vc, hw)
                small = rho <= s.eta_geometric * d_lb
                agree = small | (agree & (rho <= s.eta_value * d_lb))
            last = level == s.max_level
            if last:
                self.max_level_hits += int(np.sum(~agree))
            done = agree | last
            np.add.at(result, fp[done], refined[done])
            keep = np.repeat(~done, 4)
            if not np.any(keep):
                break
            fp, uc, vc, hw = cf[keep], cu[keep], cv[keep], ch[keep]
            parent = kids.ravel()[keep]
        return result

    # -- surface form --------------------------------------------------------
    def _surface(self, p1, Q1, p2, Q2) -> np.ndarray:
        s = self.settings
        c1, n1, a1, b1 = cube_faces(p1, Q1, self.half)
        c2, n2, a2, b2 = cube_faces(p2, Q2, self.half)
        P = len(c1)
        cos = np.einsum("pik,pjk->pij", n1, n2)                  # (P, 6, 6)
        gap = (np.linalg.norm(c1[:, :, None, :] - c2[:, None, :, :], axis=-1)
               - np.sqrt(2.0) * self.edge_m) / self.edge_m       # lower bound
        tier = np.where(gap < s.near_L, 0, np.where(gap < s.mid_L, 1, 2))
        total = np.zeros(P)
        p, i, j = np.nonzero((tier == 0) & (np.abs(cos) > 1e-15))
        if len(p):
            face = (c1[p, i], a1[p, i], b1[p, i], c2[p, j], a2[p, j], b2[p, j], n2[p, j])
            np.add.at(total, p, cos[p, i, j] * self._adaptive(face, cos[p, i, j]))
        for level, q in ((1, s.q_mid), (2, s.q_far)):
            p, i, j = np.nonzero((tier == level) & (np.abs(cos) > 1e-15))
            if len(p) == 0:
                continue
            u, v, w = self._rules[q]
            nodes = (c1[p, i][:, None, :] + u[None, :, None] * a1[p, i][:, None, :]
                     + v[None, :, None] * b1[p, i][:, None, :])       # (F, nq, 3)
            inner = self._inner(nodes, c2[p, j], a2[p, j], b2[p, j], n2[p, j])
            np.add.at(total, p, cos[p, i, j] * (inner @ w))
        return self._prefactor * total

    # -- volume form for distant pairs --------------------------------------
    def _volume(self, p1, Q1, p2, Q2) -> np.ndarray:
        g, w = self._vol_nodes, self._vol_weights
        x1 = np.einsum("pab,nb->pna", Q1, g) + p1[:, None, :]
        x2 = np.einsum("pab,nb->pna", Q2, g) + p2[:, None, :]
        d = x2[:, None, :, :] - x1[:, :, None, :]
        r2 = np.einsum("pijk,pijk->pij", d, d)
        return -(self.hamaker_J / np.pi ** 2) * np.einsum("i,pij,j->p", w, r2 ** -3, w)

    def pair_energies(self, p1, Q1, p2, Q2) -> np.ndarray:
        """Energies in J for P pairs: p1, p2 (P, 3) and Q1, Q2 (P, 3, 3)."""
        p1 = np.asarray(p1, float).reshape(-1, 3)
        p2 = np.asarray(p2, float).reshape(-1, 3)
        Q1 = np.asarray(Q1, float).reshape(-1, 3, 3)
        Q2 = np.asarray(Q2, float).reshape(-1, 3, 3)
        if len(p1) != len(p2):
            p1, p2 = np.broadcast_arrays(p1, p2)
            Q1, Q2 = np.broadcast_arrays(Q1, Q2)
        out = np.empty(len(p1))
        far = (np.linalg.norm(p2 - p1, axis=1)
               > self.settings.far_pair_L * self.edge_m)
        if np.any(far):
            out[far] = self._volume(p1[far], Q1[far], p2[far], Q2[far])
        if np.any(~far):
            near = ~far
            out[near] = self._surface(p1[near], Q1[near], p2[near], Q2[near])
        return out

    def pair_energy(self, p1, Q1, p2, Q2) -> float:
        return float(self.pair_energies(p1, Q1, p2, Q2)[0])


# --------------------------------------------------------------------------
# independent check: octree volume integration with a purely geometric rule
# --------------------------------------------------------------------------
_CHILD3 = np.stack(np.meshgrid([-0.25, 0.25], [-0.25, 0.25], [-0.25, 0.25],
                               indexing="ij"), -1).reshape(-1, 3)


def octree_energy(p1, Q1, p2, Q2, eta=0.35, order=3, edge_m=EDGE_M,
                  hamaker_J=HAMAKER_J, max_level=16, max_pairs=4_000_000):
    """Volume integral by cell-pair refinement; validation only.

    Both cubes are split into equal cells; a cell pair is integrated with an
    order^3 x order^3 Gauss rule once cell diameter <= eta x (its lower-bound
    separation), otherwise both cells are split into eight.  The rule never
    looks at integrand values, so it cannot be fooled by a narrow spike, and
    it shares no code with the surface form.
    """
    x, w = np.polynomial.legendre.leggauss(order)
    g = np.stack(np.meshgrid(x, x, x, indexing="ij"), -1).reshape(-1, 3) * 0.5
    wg = np.einsum("i,j,k->ijk", w, w, w).ravel() / 8.0
    ca, cb = np.zeros((1, 3)), np.zeros((1, 3))
    h = 1.0
    total = 0.0
    for level in range(max_level + 1):
        A = np.asarray(p1) + (ca * edge_m) @ np.asarray(Q1).T
        B = np.asarray(p2) + (cb * edge_m) @ np.asarray(Q2).T
        D = np.linalg.norm(B - A, axis=1)
        diam = np.sqrt(3.0) * h * edge_m
        d_lb = D - diam
        acc = (d_lb > 0) & (diam <= eta * d_lb)
        if level == max_level:
            if np.any(~acc):
                raise RuntimeError("octree did not converge")
        for k in np.array_split(np.nonzero(acc)[0], max(1, int(np.sum(acc)) // 20000 + 1)):
            if len(k) == 0:
                continue
            ga = (g * h * edge_m) @ np.asarray(Q1).T
            gb = (g * h * edge_m) @ np.asarray(Q2).T
            PA = A[k][:, None, :] + ga[None]
            PB = B[k][:, None, :] + gb[None]
            d = PB[:, None, :, :] - PA[:, :, None, :]
            r2 = np.einsum("nijk,nijk->nij", d, d)
            total += float(np.einsum("i,nij,j->", wg, r2 ** -3, wg)) * (h * edge_m) ** 6
        rest = ~acc
        if not np.any(rest):
            break
        ra, rb = ca[rest], cb[rest]
        if 64 * len(ra) > max_pairs:
            raise RuntimeError("octree refinement too large")
        h *= 0.5
        na = (ra[:, None, :] + _CHILD3[None] * 2 * h).reshape(-1, 3)
        nb = (rb[:, None, :] + _CHILD3[None] * 2 * h).reshape(-1, 8, 3)
        ca = np.repeat(na, 8, axis=0)
        cb = np.repeat(nb, 8, axis=0).reshape(-1, 3)      # (parent, i, j) -> B child j
    return -(hamaker_J / np.pi ** 2) * total
