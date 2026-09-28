"""V9.04 Hamiltonian: V9.03 physics with the converged sharp-cube vdW.

Only the vdW term changes.  Zeeman, first-order cubic anisotropy
(`cubic_first_raw`), point-dipole coupling and the rounded-support steric
spring are the V9.03 formulas with the V9.03 parameters, and
`tests/test_model.py` checks them term by term against the V9.03 notebook
cells.

What is new is that the target distribution is now a proper probability
distribution.  V9.03 documented (docs/MCMC_DIAGNOSTICS_ZH.md section 9) that
its position space is unbounded: a particle carried far away keeps a finite
energy, so int exp(-U/kT) d^3r diverges and no MCMC can converge to it.
`ClusterEnsemble` fixes that with the Stillinger cluster definition:

    pi(x) ~ exp(-U(x)/kT) * 1[bond graph connected] * 1[sharp cores >= D0 apart]

with a bond meaning centre distance < bond_nm.  Connectivity bounds every
distance by 26 * bond_nm, so the support is compact and pi normalisable.
The core indicator is needed twice over: the vdW integral over two
interpenetrating cubes diverges, and even for disjoint sharp cores it can
grow like 1/delta near contact (see ClusterEnsemble.core_min_nm).
Both constraints are hard walls, handled by rejecting a proposal that
violates them; how often that happens is reported as a diagnostic.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict

import numpy as np
from scipy import constants
from scipy.spatial.transform import Rotation

from . import geometry
from .vdw import HamakerEvaluator, HamakerSettings, uniform_voxel_energy

K_B = constants.k
MU0_4PI = constants.mu_0 / (4.0 * np.pi)
E111 = np.ones(3) / np.sqrt(3.0)
FIELD_HAT = np.array([1.0, 0.0, 0.0])


@dataclass(frozen=True)
class PhysicalParameters:
    """V9.03 PHYS_500G at 298.15 K; only `vdw_model` is new."""
    edge_nm: float = 16.0
    Ms_Apm: float = 285_000.0
    B_T: float = 0.05
    K_Jpm3: float = 2.0e4
    hamaker_J: float = 2.0e-20
    gap_nm: float = 1.8
    k_stiff_Jpm2: float = 1.0e8
    roundness_nm: float = 1.5
    temperature_K: float = 298.15
    # 'converged' = V10 sharp-cube integral at its voxel limit (default);
    # 'inherited' = V9.03 4^3 voxels, d^2 floor 1e-19, cutoff 2.2 a (regression only)
    vdw_model: str = "converged"
    inherited_cutoff_nm: float = 2.2 * 21.0
    vdw_settings: HamakerSettings = field(default_factory=HamakerSettings)


@dataclass(frozen=True)
class ClusterEnsemble:
    """Constraints that make the target normalisable.

    core_min_nm: sharp cores must stay at least this far apart.  0.165 nm is
    the standard Hamaker contact cutoff D0 (Israelachvili).  It is required,
    not cosmetic: the rounded-support steric wall lets the SHARP corners the
    vdW integral runs over approach each other (they protrude up to
    (sqrt3 - 1) x 1.5 = 1.1 nm beyond the rounded body), and where two nearly
    parallel faces meet near such a contact the sharp-cube attraction grows
    like 1/delta, so int exp(+c/delta) d delta diverges.  A production chain
    found exactly such a state (two cores 1.5e-5 nm apart, -458 kBT; an
    independent octree integration confirmed the value is the true integral).
    """
    bond_nm: float = 30.0
    core_min_nm: float = 0.165
    enforce_connectivity: bool = True
    enforce_core_exclusion: bool = True


class Hamiltonian:
    def __init__(self, phys: PhysicalParameters = PhysicalParameters(),
                 ensemble: ClusterEnsemble = ClusterEnsemble()):
        self.phys = phys
        self.ensemble = ensemble
        self.L = phys.edge_nm * 1e-9
        self.vol = self.L ** 3
        self.mu_mag = phys.Ms_Apm * self.vol
        self.B = phys.B_T
        self.kT = K_B * phys.temperature_K
        self.gap0 = phys.gap_nm * 1e-9
        self.r_round = min(phys.roundness_nm * 1e-9, self.L / 2 - 1e-12)
        self.dd_coeff = MU0_4PI * self.mu_mag ** 2
        self.bond_m = ensemble.bond_nm * 1e-9
        if phys.vdw_model not in ("converged", "inherited"):
            raise ValueError("vdw_model must be 'converged' or 'inherited'")
        self.hamaker = HamakerEvaluator(self.L, phys.hamaker_J, phys.vdw_settings)
        self.exact_distance_checks = 0

    # ------------------------------------------------------------------ terms
    def single_terms(self, Q, mu):
        """Zeeman and cubic anisotropy per particle, (N,) each, in J."""
        Q = np.asarray(Q).reshape(-1, 3, 3)
        mu = np.asarray(mu).reshape(-1, 3)
        zeeman = -self.mu_mag * self.B * mu[:, 0]
        mb = np.einsum("nab,na->nb", Q, mu)                 # Q^T mu
        mx, my, mz = mb.T
        aniso = -self.phys.K_Jpm3 * self.vol * ((mx * my) ** 2 + (mx * mz) ** 2
                                                + (my * mz) ** 2)
        return zeeman, aniso

    def dipole_row(self, i, pos, mu):
        """U_dd(i, j) for every j (0 at j = i), in J."""
        r = pos - pos[i]
        d2 = np.einsum("ij,ij->i", r, r)
        d2[i] = 1.0
        d = np.sqrt(d2)
        rhat = r / d[:, None]
        val = (mu @ mu[i] - 3.0 * (rhat @ mu[i]) * np.einsum("ij,ij->i", rhat, mu)) / d ** 3
        val[i] = 0.0
        return self.dd_coeff * val

    def dipole_total(self, pos, mu):
        r = pos[None, :, :] - pos[:, None, :]
        d2 = np.einsum("ijk,ijk->ij", r, r)
        np.fill_diagonal(d2, 1.0)
        d = np.sqrt(d2)
        rhat = r / d[..., None]
        mi_r = np.einsum("ijk,ik->ij", rhat, mu)
        mj_r = np.einsum("ijk,jk->ij", rhat, mu)
        val = (mu @ mu.T - 3.0 * mi_r * mj_r) / d ** 3
        np.fill_diagonal(val, 0.0)
        return 0.5 * self.dd_coeff * float(val.sum())

    def steric_pairs(self, r_vec, Qa, Qb):
        gap = geometry.support_gap(r_vec, Qa, Qb, self.L, self.r_round)
        over = np.maximum(self.gap0 - gap, 0.0)
        return 0.5 * self.phys.k_stiff_Jpm2 * over ** 2, gap

    def vdw_pairs(self, pa, Qa, pb, Qb):
        if self.phys.vdw_model == "converged":
            return self.hamaker.pair_energies(pa, Qa, pb, Qb)
        out = np.zeros(len(np.atleast_2d(pa)))
        pa, pb = np.atleast_2d(pa), np.atleast_2d(pb)
        Qa, Qb = np.reshape(Qa, (-1, 3, 3)), np.reshape(Qb, (-1, 3, 3))
        for k in range(len(out)):
            if np.linalg.norm(pb[k] - pa[k]) < self.phys.inherited_cutoff_nm * 1e-9:
                out[k] = uniform_voxel_energy(pa[k], Qa[k], pb[k], Qb[k], 4, self.L,
                                              self.phys.hamaker_J, d2_floor_m2=1e-19)
        return out

    # ------------------------------------------------------------ constraints
    def overlap_any(self, pa, Qa, pb, Qb):
        """True if any pair of sharp cores is closer than core_min (or overlaps).

        The separating-axis test clears almost every pair at once: an axis
        separating the cores by more than core_min proves the distance is at
        least that.  Only pairs it cannot clear get the exact distance.
        """
        if not self.ensemble.enforce_core_exclusion:
            return False
        r = np.atleast_2d(pb) - np.atleast_2d(pa)
        d_min = self.ensemble.core_min_nm * 1e-9
        close = np.linalg.norm(r, axis=1) < np.sqrt(3.0) * self.L + d_min
        if not np.any(close):
            return False
        Qa = np.broadcast_to(np.reshape(Qa, (-1, 3, 3)), (len(r), 3, 3))
        Qb = np.broadcast_to(np.reshape(Qb, (-1, 3, 3)), (len(r), 3, 3))
        idx = np.nonzero(close)[0]
        unclear = geometry.cores_overlap(r[idx], Qa[idx], Qb[idx], self.L, tol_m=d_min)
        for k in idx[unclear]:
            self.exact_distance_checks += 1
            if geometry.core_distance(np.zeros(3), Qa[k], r[k], Qb[k], self.L) < d_min:
                return True
        return False

    def allowed(self, pos, Q):
        """Full constraint check for a whole configuration."""
        if self.ensemble.enforce_connectivity and not geometry.connected(pos, self.bond_m):
            return False
        i, j = np.triu_indices(len(pos), 1)
        return not self.overlap_any(pos[i], Q[i], pos[j], Q[j])

    # ------------------------------------------------------------ full energy
    def pair_matrices(self, pos, Q):
        """Symmetric (N, N) vdW and steric pair energies, plus the gap matrix."""
        N = len(pos)
        i, j = np.triu_indices(N, 1)
        W = np.zeros((N, N))
        S = np.zeros((N, N))
        G = np.full((N, N), np.inf)
        w = self.vdw_pairs(pos[i], Q[i], pos[j], Q[j])
        s, g = self.steric_pairs(pos[j] - pos[i], Q[i], Q[j])
        W[i, j] = W[j, i] = w
        S[i, j] = S[j, i] = s
        G[i, j] = G[j, i] = g
        return W, S, G

    def energy_terms(self, pos, Q, mu):
        z, a = self.single_terms(Q, mu)
        W, S, _ = self.pair_matrices(pos, Q)
        terms = dict(Zeeman=float(z.sum()), Anisotropy=float(a.sum()),
                     Dipole=self.dipole_total(pos, mu),
                     VdW=0.5 * float(W.sum()), Steric=0.5 * float(S.sum()))
        terms["Total"] = sum(terms.values())
        return terms

    def describe(self):
        return dict(physical=asdict(self.phys), ensemble=asdict(self.ensemble),
                    kT_J=self.kT, mu_Am2=self.mu_mag)


# ---------------------------------------------------------------------------
# initial structures: the V9.03 builders, reproduced exactly
# ---------------------------------------------------------------------------
def rhombo_basis(a_nm, alpha_deg):
    a = a_nm * 1e-9
    al = np.radians(alpha_deg)
    cp = np.sqrt(max((2.0 * np.cos(al) + 1.0) / 3.0, 0.0))
    sp = np.sqrt(max(1.0 - cp ** 2, 0.0))
    # 0.866 (not sqrt(3)/2) is what V9.03 uses; kept for identical starts
    return np.array([[cp, sp, 0.0], [cp, -0.5 * sp, 0.866 * sp],
                     [cp, -0.5 * sp, -0.866 * sp]]) * a


# ideal body orientation in the lattice frame: body [111] along the 3-fold axis
LATTICE_BODY_REF = Rotation.align_vectors([[1.0, 0.0, 0.0]], [E111])[0].as_matrix()


def initial_state(a_nm=21.0, alpha_deg=74.2, tilt_deg=0.0, spacing_scale=1.0,
                  dipole_seed=None, dipoles="random", azimuth_deg=0.0):
    """V9.03 `make_initial_state`: 3x3x3 rhombohedral, bodies [111] || lattice axis.

    Also returns G, the lattice orientation (lattice frame -> lab), which the
    rigid-lattice runs track.  `azimuth_deg` then spins the co-rotated
    structure about the field axis (default 0 = V9.03).
    """
    Ry = Rotation.from_euler("y", tilt_deg, degrees=True).as_matrix()
    if azimuth_deg:
        Ry = Rotation.from_euler("x", azimuth_deg, degrees=True).as_matrix() @ Ry
    basis = rhombo_basis(a_nm, alpha_deg) @ Ry.T
    r = np.arange(-1, 2)
    idx = np.array(np.meshgrid(r, r, r, indexing="ij")).T.reshape(-1, 3)
    pos = idx @ basis
    dist = np.linalg.norm(pos, axis=1) * 1e9
    center = int(np.where(np.all(idx == 0, axis=1))[0][0])
    order = [center] + [k for k in np.argsort(dist) if k != center]
    idx, pos = idx[order], pos[order]
    Q = np.repeat((Ry @ LATTICE_BODY_REF)[None], len(idx), axis=0)
    if dipoles == "random":
        v = np.random.default_rng(dipole_seed).normal(size=(len(idx), 3))
        mu = v / np.linalg.norm(v, axis=1, keepdims=True)
    elif dipoles == "field":
        mu = np.tile(FIELD_HAT, (len(idx), 1))
    elif dipoles == "easy":
        mu = np.einsum("nij,j->ni", Q, E111)
    else:
        raise ValueError("dipoles must be random, field or easy")
    pos = pos[0] + spacing_scale * (pos - pos[0])
    return dict(idx=idx, pos=pos, Q=Q, mu=mu, G=Ry)


INITIAL_STRUCTURES = {
    "compact_aligned": dict(tilt_deg=0.0, spacing_scale=1.0),
    "compact_tilted40": dict(tilt_deg=40.0, spacing_scale=1.0),
    "expanded_aligned": dict(tilt_deg=0.0, spacing_scale=1.08),
    "expanded_tilted": dict(tilt_deg=40.0, spacing_scale=1.08),
}
