"""Metropolis-Hastings sampler for the V9.04 cluster ensemble.

Every move is a separate MH kernel that leaves pi invariant on its own:

  move           proposal                                  log q(x|x')/q(x'|x)
  -------------  ----------------------------------------  ------------------
  translate i    r_i + U(-s, s)^3              (i != 0)    0  (symmetric)
  rotate i       exp(theta [n]) Q_i, n ~ S^2, theta~U(-s,s) 0  (class function)
  dipole i       same small rotation applied to m_i        0
  co-tilt        rigid rotation of all r and Q about r_0   0
  twist (gamma)  all Q about their common [111] axis       0  (axis is kept)
  scale          r_0 + e^eps (r - r_0), eps ~ U(-s, s)     3 (N - 1) eps
  lattice rot.   all r about r_0, bodies fixed             0  (rigid-lattice runs)

Two schedules use these kernels:

* `cycle`         the cluster ensemble (positions free, Stillinger bonds);
* `cycle_lattice` the rigid-lattice ensemble: positions never move relative
                  to the lattice, which only rotates as a whole (co-tilt with
                  the bodies, or lattice rotation under fixed bodies).  The
                  state is (G, Q_1..Q_27, m_1..m_27) on SO(3) x SO(3)^27 x
                  (S^2)^27, compact, so pi needs no bond constraint.  This is
                  the ensemble for "given that the SL exists": the
                  colloidal/dispersed states are excluded by construction.

`temperature_factor` multiplies kT in every acceptance test.  It is 1 for
sampling; simulated annealing (v904/anneal.py) lowers it to use the same
kernels as an optimiser.

A cycle applies them in a fixed order, particles in a fresh random order
each sweep.  A composition of pi-invariant kernels is pi-invariant, and a
state-independent random order is a mixture of such compositions, so the
cycle kernel is too (it is not reversible as a whole, which MCMC does not
need).  Proposals violating a hard constraint have pi = 0 and are rejected
before any energy is computed.

Translation and body rotation may use a two-scale mixture: each attempt
takes the small step with probability `small_prob`, else the large one.  The
choice does not depend on the state, so the mixture of two symmetric
proposals is itself symmetric and the acceptance rule is unchanged.  It
hedges the tuning: large steps move loose surface particles, small steps
keep working if the cluster compacts.

Step sizes are FIXED inside a chain.  They are tuned only in a separate
pilot (scripts/pilot_tuning.py), whose samples are discarded, because
adapting on the running chain breaks the Markov property that the
convergence theory needs.  (Annealing is an optimiser, not a sampler, so it
may adapt; see v904/anneal.py.)
"""
from __future__ import annotations

import time
from dataclasses import dataclass, asdict

import numpy as np
from scipy.spatial.transform import Rotation

from . import geometry
from .model import E111, FIELD_HAT, LATTICE_BODY_REF, Hamiltonian


@dataclass(frozen=True)
class MoveSettings:
    trans_step_nm: float = 0.10
    rot_step_deg: float = 2.0
    trans_small_nm: float = 0.0     # used with probability small_prob
    rot_small_deg: float = 0.0
    small_prob: float = 0.0
    dip_step_rad: float = 1.60      # V9.03 fixed-geometry calibration, see README
    n_mag_sweeps: int = 10
    cotilt_step_deg: float = 3.0
    n_cotilt: int = 1
    gamma_step_deg: float = 2.0
    n_gamma: int = 1
    scale_log_step: float = 0.002
    n_scale: int = 1
    latrot_step_deg: float = 1.0    # rigid-lattice schedule only
    n_latrot: int = 1
    # True: re-evaluate vdW on co-tilt (bitwise-identical cache; used for the
    # cluster production run).  False: reuse it -- the quadrature is
    # rotation invariant up to round-off, which the drift check verifies.
    cotilt_recompute: bool = True


def random_unit(rng, n=None):
    v = rng.normal(size=(3,) if n is None else (n, 3))
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def small_rotation(rng, max_angle_rad):
    """exp(theta [n]); n uniform on S^2, theta uniform on [-max, max]."""
    return Rotation.from_rotvec(rng.uniform(-max_angle_rad, max_angle_rad)
                                * random_unit(rng)).as_matrix()


def mh_accept(log_alpha, rng):
    if np.isnan(log_alpha):
        raise ValueError("undefined acceptance ratio")
    return log_alpha >= 0.0 or np.log(rng.random()) < log_alpha


class Chain:
    """Chain state with cached per-particle and per-pair energies (J)."""

    def __init__(self, ham: Hamiltonian, state: dict, moves: MoveSettings, rng):
        self.h, self.moves, self.rng = ham, moves, rng
        self.pos = np.array(state["pos"], float)
        self.Q = np.array(state["Q"], float)
        self.mu = np.array(state["mu"], float)
        self.G = np.array(state.get("G", np.eye(3)), float)   # lattice orientation
        self.N = len(self.pos)
        self.temperature_factor = 1.0
        if not ham.allowed(self.pos, self.Q):
            raise ValueError("initial state violates the ensemble constraints")
        self.recompute()
        names = ["trans", "rot", "dip", "cotilt", "gamma", "scale", "latrot"]
        self.tries = dict.fromkeys(names, 0)
        self.accepts = dict.fromkeys(names, 0)
        self.blocked = dict.fromkeys(names, 0)     # rejected by a hard constraint
        self.sq_jump = dict.fromkeys(names, 0.0)   # accepted squared jump sizes
        # switch exists only so a test can show that dropping the Jacobian is caught
        self.scale_jacobian = True

    @property
    def kT_eff(self):
        return self.h.kT * self.temperature_factor

    # -- bookkeeping -----------------------------------------------------------
    def recompute(self):
        h = self.h
        self.zee, self.ani = h.single_terms(self.Q, self.mu)
        self.W, self.S, _ = h.pair_matrices(self.pos, self.Q)
        self.E_dd = h.dipole_total(self.pos, self.mu)

    def energy_terms(self):
        t = dict(Zeeman=float(self.zee.sum()), Anisotropy=float(self.ani.sum()),
                 Dipole=float(self.E_dd), VdW=0.5 * float(self.W.sum()),
                 Steric=0.5 * float(self.S.sum()))
        t["Total"] = sum(t.values())
        return t

    def drift_kBT(self):
        """Cached minus freshly computed total energy, in kBT."""
        fresh = self.h.energy_terms(self.pos, self.Q, self.mu)["Total"]
        return (self.energy_terms()["Total"] - fresh) / self.h.kT

    def _others(self, i):
        return np.arange(self.N) != i

    def _row(self, i, p_i, Q_i):
        """New vdW and steric rows of particle i (zero on the diagonal).

        Every pair is evaluated with the lower index as the first cube, as
        in `Hamiltonian.pair_matrices`.  The adaptive quadrature is exchange
        symmetric only to its tolerance, so a fixed order is what makes the
        cached energy the same deterministic function U_q(x) as a fresh one.
        """
        o = self._others(i)
        j = np.nonzero(o)[0]
        pos = self.pos.copy()
        Q = self.Q.copy()
        pos[i], Q[i] = p_i, Q_i
        a, b = np.minimum(i, j), np.maximum(i, j)
        w = np.zeros(self.N)
        s = np.zeros(self.N)
        w[o] = self.h.vdw_pairs(pos[a], Q[a], pos[b], Q[b])
        s[o], _ = self.h.steric_pairs(pos[b] - pos[a], Q[a], Q[b])
        return w, s

    def _set_row(self, i, w, s):
        self.W[i, :] = w
        self.W[:, i] = w
        self.S[i, :] = s
        self.S[:, i] = s

    # -- local moves ----------------------------------------------------------
    def translate(self, i):
        m, h, rng = self.moves, self.h, self.rng
        self.tries["trans"] += 1
        small = m.small_prob > 0 and rng.random() < m.small_prob
        step = m.trans_small_nm if small else m.trans_step_nm
        delta = rng.uniform(-1.0, 1.0, 3) * step * 1e-9
        new_p = self.pos[i] + delta
        o = self._others(i)
        trial = self.pos.copy()
        trial[i] = new_p
        if (h.ensemble.enforce_connectivity and not geometry.connected(trial, h.bond_m)) \
                or h.overlap_any(np.broadcast_to(new_p, (self.N - 1, 3)), self.Q[i],
                                 self.pos[o], self.Q[o]):
            self.blocked["trans"] += 1
            return False
        w, s = self._row(i, new_p, self.Q[i])
        dd_old = h.dipole_row(i, self.pos, self.mu)
        dd_new = h.dipole_row(i, trial, self.mu)
        dE = (w.sum() - self.W[i].sum()) + (s.sum() - self.S[i].sum()) \
            + (dd_new.sum() - dd_old.sum())
        if mh_accept(-dE / self.kT_eff, rng):
            self.pos[i] = new_p
            self._set_row(i, w, s)
            self.E_dd += dd_new.sum() - dd_old.sum()
            self.accepts["trans"] += 1
            self.sq_jump["trans"] += float(delta @ delta) * 1e18
            return True
        return False

    def rotate(self, i):
        m, h, rng = self.moves, self.h, self.rng
        self.tries["rot"] += 1
        small = m.small_prob > 0 and rng.random() < m.small_prob
        Rm = small_rotation(rng, np.radians(m.rot_small_deg if small else m.rot_step_deg))
        new_Q = Rm @ self.Q[i]
        o = self._others(i)
        if h.overlap_any(np.broadcast_to(self.pos[i], (self.N - 1, 3)), new_Q,
                         self.pos[o], self.Q[o]):
            self.blocked["rot"] += 1
            return False
        w, s = self._row(i, self.pos[i], new_Q)
        _, a_new = h.single_terms(new_Q, self.mu[i])
        dE = (w.sum() - self.W[i].sum()) + (s.sum() - self.S[i].sum()) + (a_new[0] - self.ani[i])
        if mh_accept(-dE / self.kT_eff, rng):
            self.Q[i] = new_Q
            self._set_row(i, w, s)
            self.ani[i] = a_new[0]
            self.accepts["rot"] += 1
            angle = np.linalg.norm(Rotation.from_matrix(Rm).as_rotvec())
            self.sq_jump["rot"] += float(np.degrees(angle)) ** 2
            return True
        return False

    def dipole(self, i):
        m, h, rng = self.moves, self.h, self.rng
        self.tries["dip"] += 1
        Rm = small_rotation(rng, m.dip_step_rad)
        new_mu = Rm @ self.mu[i]
        new_mu /= np.linalg.norm(new_mu)
        z_new, a_new = h.single_terms(self.Q[i], new_mu)
        dd_old = h.dipole_row(i, self.pos, self.mu).sum()
        old = self.mu[i].copy()
        self.mu[i] = new_mu
        dd_new = h.dipole_row(i, self.pos, self.mu).sum()
        dE = (z_new[0] - self.zee[i]) + (a_new[0] - self.ani[i]) + (dd_new - dd_old)
        if mh_accept(-dE / self.kT_eff, rng):
            self.zee[i], self.ani[i] = z_new[0], a_new[0]
            self.E_dd += dd_new - dd_old
            self.accepts["dip"] += 1
            self.sq_jump["dip"] += float(np.degrees(np.arccos(np.clip(old @ new_mu, -1, 1)))) ** 2
            return True
        self.mu[i] = old
        return False

    # -- collective moves -----------------------------------------------------
    def cotilt(self):
        """Rigid rotation of positions and bodies about r_0; moments stay.

        vdW, steric and both constraints are invariant under a rigid rotation
        of every position AND body together, so only anisotropy and dipole
        coupling change.  With `cotilt_recompute` the vdW is re-evaluated
        anyway, which keeps the cache bitwise equal to a fresh evaluation;
        without it the cached values are reused, exact up to round-off
        (pairs are always evaluated lower index first, so the only
        difference is floating-point noise, ~1e-13 kBT in the drift check).
        """
        m, h, rng = self.moves, self.h, self.rng
        self.tries["cotilt"] += 1
        Rg = small_rotation(rng, np.radians(m.cotilt_step_deg))
        angle_deg = float(np.degrees(np.linalg.norm(Rotation.from_matrix(Rg).as_rotvec())))
        anchor = self.pos[0]
        new_pos = anchor + (self.pos - anchor) @ Rg.T
        new_Q = np.einsum("ab,nbc->nac", Rg, self.Q)
        if m.cotilt_recompute:
            W, S, _ = h.pair_matrices(new_pos, new_Q)
        else:
            W, S = self.W, self.S
        _, a_new = h.single_terms(new_Q, self.mu)
        dd_new = h.dipole_total(new_pos, self.mu)
        dE = 0.5 * (W.sum() - self.W.sum() + S.sum() - self.S.sum())             + (a_new.sum() - self.ani.sum()) + (dd_new - self.E_dd)
        if mh_accept(-dE / self.kT_eff, rng):
            self.pos, self.Q, self.ani, self.E_dd = new_pos, new_Q, a_new, dd_new
            self.W, self.S = W, S
            self.G = Rg @ self.G
            self.accepts["cotilt"] += 1
            self.sq_jump["cotilt"] += angle_deg ** 2
            return True
        return False

    def lattice_rotate(self):
        """Rotate the lattice (all positions about r_0) under fixed bodies.

        Left multiplication of G by a small rotation: symmetric, and it
        preserves the Haar measure on G.  It changes how every cube sits in
        its cage, so vdW, steric and the core constraint are re-evaluated.
        """
        m, h, rng = self.moves, self.h, self.rng
        self.tries["latrot"] += 1
        Rg = small_rotation(rng, np.radians(m.latrot_step_deg))
        anchor = self.pos[0]
        new_pos = anchor + (self.pos - anchor) @ Rg.T
        i, j = np.triu_indices(self.N, 1)
        if h.overlap_any(new_pos[i], self.Q[i], new_pos[j], self.Q[j]):
            self.blocked["latrot"] += 1
            return False
        W, S, _ = h.pair_matrices(new_pos, self.Q)
        dd_new = h.dipole_total(new_pos, self.mu)
        dE = 0.5 * (W.sum() - self.W.sum() + S.sum() - self.S.sum()) + (dd_new - self.E_dd)
        if mh_accept(-dE / self.kT_eff, rng):
            self.pos, self.W, self.S, self.E_dd = new_pos, W, S, dd_new
            self.G = Rg @ self.G
            self.accepts["latrot"] += 1
            angle = np.degrees(np.linalg.norm(Rotation.from_matrix(Rg).as_rotvec()))
            self.sq_jump["latrot"] += float(angle) ** 2
            return True
        return False

    def gamma(self):
        """Rotate every body about the mean body-[111] axis.

        Rotating all bodies by R about axis a maps the mean easy axis to
        R a = a, so the reverse move uses the same axis with -theta: the
        proposal is symmetric even though the axis depends on the state.
        """
        m, h, rng = self.moves, self.h, self.rng
        self.tries["gamma"] += 1
        axis = np.einsum("nij,j->ni", self.Q, E111).mean(axis=0)
        n = np.linalg.norm(axis)
        axis = axis / n if n > 1e-12 else FIELD_HAT
        theta = rng.uniform(-1, 1) * np.radians(m.gamma_step_deg)
        Rg = Rotation.from_rotvec(theta * axis).as_matrix()
        new_Q = np.einsum("ab,nbc->nac", Rg, self.Q)
        i, j = np.triu_indices(self.N, 1)
        if h.overlap_any(self.pos[i], new_Q[i], self.pos[j], new_Q[j]):
            self.blocked["gamma"] += 1
            return False
        W, S, _ = h.pair_matrices(self.pos, new_Q)
        _, a_new = h.single_terms(new_Q, self.mu)
        dE = 0.5 * (W.sum() - self.W.sum() + S.sum() - self.S.sum()) + (a_new.sum() - self.ani.sum())
        if mh_accept(-dE / self.kT_eff, rng):
            self.Q, self.W, self.S, self.ani = new_Q, W, S, a_new
            self.accepts["gamma"] += 1
            self.sq_jump["gamma"] += float(np.degrees(theta)) ** 2
            return True
        return False

    def scale(self):
        """Log-uniform dilation about r_0 with its Jacobian e^{3(N-1) eps}."""
        m, h, rng = self.moves, self.h, self.rng
        self.tries["scale"] += 1
        eps = rng.uniform(-m.scale_log_step, m.scale_log_step)
        anchor = self.pos[0]
        new_pos = anchor + np.exp(eps) * (self.pos - anchor)
        i, j = np.triu_indices(self.N, 1)
        if (h.ensemble.enforce_connectivity and not geometry.connected(new_pos, h.bond_m)) \
                or h.overlap_any(new_pos[i], self.Q[i], new_pos[j], self.Q[j]):
            self.blocked["scale"] += 1
            return False
        W, S, _ = h.pair_matrices(new_pos, self.Q)
        dd_new = h.dipole_total(new_pos, self.mu)
        dE = 0.5 * (W.sum() - self.W.sum() + S.sum() - self.S.sum()) + (dd_new - self.E_dd)
        log_jac = 3 * (self.N - 1) * eps if self.scale_jacobian else 0.0
        if mh_accept(-dE / self.kT_eff + log_jac, rng):
            self.pos, self.W, self.S, self.E_dd = new_pos, W, S, dd_new
            self.accepts["scale"] += 1
            self.sq_jump["scale"] += float(eps) ** 2
            return True
        return False

    # -- one cycle ------------------------------------------------------------
    def cycle(self, move_positions=True, move_orientations=True, move_dipoles=True):
        m, rng = self.moves, self.rng
        for i in rng.permutation(self.N):
            if move_positions and i != 0:
                self.translate(i)
            if move_orientations:
                self.rotate(i)
        if move_positions and move_orientations:
            for _ in range(m.n_cotilt):
                self.cotilt()
        if move_orientations:
            for _ in range(m.n_gamma):
                self.gamma()
        if move_positions:
            for _ in range(m.n_scale):
                self.scale()
        if move_dipoles:
            for _ in range(m.n_mag_sweeps):
                for i in rng.permutation(self.N):
                    self.dipole(i)

    def cycle_lattice(self):
        """Rigid-lattice schedule: bodies, lattice orientation, moments."""
        m, rng = self.moves, self.rng
        for i in rng.permutation(self.N):
            self.rotate(i)
        for _ in range(m.n_latrot):
            self.lattice_rotate()
        for _ in range(m.n_cotilt):
            self.cotilt()
        for _ in range(m.n_gamma):
            self.gamma()
        for _ in range(m.n_mag_sweeps):
            for i in rng.permutation(self.N):
                self.dipole(i)

    def acceptance(self):
        return {k: self.accepts[k] / self.tries[k] if self.tries[k] else np.nan
                for k in self.tries}

    def reset_counters(self):
        for d in (self.tries, self.accepts, self.blocked):
            for k in d:
                d[k] = 0
        for k in self.sq_jump:
            self.sq_jump[k] = 0.0


# ---------------------------------------------------------------------------
# observables (V9.03 definitions)
# ---------------------------------------------------------------------------
def axial_tilt_deg(axis):
    n = np.linalg.norm(axis)
    if n < 1e-15:
        return np.nan
    return float(np.degrees(np.arccos(np.clip(abs(axis[0]) / n, 0.0, 1.0))))


def principal_axis(pos):
    c = pos - pos.mean(axis=0)
    vals, vecs = np.linalg.eigh(c.T @ c)
    order = np.argsort(vals)[::-1]
    vals, axis = vals[order], vecs[:, order[0]]
    if axis[0] < 0:
        axis = -axis
    denom = float(np.sum(np.maximum(vals, 0.0)))
    return axis, (0.0 if denom < 1e-30 else max(float((vals[0] - vals[1]) / denom), 0.0))


OBSERVABLES = ["energy_kBT", "vdw_kBT", "dipole_kBT", "local_beta_deg", "body_tilt_deg",
               "body_order", "sl_pca_tilt_deg", "sl_pca_order", "body_sl_difference_deg",
               "body_sl_axis_angle_deg", "magnetization", "abs_muB", "min_gap_nm",
               "p05_gap_nm", "rg_nm", "n_bonds"]
# rigid-lattice runs add these (G is meaningful only when positions move rigidly)
LATTICE_OBSERVABLES = ["sl_axis_tilt_deg", "body_rotation_deg", "body_rotation_max_deg",
                       "body_lattice_angle_deg", "anisotropy_kBT", "zeeman_kBT"]
LATTICE_AXIS = np.array([1.0, 0.0, 0.0])      # 3-fold axis of the rhombohedral basis


def observe(chain: Chain, lattice: bool = False) -> dict:
    h, pos, Q, mu = chain.h, chain.pos, chain.Q, chain.mu
    easy = np.einsum("nij,j->ni", Q, E111)
    beta = np.degrees(np.arccos(np.clip(np.sum(easy * mu, axis=1), -1, 1)))
    mean_axis = easy.mean(axis=0)
    order = float(np.linalg.norm(mean_axis))
    body_tilt = axial_tilt_deg(mean_axis)
    sl_axis, sl_order = principal_axis(pos)
    sl_tilt = axial_tilt_deg(sl_axis)
    axis_angle = float(np.degrees(np.arccos(np.clip(abs(mean_axis @ sl_axis) / max(order, 1e-15), 0, 1))))
    i, j = np.triu_indices(chain.N, 1)
    gaps = geometry.support_gap(pos[j] - pos[i], Q[i], Q[j], h.L, h.r_round) * 1e9
    if gaps.size == 0:                       # a single particle has no pairs
        gaps = np.array([np.nan])
    dist = np.linalg.norm(pos[j] - pos[i], axis=1)
    terms = chain.energy_terms()
    out = dict(
        energy_kBT=terms["Total"] / h.kT, vdw_kBT=terms["VdW"] / h.kT,
        dipole_kBT=terms["Dipole"] / h.kT, local_beta_deg=float(beta.mean()),
        body_tilt_deg=body_tilt, body_order=order, sl_pca_tilt_deg=sl_tilt,
        sl_pca_order=sl_order, body_sl_difference_deg=body_tilt - sl_tilt,
        body_sl_axis_angle_deg=axis_angle, magnetization=float(mu[:, 0].mean()),
        abs_muB=float(np.abs(mu[:, 0]).mean()), min_gap_nm=float(gaps.min()),
        p05_gap_nm=float(np.percentile(gaps, 5.0)),
        rg_nm=float(np.sqrt(np.mean(np.sum((pos - pos.mean(0)) ** 2, axis=1))) * 1e9),
        n_bonds=float(np.sum(dist < h.bond_m)))
    if lattice:
        axis = chain.G @ LATTICE_AXIS
        # each cube's rotation away from its ideal orientation in the lattice frame
        rel = np.einsum("ba,nbc,dc->nad", chain.G, Q, LATTICE_BODY_REF)
        omega = np.degrees(np.arccos(np.clip((np.trace(rel, axis1=1, axis2=2) - 1) / 2, -1, 1)))
        out.update(
            sl_axis_tilt_deg=axial_tilt_deg(axis),
            body_rotation_deg=float(omega.mean()), body_rotation_max_deg=float(omega.max()),
            body_lattice_angle_deg=float(np.degrees(np.arccos(np.clip(np.abs(easy @ axis), 0, 1))).mean()),
            anisotropy_kBT=terms["Anisotropy"] / h.kT, zeeman_kBT=terms["Zeeman"] / h.kT)
    return out


def run_chain(ham: Hamiltonian, state: dict, moves: MoveSettings, rng, n_cycles: int,
              progress=None, snapshot_every: int = 0, check_drift_every: int = 0,
              schedule: str = "cluster", temperature_factors=None):
    """Run n_cycles, recording every observable after every cycle.

    schedule 'cluster' uses Chain.cycle, 'lattice' Chain.cycle_lattice.
    temperature_factors (length n_cycles) turns the run into an annealing
    run; it is None for sampling.  Nothing is discarded here; warm-up is
    removed at analysis time so the burn-in choice can itself be checked.
    """
    if schedule not in ("cluster", "lattice"):
        raise ValueError("schedule must be 'cluster' or 'lattice'")
    lattice = schedule == "lattice"
    chain = Chain(ham, state, moves, rng)
    step = chain.cycle_lattice if lattice else chain.cycle
    keys = OBSERVABLES + (LATTICE_OBSERVABLES if lattice else [])
    traj = {k: np.empty(n_cycles) for k in keys}
    if temperature_factors is not None:
        traj["temperature_factor"] = np.asarray(temperature_factors, float)
    snaps, drifts = [], []
    t0 = time.perf_counter()
    for c in range(n_cycles):
        if temperature_factors is not None:
            chain.temperature_factor = float(temperature_factors[c])
        step()
        for k, v in observe(chain, lattice).items():
            traj[k][c] = v
        if snapshot_every and (c + 1) % snapshot_every == 0:
            snaps.append(dict(cycle=c + 1, pos=chain.pos.copy(), Q=chain.Q.copy(),
                              mu=chain.mu.copy(), G=chain.G.copy()))
        if check_drift_every and (c + 1) % check_drift_every == 0:
            drifts.append((c + 1, chain.drift_kBT()))
        if progress is not None:
            progress(c + 1, n_cycles, time.perf_counter() - t0)
    return dict(traj=traj, chain=chain, snapshots=snaps, drifts=drifts,
                elapsed_sec=time.perf_counter() - t0,
                acceptance=chain.acceptance(), tries=dict(chain.tries),
                blocked=dict(chain.blocked), sq_jump=dict(chain.sq_jump),
                moves=asdict(moves))
