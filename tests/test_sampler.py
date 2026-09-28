"""Each move samples the right distribution, checked against an independent answer.

A Markov chain can be internally consistent and still sample the wrong
distribution (a missing Jacobian, an asymmetric proposal).  So every kernel
is run on a problem whose answer is known by other means -- a closed form,
a quadrature, or i.i.d. importance sampling -- and the MCMC mean must agree
within 4 Monte Carlo standard errors (batch means).
"""
import unittest

import numpy as np

from v904.diagnostics import batch_means_mcse, integrated_time
from v904.model import ClusterEnsemble, Hamiltonian, PhysicalParameters
from v904.sampler import Chain, MoveSettings

SIGMA = 4.0


def one_particle(**phys):
    ham = Hamiltonian(PhysicalParameters(**phys))
    state = dict(pos=np.zeros((1, 3)), Q=np.eye(3)[None], mu=np.array([[0.0, 0.0, 1.0]]))
    return ham, state


class Diagnostics(unittest.TestCase):
    def test_ar1_tau_and_mcse(self):
        """AR(1) with phi: tau = (1 + phi) / (1 - phi) exactly."""
        rng = np.random.default_rng(0)
        phi, n = 0.8, 200_000
        x = np.empty(n)
        x[0] = 0
        e = rng.normal(size=n)
        for t in range(1, n):
            x[t] = phi * x[t - 1] + e[t]
        tau, _, reliable = integrated_time(x)
        self.assertTrue(reliable)
        self.assertAlmostEqual(tau / 9.0, 1.0, delta=0.08)
        exact = np.sqrt(tau / n) * np.std(x)
        self.assertAlmostEqual(batch_means_mcse(x) / exact, 1.0, delta=0.15)


class RotationKernel(unittest.TestCase):
    def test_haar_uniform_when_energy_is_flat(self):
        """K = 0, no field on the body: rotations must be Haar, E[cos w] = -1/2."""
        ham, state = one_particle(K_Jpm3=0.0)
        chain = Chain(ham, state, MoveSettings(rot_step_deg=90.0), np.random.default_rng(1))
        cos_w = np.empty(40_000)
        for t in range(len(cos_w)):
            chain.rotate(0)
            cos_w[t] = 0.5 * (np.trace(chain.Q[0]) - 1.0)
        self.assertLess(abs(cos_w.mean() + 0.5), SIGMA * batch_means_mcse(cos_w))


class DipoleKernel(unittest.TestCase):
    def test_single_moment_against_sphere_quadrature(self):
        ham, state = one_particle()
        chain = Chain(ham, state, MoveSettings(dip_step_rad=1.6), np.random.default_rng(2))
        mx = np.empty(60_000)
        for t in range(len(mx)):
            chain.dipole(0)
            mx[t] = chain.mu[0, 0]
        # quadrature: Gauss-Legendre in cos(theta) x uniform phi
        c, w = np.polynomial.legendre.leggauss(200)
        phi = np.linspace(0, 2 * np.pi, 400, endpoint=False)
        C, P = np.meshgrid(c, phi, indexing="ij")
        S = np.sqrt(1 - C ** 2)
        m = np.stack([S * np.cos(P), S * np.sin(P), C], -1).reshape(-1, 3)
        z, a = ham.single_terms(np.repeat(np.eye(3)[None], len(m), 0), m)
        weight = np.repeat(w, len(phi)) * np.exp(-(z + a - (z + a).min()) / ham.kT)
        exact = float(np.sum(weight * m[:, 0]) / np.sum(weight))
        self.assertLess(abs(mx.mean() - exact), SIGMA * batch_means_mcse(mx))


def two_cubes(bond_nm=22.0, **phys):
    ham = Hamiltonian(PhysicalParameters(**phys), ClusterEnsemble(bond_nm=bond_nm))
    state = dict(pos=np.array([[0, 0, 0], [19.5e-9, 0, 0]]),
                 Q=np.repeat(np.eye(3)[None], 2, 0),
                 mu=np.array([[1.0, 0, 0], [0.6, 0.8, 0]]))
    return ham, state


def pair_energy(ham, state, p):
    """U(r) of the frozen pair for many positions p (n, 3) of cube 1."""
    n = len(p)
    eye = np.repeat(np.eye(3)[None], n, 0)
    w = ham.vdw_pairs(np.zeros((n, 3)), eye, p, eye)
    s, _ = ham.steric_pairs(p, eye, eye)
    m0, m1 = state["mu"]
    d = np.linalg.norm(p, axis=1)
    rh = p / d[:, None]
    return w + s + ham.dd_coeff * (m0 @ m1 - 3 * (rh @ m0) * (rh @ m1)) / d ** 3


class TranslationKernel(unittest.TestCase):
    """Two cubes, bodies and moments frozen; vdW + dipole + steric + constraints.

    Interaction strengths are softened (A, Ms smaller) so uniform i.i.d.
    importance sampling is an accurate reference; the kernel does not care
    how deep the wells are.
    """

    @classmethod
    def setUpClass(cls):
        # a wider bond ball, so fewer uniform draws land on overlapping cores
        cls.ham, cls.state = two_cubes(bond_nm=26.0, hamaker_J=4e-21, Ms_Apm=1.2e5)
        rng = np.random.default_rng(3)
        n, R = 60_000, cls.ham.bond_m
        v = rng.normal(size=(n, 3))
        v *= (R * rng.random(n) ** (1 / 3) / np.linalg.norm(v, axis=1))[:, None]
        # the reference must sample the SAME target: drop every state the
        # ensemble excludes (overlap or sharp cores closer than D0), exactly
        # as the chain's constraint check does
        v = v[[not cls.ham.overlap_any(np.zeros(3), np.eye(3), p, np.eye(3)) for p in v]]
        U = np.concatenate([pair_energy(cls.ham, cls.state, c) for c in np.array_split(v, 60)])
        wgt = np.exp(-(U - U.min()) / cls.ham.kT)
        r = np.linalg.norm(v, axis=1)
        cls.exact_r = float(np.sum(wgt * r) / np.sum(wgt))
        cls.ess_is = float(wgt.sum() ** 2 / np.sum(wgt ** 2))
        # self-normalised IS standard error of the weighted mean
        cls.se_is = float(np.sqrt(np.sum(wgt ** 2 * (r - cls.exact_r) ** 2)) / wgt.sum())

    def sample(self, use_scale, steps=10_000):
        chain = Chain(self.ham, self.state, MoveSettings(trans_step_nm=0.8, scale_log_step=0.03),
                      np.random.default_rng(4))
        r = np.empty(steps)
        for t in range(steps):
            chain.translate(1)
            if use_scale:
                chain.scale()
            r[t] = np.linalg.norm(chain.pos[1])
        self.assertLess(abs(chain.drift_kBT()), 1e-6)
        return r[steps // 10:]

    def test_reference_is_usable(self):
        self.assertGreater(self.ess_is, 1500)

    def test_translation_matches_importance_sampling(self):
        r = self.sample(False)
        err = np.hypot(batch_means_mcse(r), self.se_is)
        self.assertLess(abs(r.mean() - self.exact_r), SIGMA * err)

    def test_translation_plus_scale_matches(self):
        r = self.sample(True)
        err = np.hypot(batch_means_mcse(r), self.se_is)
        self.assertLess(abs(r.mean() - self.exact_r), SIGMA * err)


class ScaleJacobian(unittest.TestCase):
    """Scale moves alone slide cube 1 along a fixed ray; exact 1D answer.

    Along the ray the target is p(s) ~ s^2 exp(-U(s e)/kT) (N = 2: one free
    particle, Jacobian s^{3(N-1)} = s^3 per log s, i.e. s^2 ds).  Softened
    wells, so the distribution spreads over the ray and the s^2 factor
    matters: with the real 7 kBT contact well the pair sits at contact and a
    missing Jacobian would shift the mean by < 1 MCSE.
    """

    @classmethod
    def setUpClass(cls):
        cls.ham, cls.state = two_cubes(hamaker_J=4e-21, Ms_Apm=1.2e5)
        e = np.array([1.0, 0.0, 0.0])
        # start 1 nm past contact: the steric wall already makes the weight
        # exp(-8e9); at zero gap the integrand is not defined at all
        s = np.linspace(17e-9, cls.ham.bond_m, 4001)
        U = pair_energy(cls.ham, cls.state, s[:, None] * e)
        cls.with_jac = cls._mean(s, s ** 2 * np.exp(-(U - U.min()) / cls.ham.kT))
        cls.without = cls._mean(s, np.exp(-(U - U.min()) / cls.ham.kT))

    @staticmethod
    def _mean(s, w):
        return float(np.trapezoid(s * w, s) / np.trapezoid(w, s))

    def sample(self, jacobian):
        chain = Chain(self.ham, self.state, MoveSettings(scale_log_step=0.08),
                      np.random.default_rng(6))
        chain.scale_jacobian = jacobian
        r = np.empty(40_000)
        for t in range(len(r)):
            chain.scale()
            r[t] = chain.pos[1, 0]
        return r[3000:]

    def test_with_jacobian_matches_quadrature(self):
        r = self.sample(True)
        self.assertLess(abs(r.mean() - self.with_jac), SIGMA * batch_means_mcse(r))

    def test_missing_jacobian_is_detected(self):
        r = self.sample(False)
        self.assertGreater(abs(r.mean() - self.with_jac), 5 * batch_means_mcse(r))
        self.assertLess(abs(r.mean() - self.without), SIGMA * batch_means_mcse(r))


if __name__ == "__main__":
    unittest.main()
