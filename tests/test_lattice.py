"""Rigid-lattice ensemble and annealing, each checked against an independent answer."""
import unittest
from dataclasses import replace

import numpy as np

from v904.anneal import anneal, geometric_schedule
from v904.diagnostics import batch_means_mcse
from v904.model import ClusterEnsemble, Hamiltonian, PhysicalParameters, initial_state
from v904.sampler import Chain, MoveSettings, observe

SIGMA = 4.0


def dipole_pair(Ms=5.6e5, d_nm=40.0):
    """Two cubes 40 nm apart, moments frozen along x, no vdW, no field.

    U depends only on the pair axis u: U = C (1 - 3 u_x^2) / d^3, so the
    exact distribution of u is known up to a 1D integral over cos(theta).
    """
    ham = Hamiltonian(PhysicalParameters(Ms_Apm=Ms, hamaker_J=0.0, B_T=0.0, K_Jpm3=0.0),
                      ClusterEnsemble(enforce_connectivity=False))
    state = dict(pos=np.array([[0, 0, 0], [0, d_nm * 1e-9, 0]]),
                 Q=np.repeat(np.eye(3)[None], 2, 0), mu=np.array([[1.0, 0, 0], [1.0, 0, 0]]))
    c = np.linspace(-1, 1, 20001)
    U = ham.dd_coeff * (1 - 3 * c ** 2) / (d_nm * 1e-9) ** 3 / ham.kT
    w = np.exp(-(U - U.min()))
    exact = float(np.trapezoid(c ** 2 * w, c) / np.trapezoid(w, c))
    return ham, state, exact


class LatticeKernels(unittest.TestCase):
    def check(self, move_name):
        ham, state, exact = dipole_pair()
        moves = MoveSettings(latrot_step_deg=60.0, cotilt_step_deg=60.0, cotilt_recompute=False)
        chain = Chain(ham, state, moves, np.random.default_rng(11))
        move = getattr(chain, move_name)
        x = np.empty(30_000)
        for t in range(len(x)):
            move()
            u = chain.pos[1] - chain.pos[0]
            x[t] = (u[0] / np.linalg.norm(u)) ** 2
        x = x[1000:]
        self.assertLess(abs(x.mean() - exact), SIGMA * batch_means_mcse(x))
        # the pair distance never changes: the lattice is rigid
        self.assertAlmostEqual(np.linalg.norm(chain.pos[1] - chain.pos[0]) * 1e9, 40.0, places=9)
        return exact

    def test_lattice_rotation_samples_the_pair_axis_distribution(self):
        exact = self.check("lattice_rotate")
        self.assertGreater(abs(exact - 1 / 3), 0.05)     # the test is not trivial

    def test_cotilt_samples_the_same_distribution(self):
        self.check("cotilt")


class RigidBookkeeping(unittest.TestCase):
    def test_compact_lattice_cycles(self):
        ham = Hamiltonian(PhysicalParameters(), ClusterEnsemble(enforce_connectivity=False))
        s = initial_state(tilt_deg=20.0, dipole_seed=4)
        moves = MoveSettings(rot_step_deg=1.0, latrot_step_deg=0.5, cotilt_step_deg=6.0,
                             gamma_step_deg=1.0, n_mag_sweeps=2, cotilt_recompute=False)
        chain = Chain(ham, s, moves, np.random.default_rng(2))
        for _ in range(3):
            chain.cycle_lattice()
        self.assertTrue(all(chain.accepts[k] > 0 for k in ("rot", "cotilt", "dip")))
        self.assertLess(abs(chain.drift_kBT()), 1e-6)
        # positions are the reference lattice rotated by the tracked G
        ref = initial_state(tilt_deg=0.0, dipole_seed=4)
        where = {tuple(r): k for k, r in enumerate(ref["idx"])}
        order = [where[tuple(r)] for r in s["idx"]]
        lattice = ref["pos"][order] - ref["pos"][order][0]
        np.testing.assert_allclose(chain.pos - chain.pos[0], lattice @ chain.G.T, atol=1e-20)
        o = observe(chain, lattice=True)
        self.assertAlmostEqual(o["sl_axis_tilt_deg"], o["sl_pca_tilt_deg"], delta=1e-6)


class Annealing(unittest.TestCase):
    def test_single_moment_reaches_the_grid_minimum(self):
        ham = Hamiltonian(PhysicalParameters(), ClusterEnsemble(enforce_connectivity=False))
        state = dict(pos=np.zeros((1, 3)), Q=np.eye(3)[None], mu=np.array([[0.0, 0.0, -1.0]]))
        # exact minimum of Zeeman + cubic anisotropy on a dense sphere grid
        c, phi = np.meshgrid(np.linspace(-1, 1, 801), np.linspace(0, 2 * np.pi, 1601), indexing="ij")
        s = np.sqrt(1 - c ** 2)
        m = np.stack([s * np.cos(phi), s * np.sin(phi), c], -1).reshape(-1, 3)
        z, a = ham.single_terms(np.repeat(np.eye(3)[None], len(m), 0), m)
        e_min = float((z + a).min() / ham.kT)
        out = anneal(ham, state, MoveSettings(n_mag_sweeps=5), np.random.default_rng(0),
                     geometric_schedule(400, 5.0, 1e-3), adapt_every=20, lattice=True)
        self.assertLess(out["best_energy_kBT"] - e_min, 1e-3)
        self.assertLess(abs(out["drift_kBT"]), 1e-9)


if __name__ == "__main__":
    unittest.main()
