"""The Hamiltonian is V9.03's except for vdW, and every move's dE is exact."""
import contextlib
import io
import unittest
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from v904 import geometry
from v904.model import (Hamiltonian, INITIAL_STRUCTURES, PhysicalParameters,
                        initial_state)
from v904.sampler import Chain, MoveSettings

HERE = Path(__file__).resolve().parents[1]
# bundled copy of the two V9.03 cells; the sibling folder is the original location
V903 = next((p for p in (HERE / "reference" / "v903_cells",
                         HERE.parent / "V9.03_N27_collective_mc_500G" / "_v903cells")
             if (p / "cell2.py").exists()), HERE / "reference" / "v903_cells")


def load_v903():
    ns = {}
    with contextlib.redirect_stdout(io.StringIO()):
        for name in ("cell1.py", "cell2.py"):
            path = V903 / name
            exec(compile(path.read_text(encoding="utf-8"), str(path), "exec"), ns)
    return ns["CollectiveNanocubeMC"](16.0, ns["PHYS_500G"], anis_model="cubic_first_raw")


def jiggle(state, rng, trans_nm=0.3, rot_deg=3.0):
    pos = state["pos"] + rng.normal(0, trans_nm * 1e-9, state["pos"].shape)
    pos[0] = state["pos"][0]
    Q = np.array([Rotation.from_rotvec(np.radians(rot_deg) * rng.normal(size=3)).as_matrix() @ q
                  for q in state["Q"]])
    mu = rng.normal(size=state["mu"].shape)
    return dict(pos=pos, Q=Q, mu=mu / np.linalg.norm(mu, axis=1, keepdims=True))


class MatchesV903(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.old = load_v903()
        cls.ham = Hamiltonian(PhysicalParameters(vdw_model="inherited"))

    def test_every_term(self):
        rng = np.random.default_rng(2)
        for name, cfg in INITIAL_STRUCTURES.items():
            s = jiggle(initial_state(**cfg, dipole_seed=3), rng, trans_nm=0.05)
            old = self.old.energy_full_cluster(s["pos"], s["Q"], s["mu"], a_nm=21.0)
            new = self.ham.energy_terms(s["pos"], s["Q"], s["mu"])
            for key in ("Zeeman", "Anisotropy", "Dipole", "VdW", "Steric", "Total"):
                a, b = new[key] / self.ham.kT, old[key] / self.old.kT
                # jiggled states can overlap (steric ~ 1e11 kBT), so compare relatively
                self.assertLessEqual(abs(a - b), 1e-9 * max(1.0, abs(b)), msg=f"{name} {key}")

    def test_initial_structures_identical(self):
        # V9.03 initial_state with random dipoles drawn from default_rng(seed)
        s = initial_state(**INITIAL_STRUCTURES["expanded_tilted"], dipole_seed=7)
        old = self.old.initial_state(21.0, 74.2, beta_deg=40.0, n=1, convention="co_rotate",
                                     dipoles="random", dipole_rng=np.random.default_rng(7))
        old_pos = old["pos"][0] + 1.08 * (old["pos"] - old["pos"][0])
        np.testing.assert_allclose(s["pos"], old_pos, atol=1e-20)
        np.testing.assert_allclose(s["Q"], old["Q"], atol=1e-14)
        np.testing.assert_allclose(s["mu"], old["mu"], atol=1e-14)

    def test_steric_gap(self):
        rng = np.random.default_rng(4)
        s = jiggle(initial_state(dipole_seed=1), rng)
        for i, j in [(0, 5), (3, 17), (8, 26)]:
            r = s["pos"][j] - s["pos"][i]
            g_old = self.old.surface_gap_pair(r, s["Q"][i], s["Q"][j])
            g_new = geometry.support_gap(r, s["Q"][i], s["Q"][j], self.ham.L,
                                         self.ham.r_round)[0]
            self.assertAlmostEqual(g_new / g_old, 1.0, places=12)


class Constraints(unittest.TestCase):
    def test_overlap_agrees_with_exact_distance(self):
        rng = np.random.default_rng(9)
        L = 16e-9
        for _ in range(300):
            Qa = Rotation.random(random_state=rng).as_matrix()
            Qb = Rotation.random(random_state=rng).as_matrix()
            r = rng.normal(size=3)
            r *= rng.uniform(0.8, 1.9) * L / np.linalg.norm(r)
            d = geometry.core_distance(np.zeros(3), Qa, r, Qb, L)
            if d > 1e-12:
                self.assertFalse(geometry.cores_overlap(r, Qa, Qb, L)[0])
            elif geometry.cores_overlap(r, Qa, Qb, L)[0] is np.False_:
                self.fail("overlap missed by the separating-axis test")

    def test_minimum_core_separation(self):
        """Cores closer than D0 = 0.165 nm are excluded, farther ones allowed."""
        from tests.test_vdw import vertex_contact
        ham = Hamiltonian()
        for delta, excluded in ((0.10, True), (0.16, True), (0.17, False), (0.40, False)):
            QB, pB = vertex_contact(delta)
            self.assertEqual(ham.overlap_any(np.zeros(3), np.eye(3), pB, QB), excluded,
                             msg=f"delta = {delta} nm")

    def test_connectivity(self):
        pos = np.array([[0, 0, 0], [25e-9, 0, 0], [50e-9, 0, 0]])
        self.assertTrue(geometry.connected(pos, 30e-9))
        self.assertFalse(geometry.connected(pos, 20e-9))


class ExactDeltas(unittest.TestCase):
    """Cached energies after accepted moves equal a fresh evaluation."""

    def test_all_moves_keep_bookkeeping_exact(self):
        ham = Hamiltonian()
        s = initial_state(**INITIAL_STRUCTURES["compact_tilted40"], dipole_seed=2)
        # large steps so every move type is accepted a few times
        moves = MoveSettings(trans_step_nm=0.3, rot_step_deg=4, cotilt_step_deg=3,
                             gamma_step_deg=3, scale_log_step=0.003, n_mag_sweeps=1)
        chain = Chain(ham, s, moves, np.random.default_rng(1))
        for _ in range(3):
            chain.cycle()
        self.assertTrue(all(chain.accepts[k] > 0 for k in ("trans", "rot", "dip", "cotilt")))
        self.assertLess(abs(chain.drift_kBT()), 1e-6)

    def test_cotilt_leaves_vdw_and_steric_unchanged(self):
        ham = Hamiltonian()
        s = jiggle(initial_state(dipole_seed=2), np.random.default_rng(3), trans_nm=0.05)
        W0, S0, _ = ham.pair_matrices(s["pos"], s["Q"])
        G = Rotation.from_rotvec([0.3, -0.2, 0.5]).as_matrix()
        pos = s["pos"][0] + (s["pos"] - s["pos"][0]) @ G.T
        Q = np.einsum("ab,nbc->nac", G, s["Q"])
        W1, S1, _ = ham.pair_matrices(pos, Q)
        self.assertLess(np.max(np.abs(W1 - W0)) / ham.kT, 1e-5)
        self.assertLess(np.max(np.abs(S1 - S0) / np.maximum(np.abs(S0), ham.kT)), 1e-9)

    def test_gamma_axis_is_preserved(self):
        """Reverse twist uses the same axis, so the proposal is symmetric."""
        s = jiggle(initial_state(dipole_seed=2), np.random.default_rng(5), trans_nm=0.0)
        e111 = np.ones(3) / np.sqrt(3)
        axis = np.einsum("nij,j->ni", s["Q"], e111).mean(0)
        axis /= np.linalg.norm(axis)
        G = Rotation.from_rotvec(0.4 * axis).as_matrix()
        new_axis = np.einsum("nij,j->ni", np.einsum("ab,nbc->nac", G, s["Q"]), e111).mean(0)
        np.testing.assert_allclose(new_axis / np.linalg.norm(new_axis), axis, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
