"""The fast evaluator is the voxel limit of the V10 sharp-cube sum."""
import unittest

import numpy as np
from scipy import integrate
from scipy.spatial.transform import Rotation

from v904.geometry import core_distance
from v904.vdw import (EDGE_M, HAMAKER_J, HamakerEvaluator, HamakerSettings,
                      octree_energy, parallel_voxel_energy, richardson,
                      square_inverse_quartic, uniform_voxel_energy)

KT = 1.380649e-23 * 298.15
I3 = np.eye(3)
FAST = HamakerEvaluator()
TIGHT = HamakerEvaluator(settings=HamakerSettings(q_near=10, tol_kBT=1e-9, q_mid=20,
                                                  q_far=12, far_pair_L=99))


def rot(vec_deg):
    return Rotation.from_rotvec(np.radians(vec_deg)).as_matrix()


def random_pair(rng, min_core_nm, max_core_nm=4.0):
    """A tilted near-contact pair whose sharp cores stay min_core_nm apart."""
    while True:
        Qa, Qb = rot(rng.normal(0, 6, 3)), rot(rng.normal(0, 6, 3))
        r = np.array([EDGE_M + rng.uniform(1.0, 5.0) * 1e-9, *rng.uniform(-4e-9, 4e-9, 2)])
        gap = core_distance(np.zeros(3), Qa, r, Qb, EDGE_M) * 1e9
        if min_core_nm <= gap <= max_core_nm:
            return Qa, r, Qb


class ClosedForm(unittest.TestCase):
    def test_inner_integral_matches_quadrature(self):
        rng = np.random.default_rng(0)
        for _ in range(5):
            u0, v0 = rng.uniform(-1.5, 1.5, 2)
            z = rng.uniform(0.05, 1.0)
            num = integrate.dblquad(
                lambda v, u: 1.0 / ((u - u0) ** 2 + (v - v0) ** 2 + z * z) ** 2,
                -1, 1, -1, 1, epsabs=1e-13, epsrel=1e-11)[0]
            self.assertAlmostEqual(square_inverse_quartic(u0, v0, z * z, 1.0) / num, 1.0,
                                   places=9)


class VoxelLimit(unittest.TestCase):
    def test_v10_numbers(self):
        """n = 4 reproduces the inherited V10 face value; the limit is the audit value."""
        r = np.array([EDGE_M + 3e-9, 0, 0])
        self.assertAlmostEqual(parallel_voxel_energy(r, 4) / -4.883e-21, 1.0, places=3)
        self.assertAlmostEqual(uniform_voxel_energy(np.zeros(3), I3, r, I3, 4,
                                                    d2_floor_m2=1e-19) / -4.883e-21, 1.0, 3)
        # V10 quotes -9.156e-21 for its refined sum; the h -> 0 limit is 0.4 % deeper
        ns = [48, 64, 96]
        vals = [parallel_voxel_energy(r, n) for n in ns]
        self.assertTrue(vals[0] > vals[1] > vals[2])            # monotone deepening
        limit = richardson(vals, ns)
        self.assertLess(abs(FAST.pair_energy(np.zeros(3), I3, r, I3) / limit - 1), 2e-5)

    def test_rotated_pair_voxel_convergence(self):
        """Arbitrary orientation: n^3 midpoint sums converge onto the fast value."""
        Qa, r, Qb = random_pair(np.random.default_rng(11), 1.5, 2.5)
        ns = [8, 12, 16, 24]
        vals = [uniform_voxel_energy(np.zeros(3), Qa, r, Qb, n) for n in ns]
        fast = FAST.pair_energy(np.zeros(3), Qa, r, Qb)
        errors = [abs(v / fast - 1) for v in vals]
        # midpoint error ~ h^2: 16 -> 24 should shrink it by (24/16)^2 = 2.25
        self.assertTrue(all(e1 > e2 for e1, e2 in zip(errors, errors[1:])))
        self.assertAlmostEqual(errors[2] / errors[3], 2.25, delta=0.3)
        rich = [abs(richardson(vals[:k + 1], ns[:k + 1]) / fast - 1) for k in (1, 2, 3)]
        self.assertTrue(rich[0] > rich[1] > rich[2])
        self.assertLess(rich[-1], 0.01)

    def test_fast_matches_tight_reference(self):
        rng = np.random.default_rng(3)
        for _ in range(12):
            Qa, r, Qb = random_pair(rng, 0.3)
            e_fast = FAST.pair_energy(np.zeros(3), Qa, r, Qb)
            e_ref = TIGHT.pair_energy(np.zeros(3), Qa, r, Qb)
            self.assertLess(abs(e_fast - e_ref) / KT, 2e-5)


def vertex_contact(delta_nm, seed=1):
    """Cube B's most exposed vertex placed delta from cube A's (+,+,+) vertex."""
    a = EDGE_M / 2
    verts = np.array([[x, y, z] for x in (-a, a) for y in (-a, a) for z in (-a, a)])
    QB = Rotation.random(random_state=np.random.default_rng(seed)).as_matrix()
    u = np.ones(3) / np.sqrt(3)
    vb = verts @ QB.T
    return QB, np.array([a, a, a]) + delta_nm * 1e-9 * u - vb[np.argmin(vb @ u)]


class NearContact(unittest.TestCase):
    """Sharp cores a few tenths of a nm apart: the regime the MC reaches."""

    def test_matches_octree_volume_integral(self):
        """Independent method: value-blind octree refinement of the volume form."""
        for delta in (0.165, 0.3):
            QB, pB = vertex_contact(delta)
            fast = HamakerEvaluator()
            e = fast.pair_energy(np.zeros(3), I3, pB, QB)
            ref = octree_energy(np.zeros(3), I3, pB, QB, eta=0.5, order=2)
            self.assertLess(abs(e / ref - 1), 2e-3)
            self.assertEqual(fast.max_level_hits, 0)

    def test_value_only_refinement_is_not_trusted(self):
        """A narrow spike could fool a value-only test; the default hybrid rule
        forces refinement near the other face whatever the values say."""
        s = HamakerSettings()
        self.assertEqual(s.scheme, "hybrid")
        self.assertLessEqual(s.eta_value, 4.0)


class Invariances(unittest.TestCase):
    def setUp(self):
        self.Qa, self.r, self.Qb = random_pair(np.random.default_rng(5), 1.0)
        self.e = FAST.pair_energy(np.zeros(3), self.Qa, self.r, self.Qb)

    def test_rigid_motion(self):
        G = rot([30, -50, 20])
        t = np.array([4e-9, -7e-9, 2e-9])
        e = FAST.pair_energy(t, G @ self.Qa, t + G @ self.r, G @ self.Qb)
        self.assertLess(abs(e / self.e - 1), 1e-9)

    def test_exchange(self):
        e = FAST.pair_energy(self.r, self.Qb, np.zeros(3), self.Qa)
        self.assertLess(abs(e / self.e - 1), 1e-6)

    def test_cube_symmetry(self):
        C = rot([90, 0, 0])                                   # maps the cube to itself
        e = FAST.pair_energy(np.zeros(3), self.Qa @ C, self.r, self.Qb)
        self.assertLess(abs(e / self.e - 1), 1e-9)

    def test_far_field(self):
        d = 20 * EDGE_M
        e = FAST.pair_energy(np.zeros(3), self.Qa, np.array([d, 0, 0]), self.Qb)
        point = -HAMAKER_J / np.pi ** 2 * EDGE_M ** 6 / d ** 6
        self.assertLess(abs(e / point - 1), 0.01)


if __name__ == "__main__":
    unittest.main()
