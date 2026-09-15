"""Short regression tests only. No production or full multi-seed runs."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NS = {}
with contextlib.redirect_stdout(io.StringIO()):
    for name in ["cell1.py", "cell2.py", "cell2a.py", "cell3.py"]:
        exec(compile((ROOT / "_v903cells" / name).read_text(encoding="utf-8"), name, "exec"), NS)
NS.update(SHOW_PROGRESS_BARS=False, SHOW_FIGURES=False, SAVE_REPORT_PDFS=False)


class SamplingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = NS["model"]
        cls.state = NS["state0"]

    def run_short(self, **kwargs):
        settings = dict(n_cycles=2, n_equil=0, n_mag_per_cycle=5,
                        n_cotilt_per_cycle=1, n_gamma_per_cycle=1,
                        rng=np.random.default_rng(17))
        settings.update(kwargs)
        return NS["run_local_collective_mc"](self.model, self.state, NS["A_NM"], **settings)

    def test_energy_model_unchanged_at_initial_state(self):
        self.assertAlmostEqual(NS["E0"]["Total"] / self.model.kT, -120.36913284974837, places=6)
        self.assertEqual(NS["PHYS_500G"]["k_stiff"], 1e8)
        self.assertEqual(NS["PHYS_500G"]["Hamaker"], 2e-20)

    def test_magnetic_delta_matches_full_for_every_particle(self):
        m, s = self.model, self.state
        old = m.energy_full_cluster(s["pos"], s["Q"], s["mu"], a_nm=NS["A_NM"])["Total"]
        rotation = NS["R"].from_euler("y", 11, degrees=True).as_matrix()
        for i in range(len(s["pos"])):
            mu = s["mu"].copy()
            e0 = m.local_magnetic_energy(i, s["pos"], s["Q"], mu)
            mu[i] = rotation @ mu[i]
            delta = m.local_magnetic_energy(i, s["pos"], s["Q"], mu) - e0
            full = m.energy_full_cluster(s["pos"], s["Q"], mu, a_nm=NS["A_NM"])["Total"] - old
            self.assertAlmostEqual((delta - full) / m.kT, 0, places=10)

    def test_complete_sweeps_and_energy_tracking(self):
        r = self.run_short()
        self.assertEqual(r["attempts"]["dipole"], 2 * 5 * 27)
        self.assertEqual(r["attempts"]["mechanical"], 2 * 27)
        exact = self.model.energy_full_cluster(r["pos"], r["Q"], r["mu"], a_nm=NS["A_NM"])
        self.assertAlmostEqual(exact["Total"] / self.model.kT, r["traj_E"][-1], places=8)
        self.assertGreater(np.min(r["traj_gap_min_nm"]), 0)
        np.testing.assert_allclose(np.linalg.norm(r["mu"], axis=1), 1, atol=1e-12)
        np.testing.assert_allclose(r["Q"] @ r["Q"].transpose(0, 2, 1),
                                   np.tile(np.eye(3), (27, 1, 1)), atol=1e-12)

    def test_switches_cover_global_moves(self):
        r = self.run_short(move_positions=False, move_orientations=False)
        np.testing.assert_array_equal(r["pos"], self.state["pos"])
        np.testing.assert_array_equal(r["Q"], self.state["Q"])
        self.assertEqual(r["attempts"]["mechanical"], 0)
        self.assertEqual(r["attempts"]["cotilt"], 0)
        self.assertEqual(r["attempts"]["gamma"], 0)
        r = self.run_short(move_positions=True, move_orientations=False)
        np.testing.assert_array_equal(r["Q"], self.state["Q"])
        self.assertEqual(r["attempts"]["cotilt"], 0)

    def test_vdw_switch_energy_consistency(self):
        r = self.run_short(include_vdw=False)
        e = self.model.energy_full_cluster(r["pos"], r["Q"], r["mu"],
                                          a_nm=NS["A_NM"], include_vdw=False)
        self.assertEqual(e["VdW"], 0)
        self.assertAlmostEqual(e["Total"] / self.model.kT, r["traj_E"][-1], places=8)

    def test_observables_have_distinct_definitions(self):
        q = self.state["Q"][:1]
        easy = q[0] @ self.model.e111
        perpendicular = q[0] @ self.model.body_u
        for mu, expected in [(easy, 0), (perpendicular, 90), (-easy, 180)]:
            actual = self.model.local_beta_angles(q, mu[None])[0]
            self.assertAlmostEqual(actual, expected, places=5)
        r = self.run_short()
        np.testing.assert_allclose(r["traj_body_sl_mismatch"],
                                   r["traj_body_tilt"] - r["traj_sl_pca_tilt"])
        self.assertTrue(np.all(r["traj_body_sl_axis_angle"] >= 0))

    def test_reproducibility_and_input_immutability(self):
        saved = {k: v.copy() for k, v in self.state.items()}
        a, b = self.run_short(), self.run_short()
        for key in NS["DIAGNOSTIC_TRAJECTORIES"].values():
            np.testing.assert_array_equal(a[key], b[key])
        for key, value in saved.items():
            np.testing.assert_array_equal(self.state[key], value)

    def test_mh_asymmetric_proposal_stationary_distribution(self):
        # pi(1)=1/3, q(1|0)=0.8 and q(0|1)=0.2. Omitting Hastings fails this test.
        rng = np.random.default_rng(19)
        state, samples = 0, []
        for i in range(16000):
            q_forward = .8 if state == 0 else .2
            if rng.random() < q_forward:
                q_reverse = .2 if state == 0 else .8
                delta = np.log(2) if state == 0 else -np.log(2)
                if NS["metropolis_hastings_accept"](
                        delta, 1, rng, np.log(q_reverse / q_forward)):
                    state = 1 - state
            if i >= 1000:
                samples.append(state)
        self.assertLess(abs(np.mean(samples) - 1 / 3), .025)

    def test_mh_extreme_ratios_and_invalid_inputs(self):
        accept = NS["metropolis_hastings_accept"]
        rng = np.random.default_rng(5)
        self.assertFalse(accept(np.inf, 1, rng))
        self.assertTrue(accept(-np.inf, 1, rng))
        with self.assertRaises(ValueError):
            accept(np.nan, 1, rng)
        with self.assertRaises(ValueError):
            self.run_short(n_cycles=0)
        with self.assertRaises(ValueError):
            self.run_short(n_equil=2)
        with self.assertRaises(ValueError):
            self.run_short(n_mag_per_cycle=1.5)

    def test_initial_structures_keep_body_and_sl_aligned(self):
        states = [NS["make_initial_state"](1, name) for name in NS["MULTI_STARTS"]]
        for s in states:
            self.assertEqual(len(s["pos"]), 27)
            self.assertGreater(self.model.surface_gap_stats(s["pos"], s["Q"])["min_gap_nm"], 0)
            self.assertLess(self.model.body_sl_mismatch_deg(s["pos"], s["Q"])["coherent_pca_deg"], 1e-5)
        d0 = np.linalg.norm(states[0]["pos"], axis=1)
        d2 = np.linalg.norm(states[2]["pos"], axis=1)
        np.testing.assert_allclose(d2, 1.08 * d0)

    def test_ensemble_and_export(self):
        b = NS["run_chain_ensemble"](NS["MULTI_STARTS"], [1, 2], 2, 0, 1)
        self.assertEqual(b["draws"]["energy_kBT"].shape, (6, 2))
        self.assertEqual(len({tuple(r["sample_entropy"]) for r in b["rows"]}), 6)
        self.assertEqual(len(b["grouped_diagnostics"]), 3)
        self.assertTrue(all(r["status"] == "insufficient_draws" for r in b["diagnostics"]))
        with tempfile.TemporaryDirectory() as tmp:
            prefix = str(Path(tmp) / "smoke")
            NS["export_chain_ensemble"](b, prefix, report_pdf=prefix + ".pdf")
            metadata = json.loads(Path(prefix + "_metadata.json").read_text())
            self.assertEqual(metadata["n_cycles"], 2)
            with np.load(prefix + "_trajectories.npz") as data:
                self.assertEqual(data["traj_E"].shape, (6, 2))
                np.testing.assert_array_equal(data["traj_E"][:, 0],
                                               [r["traj_E"][0] for r in b["results"]])
            self.assertIn("insufficient_draws", Path(prefix + "_diagnostics.csv").read_text())
            self.assertGreater(Path(prefix + ".pdf").stat().st_size, 1000)
        with self.assertRaises(ValueError):
            NS["run_chain_ensemble"](NS["MULTI_STARTS"], [1, 1], 2, 0, 1)

    def test_short_pilot_never_claims_an_optimum(self):
        pilot = NS["compare_magnetic_sweeps"](["compact_aligned"], [1], [1, 5], 2, 0)
        self.assertIsNone(pilot["recommendation"])
        self.assertTrue(all(not row["eligible"] for row in pilot["rows"]))
        self.assertEqual(NS["N_MAG_PER_CYCLE"], 5)


class DiagnosticTests(unittest.TestCase):
    def test_iid_normal_reference(self):
        x = np.random.default_rng(40).normal(size=(4, 2000))
        r = NS["diagnose_chains"]({"x": x}, 2)[0]
        self.assertLess(r["rhat_rank"], 1.01)
        self.assertGreater(r["ess_bulk"], 5000)
        self.assertAlmostEqual(r["mcse_mean"], 1 / np.sqrt(x.size), delta=.002)
        self.assertEqual(r["status"], "checks_passed")

    def test_shifted_and_correlated_chains(self):
        rng = np.random.default_rng(9)
        shifted = rng.normal(size=(4, 2000))
        shifted[0] += 3
        r = NS["diagnose_chains"]({"x": shifted})[0]
        self.assertGreater(r["rhat_rank"], 1.1)
        self.assertIn("high_rhat", r["status"])
        ar = rng.normal(size=(4, 4000))
        for i in range(1, ar.shape[1]):
            ar[:, i] += .97 * ar[:, i - 1]
        r = NS["diagnose_chains"]({"x": ar[:, 1000:]})[0]
        self.assertLess(r["ess_mean"], 1000)
        self.assertGreater(r["mcse_mean"], np.std(ar[:, 1000:]) / np.sqrt(12000) * 3)

    def test_pathological_inputs(self):
        cases = [
            (np.ones((4, 200)), "constant_chain"),
            (np.arange(800).reshape(4, 200) * np.nan, "nonfinite_data"),
            (np.ones((4, 3)), "insufficient_draws"),
            (np.ones((1, 200)), "insufficient_chains"),
        ]
        for x, expected in cases:
            r = NS["diagnose_chains"]({"x": x})[0]
            self.assertEqual(r["status"], expected)
            self.assertTrue(np.isnan(r["mcse_mean"]))
            self.assertEqual(NS["report_value"](r["mcse_mean"]), "N/A")
        with self.assertRaises(ValueError):
            NS["diagnose_chains"]({"x": np.zeros(100)})

    def test_notebook_sources_and_schema(self):
        import nbformat
        raw = json.loads((ROOT / "V9.03_N27_500G.ipynb").read_text(encoding="utf-8"))
        ids = [cell["id"] for cell in raw["cells"]]
        self.assertEqual(len(ids), len(set(ids)))
        nb = nbformat.read(ROOT / "V9.03_N27_500G.ipynb", as_version=4)
        nbformat.validate(nb)
        sources = sorted((ROOT / "_v903cells").glob("cell*.*"))
        self.assertEqual(len(nb.cells), len(sources))
        for cell, source in zip(nb.cells, sources):
            self.assertEqual(cell.source, source.read_text(encoding="utf-8"))
            if cell.cell_type == "code":
                compile(cell.source, source.name, "exec")
                self.assertEqual(cell.outputs, [])


if __name__ == "__main__":
    unittest.main()
