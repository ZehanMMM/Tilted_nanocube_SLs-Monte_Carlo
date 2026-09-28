"""Energy of a collective cube twist about the field-aligned [111] axis.

    python scripts/twist_window.py OUT_DIR

All 27 cubes are rotated together by gamma about the lattice three-fold
axis (parallel to B), positions and moments (along B) fixed.  Gamma is
measured from the V9.03 reference orientation.  The result:

* only a window of about 30 deg is sterically allowed;
* U(gamma) is mirror-symmetric about the window centre (the mirror planes of
  cube and lattice coincide there), so the two lowest-energy states at the
  two walls are mirror images;
* the centre is an energy maximum by only a few kBT for all 27 cubes.

It also overlays the 298 K distribution of single-cube twists (from the
rigid-lattice snapshots) and the annealed (T -> 0) single-cube twists.
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from run_chains import build  # noqa: E402
from v904.model import LATTICE_BODY_REF, initial_state  # noqa: E402


def twists(G, Q):
    rel = np.einsum("ba,nbc,dc->nad", G, Q, LATTICE_BODY_REF)
    return np.degrees(Rotation.from_matrix(rel).as_rotvec()[:, 0])


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "outputs" / "twist_window"
    out.mkdir(parents=True, exist_ok=True)
    ham, _ = build(json.loads((ROOT / "configs" / "lattice_production.json").read_text()))
    st = initial_state(tilt_deg=0.0, dipoles="field")
    grid = np.arange(-5.0, 35.01, 0.5)
    rows = []
    for g in grid:
        R = Rotation.from_rotvec(np.radians(g) * np.array([1.0, 0.0, 0.0])).as_matrix()
        Q = np.repeat((R @ LATTICE_BODY_REF)[None], 27, 0)
        allowed = ham.allowed(st["pos"], Q)
        t = ham.energy_terms(st["pos"], Q, st["mu"])
        steric_ok = t["Steric"] / ham.kT < 1e-6
        rows.append(dict(gamma_deg=float(g), allowed=bool(allowed and steric_ok),
                         U_kBT=t["Total"] / ham.kT, vdW_kBT=t["VdW"] / ham.kT))
    ok = [r for r in rows if r["allowed"]]
    lo, hi = ok[0]["gamma_deg"], ok[-1]["gamma_deg"]
    U = np.array([r["U_kBT"] for r in ok])
    centre = 0.5 * (lo + hi)
    summary = dict(window_deg=[lo, hi], centre_deg=centre,
                   U_wall_kBT=float(U.min()), U_centre_kBT=float(U.max()),
                   barrier_all_27_kBT=float(U.max() - U.min()),
                   barrier_per_cube_kBT=float((U.max() - U.min()) / 27),
                   symmetric=bool(np.allclose(U, U[::-1], atol=1e-6)))
    samp = []
    for f in glob.glob(str(ROOT / "outputs" / "lattice_production_cont" / "*" / "snapshots.npz")):
        with np.load(f) as d:
            for k in range(len(d["cycle"])):
                samp.append(twists(d["G"][k], d["Q"][k]))
    samp = np.concatenate(samp)
    ann = np.concatenate([twists(np.load(f)["G"], np.load(f)["Q"]) for f in
                          glob.glob(str(ROOT / "outputs" / "lattice_anneal" / "*" / "best_state.npz"))])
    summary.update(sampled_mean_deg=float(samp.mean()), sampled_sd_deg=float(samp.std()),
                   sampled_q025_q975_deg=[float(np.quantile(samp, 0.025)), float(np.quantile(samp, 0.975))],
                   annealed_fraction_near_walls=float(np.mean((ann < lo + 8) | (ann > hi - 8))))
    (out / "twist_window.json").write_text(json.dumps(dict(summary=summary, scan=rows), indent=2),
                                           encoding="utf-8")
    print(json.dumps(summary, indent=2))

    fig, ax1 = plt.subplots(figsize=(7.5, 4.3))
    ax1.plot([r["gamma_deg"] for r in ok], U - U.min(), "k-", lw=2,
             label="collective twist: U - U_min (all 27 cubes)")
    ax1.axvspan(grid[0], lo, color="0.85")
    ax1.axvspan(hi, grid[-1], color="0.85", label="sterically forbidden")
    ax1.axvline(centre, color="k", ls=":", lw=1)
    ax1.set(xlabel="cube twist about the field-aligned [111] axis / deg (0 = V9.03 reference)",
            ylabel="U - U_min / kBT", xlim=(grid[0], grid[-1]))
    ax2 = ax1.twinx()
    bins = np.arange(grid[0], grid[-1] + 0.1, 1.5)
    ax2.hist(samp, bins=bins, density=True, color="#2a6f97", alpha=0.45, label="298 K single cubes")
    ax2.hist(ann, bins=bins, density=True, histtype="step", color="#c44e52", lw=1.6,
             label="annealed single cubes (T -> 0)")
    ax2.set_ylabel("probability density")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="upper center")
    ax1.set_title(f"twist window {lo:.1f}-{hi:.1f} deg, mirror-symmetric about {centre:.1f} deg")
    fig.tight_layout()
    fig.savefig(out / "twist_window.png", dpi=160)


if __name__ == "__main__":
    main()
