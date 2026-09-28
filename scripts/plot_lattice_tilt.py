"""SL-axis tilt at 298 K: distribution, free energy per solid angle, annealed minimum.

    python scripts/plot_lattice_tilt.py outputs/lattice_production 500 outputs/lattice_anneal OUT_DIR

A mean tilt above zero is NOT evidence of a tilted optimum: for a uniformly
random axis p(theta) ~ sin(theta), so <theta> = 57.3 deg with no energy at
all.  The free energy per unit solid angle,

    F(theta) = -kT ln[ p(theta) / sin(theta) ] + const,

removes that geometric factor; its minimum is the most probable axis
direction.  Error bands come from the spread of the per-chain histograms.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main():
    run, equil, anneal_dir, out = Path(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
    out.mkdir(parents=True, exist_ok=True)
    per_chain = []
    for folder in sorted(run.iterdir()):
        if (folder / "complete.json").exists():
            with np.load(folder / "trajectory.npz") as d:
                per_chain.append(d["sl_axis_tilt_deg"][equil:])
    theta = np.concatenate(per_chain)
    edges = np.linspace(0, 90, 46)
    centres = 0.5 * (edges[1:] + edges[:-1])
    sin = np.sin(np.radians(centres))

    def free_energy(x):
        h, _ = np.histogram(x, bins=edges, density=True)
        with np.errstate(divide="ignore"):
            f = -np.log(h / sin)
        return f - np.nanmin(f[np.isfinite(f)])

    F = free_energy(theta)
    Fc = np.array([free_energy(x) for x in per_chain])
    band = np.nanstd(np.where(np.isfinite(Fc), Fc, np.nan), axis=0) / np.sqrt(len(per_chain))

    best = None
    summary_path = anneal_dir / "anneal_summary.json"
    if summary_path.exists():
        best = json.loads(summary_path.read_text(encoding="utf-8"))["summary"]["best"]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].hist(theta, bins=edges, density=True, color="#2a6f97", alpha=0.75, label="298 K sampling")
    iso = sin / np.trapezoid(sin, centres)
    axes[0].plot(centres, iso, "k--", lw=1, label="random axis (sin theta)")
    axes[0].set(xlabel="SL axis tilt theta / deg", ylabel="p(theta)",
                title=f"mean {theta.mean():.2f} deg, median {np.median(theta):.2f} deg")
    ok = np.isfinite(F)
    axes[1].plot(centres[ok], F[ok], "o-", color="#c44e52", ms=3, label="F(theta) / kT")
    axes[1].fill_between(centres[ok], (F - band)[ok], (F + band)[ok], color="#c44e52", alpha=0.2)
    if best is not None:
        axes[1].axvline(best["sl_axis_tilt_deg"], color="k", ls=":", label="annealed minimum of U")
        axes[0].axvline(best["sl_axis_tilt_deg"], color="k", ls=":")
    axes[1].set(xlabel="SL axis tilt theta / deg", ylabel="F / kBT (per solid angle)",
                title="free energy of the SL axis direction")
    axes[1].set_ylim(-0.3, min(12, np.nanmax(F[ok]) + 0.5))
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "sl_tilt_free_energy.png", dpi=150)
    res = dict(n=int(len(theta)), mean_deg=float(theta.mean()), median_deg=float(np.median(theta)),
               q05_deg=float(np.quantile(theta, 0.05)), q95_deg=float(np.quantile(theta, 0.95)),
               F_min_at_deg=float(centres[ok][np.argmin(F[ok])]),
               F_kBT=dict(zip([f"{c:.1f}" for c in centres[ok]], F[ok].round(3).tolist())),
               annealed_tilt_deg=None if best is None else best["sl_axis_tilt_deg"])
    (out / "sl_tilt_free_energy.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps({k: v for k, v in res.items() if k != "F_kBT"}, indent=2))


if __name__ == "__main__":
    main()
