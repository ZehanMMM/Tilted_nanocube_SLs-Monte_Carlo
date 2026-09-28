"""Full history of every chain: v1 relaxation (3000) followed by v2 (4000).

    python scripts/plot_history.py outputs/production_v1_energy_hole outputs/production OUT.png

v1 sampled a slightly different target (no D0 core minimum), so its part is
shown only as the relaxation history; statistics use v2 alone.  The chain
that fell into the energy hole restarts from its cycle-1300 snapshot, so its
v1 segment is cut there.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

KEYS = [("rg_nm", "Rg / nm"), ("vdw_kBT", "vdW / kBT"), ("n_bonds", "bonds (< 30 nm)"),
        ("energy_kBT", "U / kBT"), ("sl_pca_tilt_deg", "SL PCA tilt / deg"),
        ("body_tilt_deg", "coherent body tilt / deg"), ("magnetization", "m_x")]
COLORS = {"compact_aligned": "#2a6f97", "compact_tilted40": "#c44e52",
          "expanded_aligned": "#55a868", "expanded_tilted": "#dd8452"}


def main():
    v1, v2, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    sources = json.loads((Path("outputs/v2_sources") / "SOURCES.json").read_text())
    fig, axes = plt.subplots(len(KEYS), 1, figsize=(12, 2.1 * len(KEYS)), sharex=True)
    for folder in sorted(v2.iterdir()):
        if not (folder / "complete.json").exists():
            continue
        info = json.loads((folder / "chain.json").read_text())
        start, seed = info["start"], info["seed"]
        old = v1 / f"{start}_s{seed - 10}"
        cut = int(sources[old.name]["source"].rsplit("_", 1)[1])
        with np.load(old / "trajectory.npz") as a, np.load(folder / "trajectory.npz") as b:
            for ax, (k, label) in zip(axes, KEYS):
                x = np.r_[a[k][:cut], b[k]]
                t = np.r_[np.arange(1, cut + 1), 3000 + np.arange(1, len(b[k]) + 1)]
                if cut < 3000:        # hole chain: show the discarded part faintly
                    ax.plot(np.arange(cut + 1, 3001), a[k][cut:], color="0.75", lw=0.4)
                ax.plot(t, x, lw=0.45, color=COLORS[start], alpha=0.8)
                ax.set_ylabel(label, fontsize=8)
    for ax in axes:
        ax.axvline(3000, color="k", lw=0.8)
        ax.axvline(4000, color="k", lw=0.8, ls="--")
    axes[0].set_title("v1 relaxation (cycles 1-3000) | v2 production; dashed: end of v2 warm-up",
                      fontsize=10)
    axes[-1].set_xlabel("MC cycle")
    handles = [plt.Line2D([], [], color=c, label=s) for s, c in COLORS.items()]
    axes[0].legend(handles=handles, fontsize=7, ncol=4, loc="lower right")
    fig.tight_layout()
    fig.savefig(out, dpi=130)


if __name__ == "__main__":
    main()
