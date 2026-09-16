"""Plot saved joint chains without rerunning MC or deleting rejected states."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    chains = []
    for folder in sorted(args.input.iterdir()):
        if not folder.is_dir() or not (folder / "complete.json").exists():
            continue
        info = json.loads((folder / "result.json").read_text(encoding="utf-8"))
        if "job" in info:
            identity = info["job"]
            equil = info["sampling"]["n_equil"]
            if not info["sampling"]["move_orientations"]:
                raise ValueError("This plot is for joint sampling, not frozen bodies.")
        else:
            if info["run"]["kind"] != "multi":
                continue
            identity = info["run"]["row"]
            equil = info["run"]["metadata"]["n_equil"]
        with np.load(folder / "result.npz") as raw:
            data = {k: raw[k].copy() for k in
                    ["traj_body_tilt", "traj_sl_pca_tilt", "traj_rg_nm", "traj_E"]}
        chains.append(dict(start=identity["start"], seed=identity["seed"],
                           sweeps=identity.get("sweeps", 5), equil=equil, data=data))
    if not chains:
        raise ValueError("No completed joint chains found.")
    if len({c["sweeps"] for c in chains}) != 1:
        raise ValueError("Plot each sweep candidate separately.")
    groups = list(dict.fromkeys(c["start"] for c in chains))
    # Separate body and SL figures. Share axes only within each variable.
    for key, label in [("traj_body_tilt", "Coherent body tilt (deg)"),
                       ("traj_sl_pca_tilt", "SL PCA tilt (deg)"),
                       ("traj_rg_nm", "Radius of gyration (nm)"),
                       ("traj_E", "Energy / kBT")]:
        fig, axes = plt.subplots(len(groups), 1, figsize=(9, 2.3 * len(groups)),
                                 sharex=True, sharey=True, constrained_layout=True)
        for ax, group in zip(np.atleast_1d(axes), groups):
            for chain in chains:
                if chain["start"] == group:
                    y = chain["data"][key]
                    ax.plot(np.arange(1, len(y) + 1), y, lw=0.6, label=f"seed {chain['seed']}")
            ax.axvline(chain["equil"], color="black", ls="--", lw=0.8)
            ax.set_title(group.replace("_", " "), fontsize=10)
            ax.set_ylabel(label, fontsize=9)
            ax.grid(alpha=0.25)
            ax.legend(fontsize=8, ncol=4, frameon=False)
        np.atleast_1d(axes)[-1].set_xlabel("MC cycle (not physical time)")
        fig.savefig(args.output / f"{key}.png", dpi=160)
        plt.close(fig)
    fig, axes = plt.subplots(len(groups), 1, figsize=(8, 2.3 * len(groups)),
                             sharex=True, sharey=True, constrained_layout=True)
    for ax, group in zip(np.atleast_1d(axes), groups):
        values = [c["data"]["traj_body_tilt"][c["equil"]:] for c in chains if c["start"] == group]
        for half, label in [(0, "First retained half"), (1, "Second retained half")]:
            data = np.concatenate([np.array_split(v, 2)[half] for v in values])
            x = np.sort(data)
            ax.plot(x, np.arange(1, len(x) + 1) / len(x), label=label)
        ax.set_title(group.replace("_", " "), fontsize=10)
        ax.set_ylabel("Empirical CDF")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8, frameon=False)
    np.atleast_1d(axes)[-1].set_xlabel("Coherent body tilt (deg)")
    fig.savefig(args.output / "body_tilt_window_cdf.png", dpi=160)
    plt.close(fig)
    print(f"Saved five stability figures for {len(chains)} chains to {args.output}")


if __name__ == "__main__":
    main()
