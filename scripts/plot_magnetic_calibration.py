"""Compare saved magnetic pilots, without pooling different frozen geometries."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=args.overwrite)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for root in args.inputs:
        config = json.loads((root / "manifest.json").read_text())["config"]
        step = config.get("dip_step_rad", 0.30)
        rows = list(csv.DictReader((root / "candidates.csv").open()))
        sweeps = [int(r["sweeps"]) for r in rows]
        axes[0].plot(sweeps, [float(r["max_rhat"]) for r in rows], "o-", label=f"step {step:.2f} rad")
        axes[1].plot(sweeps, [float(r["min_bulk_ess_per_sec"]) for r in rows], "o-", label=f"step {step:.2f} rad")
        for x, row in zip(sweeps, rows):
            axes[1].annotate(f"{row['failed_checks']} flags", (x, float(row["min_bulk_ess_per_sec"])),
                             xytext=(0, 8), textcoords="offset points", fontsize=8, ha="center")
        trace, grid = plt.subplots(len(sweeps), 2, figsize=(11, 3 * len(sweeps)), constrained_layout=True)
        grid = np.asarray(grid).reshape(len(sweeps), 2)
        for i, count in enumerate(sweeps):
            retained = []
            for seed in config["seeds"]:
                folder = root / f"m{count}_compact_tilted40_s{seed}"
                with np.load(folder / "result.npz") as data:
                    y = data["traj_particle_beta"][config["equil"]:, 12]
                retained.append(y)
                grid[i, 0].plot(np.arange(config["equil"] + 1, config["cycles"] + 1), y,
                                lw=0.4, alpha=0.75, label=f"seed {seed}")
                grid[i, 1].hist(y, bins=np.arange(0, 181, 3), histtype="step", density=True,
                                label=f"seed {seed}")
            grid[i, 0].set(title=f"{count} sweeps: particle 12, compact 40 deg",
                           xlabel="MC cycle", ylabel="Directed local beta (deg)")
            grid[i, 1].set(xlabel="Directed local beta (deg)", ylabel="Sample density")
            grid[i, 0].legend(fontsize=7, frameon=False, ncol=4)
            grid[i, 1].legend(fontsize=7, frameon=False)
        trace.suptitle(f"Magnetic proposal step {step:.2f} rad: retained samples, not an equilibrium claim")
        trace.savefig(args.output / f"particle_beta_step_{step:.2f}.png", dpi=150)
        plt.close(trace)
    axes[0].axhline(1.01, color="black", ls="--", lw=0.8, label="R-hat threshold")
    axes[0].set(xlabel="Magnetic sweeps per cycle", ylabel="Maximum rank R-hat")
    axes[1].set(xlabel="Magnetic sweeps per cycle", ylabel="Worst bulk ESS / sampling second")
    axes[1].margins(x=0.12, y=0.20)
    for ax in axes:
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    fig.savefig(args.output / "calibration_comparison.png", dpi=170)
    plt.close(fig)
    print(f"Saved calibration figures to {args.output}")


if __name__ == "__main__":
    main()
