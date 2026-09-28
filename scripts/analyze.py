"""Diagnose a finished run and write REPORT.md, CSV tables and figures.

    python scripts/analyze.py RUN_DIR --equil 1000

The warm-up length is fixed BEFORE looking at the output (--equil); the
script then checks that choice three ways instead of tuning it to the
result: MSER-5 truncation per chain, Geweke z on the kept window, and
first-half vs second-half means per start.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from v904.diagnostics import autocorrelation, diagnose, mser_burn_in  # noqa: E402
from v904.sampler import LATTICE_OBSERVABLES, OBSERVABLES  # noqa: E402

KEY = ["energy_kBT", "vdw_kBT", "body_tilt_deg", "sl_pca_tilt_deg", "rg_nm",
       "magnetization", "local_beta_deg", "body_order", "n_bonds", "min_gap_nm"]
# rigid-lattice runs: Rg and bonds are constant, the lattice observables matter
KEY_LATTICE = ["energy_kBT", "sl_axis_tilt_deg", "body_rotation_deg", "body_tilt_deg",
               "body_lattice_angle_deg", "magnetization", "local_beta_deg", "dipole_kBT",
               "anisotropy_kBT", "vdw_kBT"]


def key_list(chains):
    return KEY_LATTICE if "sl_axis_tilt_deg" in chains[0]["traj"] else KEY
MCSE_TARGETS = {"body_tilt_deg": 0.5, "sl_pca_tilt_deg": 0.5, "local_beta_deg": 0.25,
                "magnetization": 0.002}
COLORS = {"compact_aligned": "#2a6f97", "compact_tilted40": "#c44e52",
          "expanded_aligned": "#55a868", "expanded_tilted": "#dd8452",
          "lat_tilt00": "#2a6f97", "lat_tilt20": "#55a868", "lat_tilt40": "#c44e52",
          "lat_tilt60": "#dd8452"}


def load(run):
    chains = []
    for folder in sorted(Path(run).iterdir()):
        if not (folder / "complete.json").exists():
            continue
        info = json.loads((folder / "chain.json").read_text(encoding="utf-8"))
        with np.load(folder / "trajectory.npz") as d:
            traj = {k: d[k].copy() for k in d.files}
        chains.append(dict(name=folder.name, info=info, traj=traj))
    return chains


def fmt(x, nd=4):
    if isinstance(x, str):
        return x
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    return f"{x:.{nd}g}"


def table(rows, cols):
    head = "| " + " | ".join(cols) + " |\n|" + "---|" * len(cols) + "\n"
    return head + "".join("| " + " | ".join(fmt(r.get(c, "")) for c in cols) + " |\n" for r in rows)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run")
    p.add_argument("--equil", type=int, required=True)
    p.add_argument("--title", default=None)
    a = p.parse_args()
    run = Path(a.run)
    chains = load(run)
    if not chains:
        raise SystemExit("no complete chains")
    n = min(len(c["traj"]["energy_kBT"]) for c in chains)
    starts = sorted({c["info"]["start"] for c in chains})
    E = a.equil
    names = [k for k in OBSERVABLES + LATTICE_OBSERVABLES if k in chains[0]["traj"]]
    draws = {k: np.stack([c["traj"][k][E:n] for c in chains]) for k in names}
    KEYS = key_list(chains)
    # constant observables (Rg in a rigid lattice) carry no diagnostic information
    draws = {k: v for k, v in draws.items() if np.ptp(v) > 1e-9 * max(1.0, np.abs(v).max())}

    pooled = diagnose(draws, mcse_targets=MCSE_TARGETS)
    for r in pooled:
        r["scope"] = "pooled"
    per_start = []
    for s in starts:
        idx = [i for i, c in enumerate(chains) if c["info"]["start"] == s]
        for r in diagnose({k: v[idx] for k, v in draws.items()}, mcse_targets=MCSE_TARGETS):
            r["scope"] = s
            per_start.append(r)
    all_rows = pooled + per_start
    with open(run / "diagnostics.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(all_rows[0]))
        w.writeheader()
        w.writerows(all_rows)

    # burn-in checks and window stability
    burn, windows = [], []
    for c in chains:
        row = dict(chain=c["name"])
        for k in [k for k in ("energy_kBT", "rg_nm", "body_tilt_deg", "vdw_kBT",
                              "sl_axis_tilt_deg") if k in draws]:
            row[f"mser_{k}"] = mser_burn_in(c["traj"][k][:n])
        burn.append(row)
    half = E + (n - E) // 2
    for s in starts:
        idx = [i for i, c in enumerate(chains) if c["info"]["start"] == s]
        for k in KEYS:
            x1 = np.mean([chains[i]["traj"][k][E:half] for i in idx])
            x2 = np.mean([chains[i]["traj"][k][half:n] for i in idx])
            windows.append(dict(start=s, variable=k, first_half=x1, second_half=x2,
                                change=x2 - x1))

    acc_rows = []
    for c in chains:
        info = c["info"]
        acc_rows.append(dict(chain=c["name"], minutes=info["elapsed_sec"] / 60,
                             **{f"acc_{m}": v for m, v in info["acceptance"].items()},
                             **{f"blocked_{m}": info["blocked"][m] / max(info["tries"][m], 1)
                                for m in info["blocked"]},
                             drift_kBT=info["max_abs_energy_drift_kBT"]))

    # ------------------------------------------------------------------ figures
    fig, axes = plt.subplots(len(KEYS), 1, figsize=(11, 2.0 * len(KEYS)), sharex=True)
    for ax, k in zip(axes, KEYS):
        for c in chains:
            ax.plot(c["traj"][k][:n], lw=0.5, color=COLORS.get(c["info"]["start"], "k"), alpha=0.7)
        ax.axvline(E, color="k", ls="--", lw=0.8)
        ax.set_ylabel(k, fontsize=8)
    axes[-1].set_xlabel("MC cycle (dashed: end of fixed warm-up)")
    handles = [plt.Line2D([], [], color=v, label=s) for s, v in COLORS.items() if s in starts]
    axes[0].legend(handles=handles, fontsize=7, ncol=4)
    fig.tight_layout()
    fig.savefig(run / "traces.png", dpi=130)
    plt.close(fig)

    fig, axes = plt.subplots(2, 5, figsize=(16, 6))
    for ax, k in zip(axes.ravel(), KEYS):
        for s in starts:
            vals = np.concatenate([c["traj"][k][E:n] for c in chains if c["info"]["start"] == s])
            ax.hist(vals, bins=40, density=True, histtype="step", color=COLORS.get(s), label=s)
        ax.set_title(k, fontsize=9)
    axes[0, 0].legend(fontsize=6)
    fig.suptitle("post-warm-up distributions by initial structure (should coincide)")
    fig.tight_layout()
    fig.savefig(run / "distributions_by_start.png", dpi=130)
    plt.close(fig)

    fig, axes = plt.subplots(2, 5, figsize=(16, 6))
    for ax, k in zip(axes.ravel(), KEYS):
        for c in chains:
            rho = autocorrelation(c["traj"][k][E:n])
            ax.plot(rho[: min(400, len(rho))], lw=0.6, color=COLORS.get(c["info"]["start"], "k"))
        ax.axhline(0, color="k", lw=0.5)
        ax.set_title(k, fontsize=9)
        ax.set_xlabel("lag (cycles)", fontsize=8)
    fig.tight_layout()
    fig.savefig(run / "autocorrelation.png", dpi=130)
    plt.close(fig)

    # ------------------------------------------------------------------ report
    cols = ["variable", "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean", "mcse_batch",
            "tau_int_max", "geweke_max_abs", "status"]
    key_pooled = [r for r in pooled if r["variable"] in KEYS]
    text = [f"# {a.title or run.name}\n",
            f"{len(chains)} chains x {n} cycles; fixed warm-up {E} cycles; "
            f"{n - E} retained draws per chain, no thinning.\n",
            "## Pooled diagnostics (all starts)\n", table(key_pooled, cols),
            "\n## Per-start diagnostics\n"]
    for s in starts:
        text += [f"\n### {s}\n", table([r for r in per_start if r["scope"] == s and r["variable"] in KEYS],
                                      cols)]
    text += ["\n## Warm-up check: MSER-5 truncation point per chain (cycles)\n",
             table(burn, list(burn[0])),
             "\n## Stationarity: first vs second half of the kept window (chain-averaged)\n",
             table(windows, ["start", "variable", "first_half", "second_half", "change"]),
             "\n## Acceptance, constraint rejections, energy drift\n",
             table(acc_rows, list(acc_rows[0])),
             "\nFigures: traces.png, distributions_by_start.png, autocorrelation.png\n"]
    (run / "REPORT.md").write_text("".join(text), encoding="utf-8")
    summary = {r["variable"]: dict(mean=r["mean"], mcse=r["mcse_mean"], rhat=r["rhat"],
                                   ess_bulk=r["ess_bulk"], status=r["status"]) for r in pooled}
    (run / "summary.json").write_text(json.dumps(dict(equil=E, cycles=n, chains=len(chains),
                                                      pooled=summary), indent=2, default=float),
                                      encoding="utf-8")
    print("".join(text[:4]))


if __name__ == "__main__":
    main()
