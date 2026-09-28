"""Stability tests: compare variant runs with the reference, observable by observable.

    python scripts/compare_runs.py OUT.md REF_DIR:EQUIL VARIANT_DIR:EQUIL [...]

For each observable: pooled post-warm-up mean and ArviZ MCSE of every run,
and z = (mean_variant - mean_ref) / sqrt(MCSE_ref^2 + MCSE_variant^2).

How to read it depends on what the variant changes:

* sampler-only changes (step sizes)  leave the target unchanged, so every
  |z| should look like a draw from N(0, 1); a systematic |z| >> 2 means one
  of the two runs has not converged (or a kernel is wrong).
* target changes (cluster definition, vdW model) are MEANT to move the
  answer; z then measures how strongly the result depends on that choice.

z is only meaningful when both runs pass their own convergence checks.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from analyze import key_list, load  # noqa: E402
from v904.diagnostics import require_arviz  # noqa: E402


def pooled(run, equil):
    az = require_arviz()
    chains = load(run)
    n = min(len(c["traj"]["energy_kBT"]) for c in chains)
    out = {}
    for k in key_list(chains):
        x = np.stack([c["traj"][k][equil:n] for c in chains])
        out[k] = dict(mean=float(x.mean()), mcse=float(az.mcse(x, method="mean")),
                      rhat=float(az.rhat(x, method="rank")), chains=len(chains), draws=n - equil)
    return out


def main():
    out_md = Path(sys.argv[1])
    specs = [s.rsplit(":", 1) for s in sys.argv[2:]]
    runs = [(Path(p), int(e)) for p, e in specs]
    stats = [pooled(p, e) for p, e in runs]
    ref_name = runs[0][0].name
    lines = [f"# Stability comparison against `{ref_name}`\n\n",
             "z = (variant - reference) / sqrt(MCSE_ref^2 + MCSE_var^2); R-hat of each run in brackets.\n\n",
             "| observable | " + " | ".join(f"{p.name}" for p, _ in runs) + " | "
             + " | ".join(f"z {p.name}" for p, _ in runs[1:]) + " |\n",
             "|---|" + "---|" * (2 * len(runs) - 1) + "\n"]
    for k in stats[0]:
        ref = stats[0][k]
        cells = [f"{s[k]['mean']:.4g} ± {s[k]['mcse']:.2g} [{s[k]['rhat']:.3f}]" for s in stats]
        zs = [(s[k]["mean"] - ref["mean"]) / np.hypot(s[k]["mcse"], ref["mcse"]) for s in stats[1:]]
        lines.append(f"| {k} | " + " | ".join(cells) + " | " + " | ".join(f"{z:+.2f}" for z in zs)
                     + " |\n")
    out_md.write_text("".join(lines), encoding="utf-8")
    print("".join(lines))


if __name__ == "__main__":
    main()
