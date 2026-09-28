"""Pilot 2: step-size ladder, scored by mean squared jump distance per attempt.

    python scripts/pilot_tuning.py outputs/pilot_relax outputs/pilot_tuning [--cycles 120]
        [--factors 0.25,0.5,1,2,4] [--tag round1]

If the best factor of a move sits on the edge of the ladder, the optimum
has not been bracketed: run again with factors extending past that edge
(a new --tag).  tuned_moves.json always takes, per move, the argmax over
every ladder round found in the output folder.

Starting from the relaxed pilot-1 states (never from the lattice, whose
first cycles are a transient), every step size is scaled together by
f in {1/4, 1/2, 1, 2, 4}.  For each move type separately:

    MSJD = sum(accepted jump^2) / attempts  = acceptance x E[jump^2 | accepted]

For a random-walk proposal MSJD is the lag-1 quantity that ESS per attempt
tracks: too small a step is accepted but goes nowhere, too large a step is
rejected.  Cost per attempt hardly depends on the step, so the argmax of
MSJD is the step to use.  The winner of each move is written to
tuned_moves.json and then FROZEN for production; pilot samples are
discarded.  Magnetic steps are not re-tuned here: the conditional
distribution of the moments given the geometry contains no vdW, so the
V9.03 fixed-geometry calibration (1.60 rad, 10 sweeps) still applies.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from v904.model import Hamiltonian  # noqa: E402
from v904.sampler import Chain, MoveSettings  # noqa: E402

BASE = MoveSettings()
FACTORS = [0.25, 0.5, 1.0, 2.0, 4.0]
STEP_FIELD = dict(trans="trans_step_nm", rot="rot_step_deg", cotilt="cotilt_step_deg",
                  gamma="gamma_step_deg", scale="scale_log_step")


def ladder_moves(f):
    return replace(BASE, **{field: getattr(BASE, field) * f for field in STEP_FIELD.values()})


def run_level(args):
    src, start, factor, cycles, seed, tag = args
    with np.load(Path(src) / f"{start}_s101" / "final_state.npz") as d:
        state = dict(pos=d["pos"], Q=d["Q"], mu=d["mu"])
    moves = ladder_moves(factor)
    chain = Chain(Hamiltonian(), state, moves, np.random.default_rng(seed))
    t0 = time.perf_counter()
    for _ in range(cycles):
        chain.cycle()
    elapsed = time.perf_counter() - t0
    rows = []
    for move, field in STEP_FIELD.items():
        tries = chain.tries[move]
        rows.append(dict(tag=tag, start=start, factor=factor, move=move, step=getattr(moves, field),
                         tries=tries, acceptance=chain.accepts[move] / max(tries, 1),
                         blocked=chain.blocked[move] / max(tries, 1),
                         msjd=chain.sq_jump[move] / max(tries, 1),
                         sec_per_cycle=elapsed / cycles))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("source")
    p.add_argument("output")
    p.add_argument("--cycles", type=int, default=120)
    p.add_argument("--workers", type=int, default=10)
    p.add_argument("--factors", default=",".join(str(f) for f in FACTORS))
    p.add_argument("--tag", default="round1")
    a = p.parse_args()
    factors = [float(f) for f in a.factors.split(",")]
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    starts = ["compact_aligned", "compact_tilted40", "expanded_aligned", "expanded_tilted"]
    base_seed = 7000 + 1000 * len(list(out.glob("ladder_*.csv")))
    jobs = [(a.source, s, f, a.cycles, base_seed + 10 * k + n, a.tag)
            for n, s in enumerate(starts) for k, f in enumerate(factors)]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        new = [r for block in pool.map(run_level, jobs) for r in block]
    with open(out / f"ladder_{a.tag}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(new[0]))
        w.writeheader()
        w.writerows(new)
    rows = []
    for path in sorted(out.glob("ladder_*.csv")):
        with open(path, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                for key in ("factor", "step", "acceptance", "blocked", "msjd", "sec_per_cycle"):
                    r[key] = float(r[key])
                rows.append(r)
    all_factors = sorted({r["factor"] for r in rows})
    tuned, table = {}, {}
    for move, field in STEP_FIELD.items():
        # a factor measured in two rounds is averaged over both (8 chains)
        scores = {f: np.mean([r["msjd"] for r in rows if r["move"] == move and r["factor"] == f])
                  for f in all_factors}
        acc = {f: np.mean([r["acceptance"] for r in rows if r["move"] == move and r["factor"] == f])
               for f in all_factors}
        best = max(scores, key=scores.get)
        tuned[field] = getattr(BASE, field) * best
        table[move] = dict(best_factor=best, msjd=scores, acceptance=acc,
                           edge_of_ladder=best in (all_factors[0], all_factors[-1]))
    final = asdict(replace(BASE, **tuned))
    (out / "tuned_moves.json").write_text(json.dumps(dict(moves=final, ladder=table), indent=2,
                                                     default=float), encoding="utf-8")
    print(json.dumps(table, indent=2, default=float))
    print("tuned:", final)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 5, figsize=(16, 6), sharex=True)
    for c, move in enumerate(STEP_FIELD):
        for s in starts:
            sub = sorted([r for r in rows if r["move"] == move and r["start"] == s],
                         key=lambda r: r["factor"])
            axes[0, c].plot([r["factor"] for r in sub], [r["msjd"] for r in sub], "o", label=s, ms=4)
            axes[1, c].plot([r["factor"] for r in sub], [r["acceptance"] for r in sub], "o", ms=4)
        axes[0, c].axvline(table[move]["best_factor"], color="k", ls=":", lw=1)
        axes[0, c].set_title(f"{move}  (step x f)")
        axes[1, c].set_xlabel("step factor f")
        axes[0, c].set_xscale("log", base=2)
    axes[0, 0].set_ylabel("MSJD per attempt")
    axes[1, 0].set_ylabel("acceptance")
    axes[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "ladder.png", dpi=150)


if __name__ == "__main__":
    main()
