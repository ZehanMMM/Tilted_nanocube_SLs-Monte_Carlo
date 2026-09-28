"""Step-size ladder for the rigid-lattice schedule, scored by MSJD per attempt.

    python scripts/pilot_tuning_lattice.py configs/lattice_production.json outputs/pilot_tuning_lattice
        [--factors 0.25,0.5,1,2,4,8] [--warm 20] [--cycles 60]

Same logic as scripts/pilot_tuning.py: all tuned steps are scaled together by
f, each move is scored separately by accepted-jump^2 per attempt, the best
f per move is written to tuned_moves.json and then frozen for production;
pilot samples are discarded.  Starts are the configured lattice structures,
the first `warm` cycles of every level are dropped before counting (the
random initial moments relax within a few cycles).
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
sys.path.insert(0, str(ROOT / "scripts"))

import numpy as np  # noqa: E402

from run_chains import build  # noqa: E402
from v904.model import initial_state  # noqa: E402
from v904.sampler import Chain  # noqa: E402

STEP_FIELD = dict(rot="rot_step_deg", latrot="latrot_step_deg", cotilt="cotilt_step_deg",
                  gamma="gamma_step_deg")


def run_level(args):
    config, start, factor, warm, cycles, seed = args
    ham, base = build(config)
    moves = replace(base, **{f: getattr(base, f) * factor for f in STEP_FIELD.values()})
    state = initial_state(**config["structures"][start], dipole_seed=seed)
    chain = Chain(ham, state, moves, np.random.default_rng(seed))
    for _ in range(warm):
        chain.cycle_lattice()
    chain.reset_counters()
    t0 = time.perf_counter()
    for _ in range(cycles):
        chain.cycle_lattice()
    elapsed = time.perf_counter() - t0
    return [dict(start=start, factor=factor, move=m, step=getattr(moves, f),
                 tries=chain.tries[m], acceptance=chain.accepts[m] / max(chain.tries[m], 1),
                 blocked=chain.blocked[m] / max(chain.tries[m], 1),
                 msjd=chain.sq_jump[m] / max(chain.tries[m], 1), sec_per_cycle=elapsed / cycles)
            for m, f in STEP_FIELD.items()]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("config")
    p.add_argument("output")
    p.add_argument("--factors", default="0.25,0.5,1,2,4,8")
    p.add_argument("--warm", type=int, default=20)
    p.add_argument("--cycles", type=int, default=60)
    p.add_argument("--workers", type=int, default=16)
    a = p.parse_args()
    config = json.loads(Path(a.config).read_text(encoding="utf-8"))
    factors = [float(f) for f in a.factors.split(",")]
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=True)
    starts = list(config["structures"])
    jobs = [(config, s, f, a.warm, a.cycles, 8000 + 10 * k + n)
            for n, s in enumerate(starts) for k, f in enumerate(factors)]
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        rows = [r for block in pool.map(run_level, jobs) for r in block]
    with open(out / "ladder.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    _, base = build(config)
    tuned, table = {}, {}
    for move, field in STEP_FIELD.items():
        scores = {f: float(np.mean([r["msjd"] for r in rows if r["move"] == move and r["factor"] == f]))
                  for f in factors}
        acc = {f: float(np.mean([r["acceptance"] for r in rows if r["move"] == move and r["factor"] == f]))
               for f in factors}
        best = max(scores, key=scores.get)
        tuned[field] = getattr(base, field) * best
        table[move] = dict(best_factor=best, msjd=scores, acceptance=acc,
                           edge_of_ladder=best in (factors[0], factors[-1]))
    final = asdict(replace(base, **tuned))
    (out / "tuned_moves.json").write_text(json.dumps(dict(moves=final, ladder=table), indent=2),
                                          encoding="utf-8")
    print(json.dumps(table, indent=2))
    print("tuned:", {k: final[k] for k in STEP_FIELD.values()})
    print("sec/cycle:", np.mean([r["sec_per_cycle"] for r in rows]))


if __name__ == "__main__":
    main()
