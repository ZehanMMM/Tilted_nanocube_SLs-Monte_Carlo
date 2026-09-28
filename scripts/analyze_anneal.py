"""Summarise annealing runs: best states, agreement between restarts, local test.

    python scripts/analyze_anneal.py CONFIG.json RUN_DIR

For every restart: best total energy and its structure.  The evidence for a
global minimum is agreement of independent restarts, not any single run.
The overall best state is then tested for local optimality: 400 tiny
random perturbations (single-cube rotations, single-moment rotations,
lattice rotation, co-tilt, twist), each evaluated from scratch.  At a local
minimum none of them should lower the energy by more than round-off.
"""
from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from run_chains import build  # noqa: E402
from v904.sampler import Chain, observe  # noqa: E402

SHOW = ["energy_kBT", "sl_axis_tilt_deg", "body_tilt_deg", "body_rotation_deg",
        "body_rotation_max_deg", "body_lattice_angle_deg", "magnetization", "local_beta_deg",
        "zeeman_kBT", "anisotropy_kBT", "dipole_kBT", "vdw_kBT"]


def local_test(ham, moves, state, rng, n=400, step_deg=0.2):
    tiny = replace(moves, rot_step_deg=step_deg, latrot_step_deg=step_deg,
                   cotilt_step_deg=step_deg, gamma_step_deg=step_deg,
                   dip_step_rad=np.radians(step_deg), small_prob=0.0)
    chain = Chain(ham, state, tiny, rng)
    chain.temperature_factor = 1e-12            # accept only downhill moves
    E0 = chain.energy_terms()["Total"]
    kinds = ["rot", "dip", "latrot", "cotilt", "gamma"]
    drops = []
    for t in range(n):
        kind = kinds[t % len(kinds)]
        before = chain.energy_terms()["Total"]
        if kind == "rot":
            chain.rotate(int(rng.integers(chain.N)))
        elif kind == "dip":
            chain.dipole(int(rng.integers(chain.N)))
        elif kind == "latrot":
            chain.lattice_rotate()
        elif kind == "cotilt":
            chain.cotilt()
        else:
            chain.gamma()
        after = chain.energy_terms()["Total"]
        drops.append((before - after) / ham.kT)
    drops = np.array(drops)
    return dict(attempts=n, step_deg=step_deg, accepted_downhill=int(np.sum(drops > 0)),
                total_drop_kBT=float((E0 - chain.energy_terms()["Total"]) / ham.kT),
                largest_single_drop_kBT=float(drops.max()), drift_kBT=chain.drift_kBT())


def main():
    config = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    run = Path(sys.argv[2])
    ham, moves = build(config)
    rows, trajs = [], {}
    for folder in sorted(run.iterdir()):
        if not (folder / "complete.json").exists():
            continue
        info = json.loads((folder / "chain.json").read_text(encoding="utf-8"))
        with np.load(folder / "best_state.npz") as d:
            best = {k: d[k] for k in d.files}
        chain = Chain(ham, dict(pos=best["pos"], Q=best["Q"], mu=best["mu"], G=best["G"]),
                      moves, np.random.default_rng(0))
        obs = observe(chain, lattice=True)
        rows.append(dict(run=folder.name, start=info["start"], best_cycle=int(best["cycle"]),
                         **{k: obs[k] for k in SHOW}))
        with np.load(folder / "trajectory.npz") as d:
            trajs[folder.name] = {k: d[k].copy() for k in ("energy_kBT", "sl_axis_tilt_deg",
                                                          "temperature_factor")}
    rows.sort(key=lambda r: r["energy_kBT"])
    E = np.array([r["energy_kBT"] for r in rows])
    best_row = rows[0]
    with np.load(run / best_row["run"] / "best_state.npz") as d:
        best = {k: d[k] for k in d.files}
    lt = local_test(ham, moves, dict(pos=best["pos"], Q=best["Q"], mu=best["mu"], G=best["G"]),
                    np.random.default_rng(1))
    summary = dict(restarts=len(rows), best=best_row, energy_spread_kBT=float(E.max() - E.min()),
                   energies_kBT=E.tolist(), local_test=lt)
    (run / "anneal_summary.json").write_text(json.dumps(dict(summary=summary, rows=rows), indent=2),
                                             encoding="utf-8")
    cols = ["run", "best_cycle"] + SHOW
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        lines.append("| " + " | ".join(r[c] if isinstance(r[c], str) else f"{r[c]:.4g}"
                                       for c in cols) + " |")
    text = ("# Simulated annealing: best state of every restart (sorted by energy)\n\n"
            + "\n".join(lines) + "\n\n"
            + f"Energy spread over restarts: {summary['energy_spread_kBT']:.3f} kBT\n\n"
            + "Local test of the best state (tiny downhill-only perturbations): "
            + json.dumps(lt) + "\n")
    (run / "ANNEAL_REPORT.md").write_text(text, encoding="utf-8")
    print(text)

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for name, t in trajs.items():
        axes[0].plot(t["energy_kBT"], lw=0.7, label=name)
        axes[1].plot(t["sl_axis_tilt_deg"], lw=0.7)
        axes[2].semilogy(t["temperature_factor"], lw=0.7)
    axes[0].set_ylabel("U / kBT")
    axes[1].set_ylabel("SL axis tilt / deg")
    axes[2].set_ylabel("kT_anneal / kT")
    axes[2].set_xlabel("annealing cycle")
    axes[0].legend(fontsize=6, ncol=2)
    fig.tight_layout()
    fig.savefig(run / "anneal_traces.png", dpi=140)


if __name__ == "__main__":
    main()
