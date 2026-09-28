"""Run a set of independent chains in parallel, resumably.

    python scripts/run_chains.py CONFIG.json OUTPUT_DIR [--workers 10]

CONFIG keys (all optional except starts/seeds/cycles):
    starts          list of names in model.INITIAL_STRUCTURES
    seeds           list of ints
    cycles          int
    moves           MoveSettings overrides
    physics         PhysicalParameters overrides (vdw_settings as a dict)
    ensemble        ClusterEnsemble overrides
    init_from       optional folder of a previous run; its final states seed these chains
    init_from_seed  which seed's final state to take from init_from (default: same seed)
    init_seed_offset  alternatively take seed - offset, so fresh sampling seeds can
                    continue chosen source chains (e.g. seeds 21, 22 <- sources 1, 2)
    snapshot_every  save full configurations every k cycles (default 100)
    drift_every     recompute the energy from scratch every k cycles (default 250)
    schedule        'cluster' (default) or 'lattice' (rigid lattice: bodies, lattice
                    orientation and moments only)
    structures      custom starting structures {name: initial_state kwargs};
                    default model.INITIAL_STRUCTURES
    anneal          {"t_start": .., "t_end": .., "adapt_every": 20}: simulated
                    annealing along a geometric kT schedule instead of sampling;
                    the best state found is saved as best_state.npz

Each chain gets two independent streams from SeedSequence([seed, start_id,
stream]): stream 0 for the initial moments, stream 1 for the sampler.  A
chain is complete only when complete.json exists, so a killed run can be
restarted with the same command and only unfinished chains are rerun.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402

from v904.anneal import anneal, geometric_schedule  # noqa: E402
from v904.model import (ClusterEnsemble, Hamiltonian, INITIAL_STRUCTURES,  # noqa: E402
                        PhysicalParameters, initial_state)
from v904.sampler import MoveSettings, run_chain  # noqa: E402
from v904.vdw import HamakerSettings  # noqa: E402


def structures(config):
    return config.get("structures") or INITIAL_STRUCTURES


def build(config):
    phys = dict(config.get("physics", {}))
    if "vdw_settings" in phys:
        phys["vdw_settings"] = HamakerSettings(**phys["vdw_settings"])
    ham = Hamiltonian(PhysicalParameters(**phys), ClusterEnsemble(**config.get("ensemble", {})))
    moves = replace(MoveSettings(), **config.get("moves", {}))
    return ham, moves


def source_hash():
    h = hashlib.sha256()
    for path in sorted((ROOT / "v904").glob("*.py")):
        h.update(path.read_bytes())
    return h.hexdigest()


def run_one(config, start, seed, folder):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    ham, moves = build(config)
    start_ids = {name: k for k, name in enumerate(structures(config))}
    seq = np.random.SeedSequence([int(seed), start_ids[start]])
    init_seq, sample_seq = seq.spawn(2)
    if config.get("init_from"):
        src_seed = config.get("init_from_seed", seed - config.get("init_seed_offset", 0))
        src = Path(config["init_from"]) / f"{start}_s{src_seed}" / "final_state.npz"
        with np.load(src) as d:
            state = dict(pos=d["pos"], Q=d["Q"], mu=d["mu"])
            if "G" in d.files:
                state["G"] = d["G"]
    else:
        state = initial_state(**structures(config)[start],
                              dipole_seed=np.random.default_rng(init_seq).integers(2 ** 32))
    log = (folder / "progress.txt").open("w", encoding="utf-8", buffering=1)

    def progress(done, total, elapsed):
        if done == 1 or done % 50 == 0 or done == total:
            log.write(f"{done}/{total} cycles, {elapsed / 60:.1f} min, "
                      f"ETA {(total - done) * elapsed / done / 60:.1f} min\n")

    schedule = config.get("schedule", "cluster")
    extra = {}
    if config.get("anneal"):
        a = config["anneal"]
        factors = geometric_schedule(config["cycles"], a["t_start"], a["t_end"])
        out = anneal(ham, state, moves, np.random.default_rng(sample_seq), factors,
                     adapt_every=a.get("adapt_every", 20), progress=progress,
                     lattice=schedule == "lattice")
        out.update(snapshots=[], drifts=[(config["cycles"], out["drift_kBT"])],
                   acceptance=out["chain"].acceptance(), tries=dict(out["chain"].tries),
                   blocked=dict(out["chain"].blocked), sq_jump=dict(out["chain"].sq_jump),
                   moves=out["final_moves"])
        best = out["best"]
        np.savez_compressed(folder / "best_state.npz", **{k: v for k, v in best.items()})
        extra = dict(best_energy_kBT=out["best_energy_kBT"], best_cycle=best["cycle"],
                     best_allowed=bool(ham.allowed(best["pos"], best["Q"])))
    else:
        out = run_chain(ham, state, moves, np.random.default_rng(sample_seq), config["cycles"],
                        progress=progress, snapshot_every=config.get("snapshot_every", 100),
                        check_drift_every=config.get("drift_every", 250), schedule=schedule)
    chain = out["chain"]
    drifts = np.array(out["drifts"]) if out["drifts"] else np.zeros((0, 2))
    np.savez_compressed(folder / "trajectory.npz", **out["traj"])
    np.savez_compressed(folder / "final_state.npz", pos=chain.pos, Q=chain.Q, mu=chain.mu,
                        G=chain.G)
    if out["snapshots"]:
        np.savez_compressed(folder / "snapshots.npz",
                            cycle=np.array([s["cycle"] for s in out["snapshots"]]),
                            pos=np.array([s["pos"] for s in out["snapshots"]]),
                            Q=np.array([s["Q"] for s in out["snapshots"]]),
                            mu=np.array([s["mu"] for s in out["snapshots"]]),
                            G=np.array([s.get("G", np.eye(3)) for s in out["snapshots"]]))
    max_drift = float(np.max(np.abs(drifts[:, 1]))) if len(drifts) else 0.0
    info = dict(start=start, seed=int(seed), cycles=config["cycles"],
                elapsed_sec=out["elapsed_sec"], acceptance=out["acceptance"],
                tries=out["tries"], blocked=out["blocked"], sq_jump=out["sq_jump"],
                moves=out["moves"], hamiltonian=ham.describe(),
                seed_entropy=[int(x) for x in np.atleast_1d(seq.entropy)],
                spawn_key=[int(x) for x in seq.spawn_key],
                max_abs_energy_drift_kBT=max_drift,
                quadrature_max_level_hits=int(ham.hamaker.max_level_hits),
                exact_core_distance_checks=int(ham.exact_distance_checks),
                final_allowed=bool(ham.allowed(chain.pos, chain.Q)),
                all_finite=bool(all(np.all(np.isfinite(v)) for v in out["traj"].values())),
                schedule=schedule, **extra)
    (folder / "chain.json").write_text(json.dumps(info, indent=2, default=str), encoding="utf-8")
    if (max_drift > 1e-6 or not info["final_allowed"] or not info["all_finite"]
            or info["quadrature_max_level_hits"] > 0 or not info.get("best_allowed", True)):
        raise RuntimeError(f"verification failed for {folder.name}: {info}")
    (folder / "complete.json").write_text(json.dumps(dict(max_drift=max_drift)), encoding="utf-8")
    log.close()
    return folder.name


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config")
    parser.add_argument("output")
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()
    config = json.loads(Path(args.config).read_text(encoding="utf-8"))
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    ham, moves = build(config)
    manifest = dict(config=config, source_sha256=source_hash(), moves=asdict(moves),
                    hamiltonian=ham.describe(), started=time.strftime("%Y-%m-%d %H:%M:%S"))
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, default=str),
                                       encoding="utf-8")
    jobs = [(s, k) for s in config["starts"] for k in config["seeds"]
            if not (out / f"{s}_s{k}" / "complete.json").exists()]
    print(f"{len(jobs)} chains to run in {out}", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_one, config, s, k, out / f"{s}_s{k}"): (s, k) for s, k in jobs}
        for f in as_completed(futures):
            print(f"done {f.result()} after {(time.time() - t0) / 60:.1f} min", flush=True)
    (out / "finished.json").write_text(json.dumps(dict(wall_min=(time.time() - t0) / 60)),
                                       encoding="utf-8")


if __name__ == "__main__":
    main()
