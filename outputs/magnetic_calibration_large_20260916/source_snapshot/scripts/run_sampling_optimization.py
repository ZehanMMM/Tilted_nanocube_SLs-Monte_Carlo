"""Checkpointed magnetic calibration and joint multi-start validation.

Conditional runs compare seeds only within the same frozen geometry.
Joint runs monitor finite-run mixing of the existing unconfined model.
Neither changes the interaction energies or adds a spatial constraint.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import random
import shutil
import subprocess
import time

from run_sampling_experiment import ROOT, load_cells, verify_result, write_json

FACTORIAL_STARTS = ["compact_aligned", "compact_tilted40",
                    "expanded_aligned", "expanded_tilted"]
MAGNETIC_VARIABLES = ["energy_kBT", "local_beta_deg", "magnetization", "abs_muB"]


def environment(config):
    ns = load_cells(config["cycles"], config["equil"])
    ns["GLOBAL_COTILT_STEP_DEG"] = config["cotilt_step_deg"]
    ns["DIP_STEP_RAD"] = config["dip_step_rad"]
    ns["GLOBAL_SCALE_MOVES_PER_CYCLE"] = config["scale_moves"]
    ns["DIAGNOSTIC_MCSE_TARGETS"].update(local_beta_deg=0.25, magnetization=0.002)
    return ns


def worker(job, output, config):
    folder = Path(output) / job["name"]
    folder.mkdir(exist_ok=True)
    with (folder / "run.log").open("w", encoding="utf-8", buffering=1) as log:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            with contextlib.redirect_stdout(io.StringIO()):
                ns = environment(config)
            np, model = ns["np"], ns["model"]
            start_id = list(ns["INITIAL_STRUCTURES"]).index(job["start"])
            # Separate streams for calibration and joint validation, including
            # when seed labels coincide. Candidates share initial dipoles.
            phase = 2 if config["mode"] == "conditional" else 3
            init_entropy = [job["seed"], start_id, phase, 0]
            sample_entropy = [job["seed"], start_id, phase, 1]
            state = ns["make_initial_state"](np.random.SeedSequence(init_entropy), job["start"])
            if config["mode"] == "conditional":
                mode = config["seeds"].index(job["seed"]) % 4
                if mode == 1:
                    state["mu"][:] = model.Bhat
                elif mode == 2:
                    state["mu"][:] = -model.Bhat
                elif mode == 3:
                    state["mu"][:] = np.einsum("nij,j->ni", state["Q"], model.e111)
                init_mode = ["random", "field", "anti_field", "body_easy"][mode]
            else:
                init_mode = "random"
            options = ns["sampling_options"]()
            options.update(n_cycles=config["cycles"], n_equil=config["equil"],
                           n_mag_per_cycle=job["sweeps"],
                           record_magnetic_detail=config["mode"] == "conditional")
            if config["mode"] == "conditional":
                options.update(move_positions=False, move_orientations=False)
            rng = np.random.default_rng(np.random.SeedSequence(sample_entropy))
            initial_metrics = dict(
                body_tilt=model.body_axis_metrics(state["Q"])["coherent_tilt_deg"],
                sl_pca_tilt=model.sl_tilt_metrics(state["pos"])["pca_tilt_deg"],
                min_gap_nm=model.surface_gap_stats(state["pos"], state["Q"])["min_gap_nm"],
            )
            print("Actual job:", job, "magnetic initialization:", init_mode)
            print("Initial geometry:", initial_metrics)
            t0 = time.perf_counter()
            result = ns["run_local_collective_mc"](
                model, state, ns["A_NM"], **options, rng=rng, progress_label=job["name"])
            elapsed = time.perf_counter() - t0
            checks = verify_result(ns, result)
            if config["mode"] == "conditional":
                for key in ["pos", "Q"]:
                    np.testing.assert_array_equal(result[key], state[key])
                checks["frozen_geometry_unchanged"] = True
            np.testing.assert_allclose(np.linalg.norm(result["mu"], axis=1), 1, atol=1e-11)
            np.testing.assert_allclose(result["Q"] @ result["Q"].transpose(0, 2, 1),
                                       np.tile(np.eye(3), (27, 1, 1)), atol=1e-11)
            arrays = {k: v for k, v in result.items() if isinstance(v, np.ndarray)}
            scalars = {k: v for k, v in result.items() if not isinstance(v, np.ndarray)}
            np.savez_compressed(folder / "result.npz", **arrays)
            write_json(folder / "result.json", {
                "job": job, "sampling": options, "elapsed_sec": elapsed, "scalars": scalars,
                "init_entropy": init_entropy, "sample_entropy": sample_entropy,
                "dipole_initialization": init_mode, "final_rng_state": rng.bit_generator.state,
                "initial_geometry": initial_metrics,
            })
            write_json(folder / "complete.json", checks)
    return job["name"]


def aggregate(jobs, output, config, wall_sec):
    with contextlib.redirect_stdout(io.StringIO()):
        ns = environment(config)
    np = ns["np"]
    all_diags, summaries, chain_rows, window_rows = [], [], [], []
    for sweeps in config["sweeps"]:
        selected = [job for job in jobs if job["sweeps"] == sweeps]
        results, infos = [], []
        for job in selected:
            folder = output / job["name"]
            with np.load(folder / "result.npz") as data:
                results.append({key: data[key].copy() for key in data.files})
            info = json.loads((folder / "result.json").read_text(encoding="utf-8"))
            infos.append(info)
            stats = info["scalars"]
            chain_rows.append(dict(**job, elapsed_sec=info["elapsed_sec"],
                                   init_mode=info["dipole_initialization"],
                                   **{k: stats[k] for k in ["E_mean", "body_tilt_mean",
                                       "sl_pca_tilt_mean", "rg_nm_mean", "acc_dip", "acc_cotilt", "acc_scale"]}))
            for lo, hi in [(config["equil"], (config["equil"] + config["cycles"]) // 2),
                           ((config["equil"] + config["cycles"]) // 2, config["cycles"])]:
                window_rows.append(dict(**job, cycle_first=lo + 1, cycle_last=hi,
                    **{name: float(results[-1][key][lo:hi].mean())
                       for name, key in ns["DIAGNOSTIC_TRAJECTORIES"].items()}))
        groups = list(config["starts"])
        if config["mode"] == "joint":
            groups = ["all_starts"] + groups
        candidate_diags = []
        for group in groups:
            ids = [i for i, job in enumerate(selected)
                   if group == "all_starts" or job["start"] == group]
            names = MAGNETIC_VARIABLES if config["mode"] == "conditional" else list(ns["DIAGNOSTIC_TRAJECTORIES"])
            draws = {name: np.stack([results[i][ns["DIAGNOSTIC_TRAJECTORIES"][name]][config["equil"]:]
                                     for i in ids]) for name in names}
            if config["mode"] == "conditional":
                for field, prefix in [("traj_particle_muB", "particle_muB"),
                                      ("traj_particle_beta", "particle_beta")]:
                    for particle in range(27):
                        draws[f"{prefix}_{particle}"] = np.stack([
                            results[i][field][config["equil"]:, particle] for i in ids])
            elapsed = sum(infos[i]["elapsed_sec"] for i in ids)
            diags = ns["diagnose_chains"](draws, elapsed)
            all_diags.extend(dict(sweeps=sweeps, scope=group, **row) for row in diags)
            candidate_diags.extend(diags)
        rates = [d["ess_bulk_per_sec"] for d in candidate_diags]
        summaries.append(dict(sweeps=sweeps,
            elapsed_sec=sum(i["elapsed_sec"] for i in infos),
            all_checks_passed=all(d["status"] == "checks_passed" for d in candidate_diags),
            failed_checks=sum(d["status"] != "checks_passed" for d in candidate_diags),
            max_rhat=max(d["rhat_rank"] for d in candidate_diags),
            min_bulk_ess=min(d["ess_bulk"] for d in candidate_diags),
            min_bulk_ess_per_sec=min(rates)))
    eligible = [s for s in summaries if s["all_checks_passed"] and np.isfinite(s["min_bulk_ess_per_sec"])]
    recommendation = max(eligible, key=lambda s: s["min_bulk_ess_per_sec"])["sweeps"] if eligible else None
    ns["write_csv_rows"](output / "diagnostics.csv", all_diags)
    ns["write_csv_rows"](output / "candidates.csv", summaries)
    ns["write_csv_rows"](output / "chains.csv", chain_rows)
    ns["write_csv_rows"](output / "windows.csv", window_rows)
    write_json(output / "finished.json", {
        "completed_utc": datetime.now(timezone.utc).isoformat(), "jobs": len(jobs),
        "wall_sec_this_invocation": wall_sec, "recommended_sweeps": recommendation,
        "recommendation_scope": "fixed-geometry magnetic calibration only" if config["mode"] == "conditional"
                                else "finite-run joint diagnostics; no unrestricted equilibrium claim",
    })
    print("Candidate results:", summaries, flush=True)
    print("Recommended sweeps:", recommendation, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["conditional", "joint"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cycles", type=int, default=2000)
    parser.add_argument("--equil", type=int, default=500)
    parser.add_argument("--sweeps", nargs="+", type=int, default=[1, 5, 10])
    parser.add_argument("--starts", nargs="+", default=FACTORIAL_STARTS)
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 12, 13, 14])
    parser.add_argument("--cotilt-step-deg", type=float, default=0.75)
    parser.add_argument("--dip-step-rad", type=float, default=0.30)
    parser.add_argument("--scale-moves", type=int, default=1)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.equil < args.cycles - 1 or args.workers < 1:
        parser.error("Require at least two post-warm-up draws and positive workers.")
    for values, lower in [(args.sweeps, 1), (args.seeds, 0)]:
        if len(set(values)) != len(values) or min(values) < lower:
            parser.error("Sweeps/seeds must be unique and in range.")
    config = {k: v for k, v in vars(args).items() if k not in ["output", "workers", "resume"]}
    with contextlib.redirect_stdout(io.StringIO()):
        ns = environment(config)
    if len(set(args.starts)) != len(args.starts) or any(s not in ns["INITIAL_STRUCTURES"] for s in args.starts):
        parser.error("Starts must be unique known structures.")
    jobs = [dict(start=start, seed=seed, sweeps=sweeps, name=f"m{sweeps}_{start}_s{seed}")
            for sweeps in args.sweeps for start in args.starts for seed in args.seeds]
    files = sorted((ROOT / "_v903cells").glob("cell*.*")) + [Path(__file__), ROOT / "scripts/run_sampling_experiment.py"]
    config["source_sha256"] = {p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    output = args.output.resolve()
    if args.resume:
        if json.loads((output / "manifest.json").read_text())["config"] != config:
            raise RuntimeError("Resume configuration/source mismatch.")
    else:
        output.mkdir(parents=True, exist_ok=False)
        write_json(output / "manifest.json", dict(config=config, workers=args.workers,
            started_utc=datetime.now(timezone.utc).isoformat(),
            base_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            structures={s: ns["INITIAL_STRUCTURES"][s] for s in args.starts},
            physical_parameters=ns["model"].phys,
            model_settings={k: ns[k] for k in ["T_K", "L_NM", "A_NM", "ALPHA_DEG", "ANIS_MODEL"]},
            interpretation="No spatial confinement. Conditional magnetic targets are proper. Joint runs do not establish unrestricted cluster equilibrium."))
        for source in files:
            target = output / "source_snapshot" / source.relative_to(ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    pending = [job for job in jobs if not (output / job["name"] / "complete.json").exists()]
    random.Random(903).shuffle(pending)
    print(f"Starting {len(pending)} {args.mode} jobs with {args.workers} workers.", flush=True)
    t0 = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(worker, job, str(output), config) for job in pending]
        for f in as_completed(futures):
            print(f"FINISHED {f.result()}, elapsed {(time.perf_counter() - t0) / 60:.1f} min", flush=True)
    aggregate(jobs, output, config, time.perf_counter() - t0)


if __name__ == "__main__":
    main()
