"""Run the notebook's full single chain and multi-start ensemble with checkpoints."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import contextlib
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import time

# Avoid nested BLAS pools when independent chains use separate processes.
for name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ[name] = "1"
os.environ["MPLBACKEND"] = "Agg"
ROOT = Path(__file__).resolve().parents[1]


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False), encoding="utf-8")


def load_cells(cycles, equil):
    ns = {}
    for name in ["cell1.py", "cell2.py", "cell2a.py", "cell3.py"]:
        path = ROOT / "_v903cells" / name
        exec(compile(path.read_text(encoding="utf-8"), str(path), "exec"), ns)
    ns.update(N_CYCLES=cycles, N_EQUIL=equil, MULTI_CYCLES=cycles, MULTI_EQUIL=equil,
              SHOW_FIGURES=False, SHOW_PROGRESS_BARS=True, PROGRESS_EVERY_CYCLES=100,
              RUN_MAGNETIC_PILOT=False)
    return ns


def verify_result(ns, result):
    np, model = ns["np"], ns["model"]
    for key in ns["DIAGNOSTIC_TRAJECTORIES"].values():
        if len(result[key]) != ns["N_CYCLES"] or not np.all(np.isfinite(result[key])):
            raise RuntimeError(f"Invalid trajectory: {key}")
    exact = model.energy_full_cluster(result["pos"], result["Q"], result["mu"],
                                     a_nm=ns["A_NM"], anis_model=ns["ANIS_MODEL"])["Total"]
    drift = float(exact / model.kT - result["traj_E"][-1])
    if abs(drift) > 1e-6:
        raise RuntimeError(f"Energy tracking drift: {drift} kBT")
    if np.min(result["traj_gap_min_nm"]) <= 0:
        raise RuntimeError("A sampled state failed the existing positive-gap check.")
    return {"full_minus_tracked_energy_kBT": drift,
            "min_recorded_gap_nm": float(np.min(result["traj_gap_min_nm"])),
            "all_trajectories_finite": True}


def run_job(job, output, cycles, equil):
    folder = Path(output) / job["name"]
    folder.mkdir(exist_ok=True)
    # A success marker is written last. Failed partial jobs are safe to rerun.
    with (folder / "run.log").open("w", encoding="utf-8", buffering=1) as log:
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
            ns = load_cells(cycles, equil)
            np = ns["np"]
            if job["kind"] == "single":
                ns.update(
                    FULL_RUN_FIG_PDF=str(folder / f"V903_Fig1_N27_{cycles}cycles.pdf"),
                    FULL_RUN_REPORT_PDF=str(folder / "V903_FullRun_Report.pdf"),
                )
                source = ROOT / "_v903cells" / "cell5.py"
                exec(compile(source.read_text(encoding="utf-8"), str(source), "exec"), ns)
                result = ns["res"]
                info = {
                    "kind": "single", "init_seed": ns["INIT_SEED"],
                    "sample_seed": ns["FULL_RUN_SEED"], "n_cycles": cycles, "n_equil": equil,
                    "sampling": ns["sampling_options"](), "elapsed_sec": ns["elapsed"],
                    "physical_parameters": dict(ns["model"].phys), "temperature_K": ns["T_K"],
                    "L_nm": ns["L_NM"], "A_nm": ns["A_NM"], "alpha_deg": ns["ALPHA_DEG"],
                }
            else:
                bundle = ns["run_chain_ensemble"]([job["start"]], [job["seed"]], cycles, equil)
                result = bundle["results"][0]
                info = {"kind": "multi", "row": bundle["rows"][0],
                        "metadata": bundle["metadata"]}
            verification = verify_result(ns, result)
            arrays = {k: v for k, v in result.items() if isinstance(v, np.ndarray)}
            scalars = {k: v for k, v in result.items() if not isinstance(v, np.ndarray)}
            np.savez_compressed(folder / "result.npz", **arrays)
            write_json(folder / "result.json", {"scalars": scalars, "run": info})
            write_json(folder / "complete.json", verification)
            print("Completed and verified:", verification)
    return job["name"]


def assemble(jobs, output, cycles, equil, wall_sec):
    with contextlib.redirect_stdout(io.StringIO()):
        ns = load_cells(cycles, equil)
    np = ns["np"]
    results, rows, worker_metadata = [], [], []
    for job in jobs:
        if job["kind"] != "multi":
            continue
        folder = output / job["name"]
        info = json.loads((folder / "result.json").read_text(encoding="utf-8"))
        with np.load(folder / "result.npz") as data:
            result = {key: data[key].copy() for key in data.files}
        result.update(info["scalars"])
        results.append(result)
        row = info["run"]["row"]
        row["chain"] = len(rows)
        rows.append(row)
        worker_metadata.append(info["run"]["metadata"])
    draws = {name: np.stack([r[key][equil:] for r in results])
             for name, key in ns["DIAGNOSTIC_TRAJECTORIES"].items()}
    elapsed = sum(row["elapsed_sec"] for row in rows)
    grouped = {}
    for start in ns["MULTI_STARTS"]:
        ids = [i for i, row in enumerate(rows) if row["start"] == start]
        grouped[start] = ns["diagnose_chains"](
            {name: x[ids] for name, x in draws.items()},
            sum(rows[i]["elapsed_sec"] for i in ids))
    metadata = worker_metadata[0]
    metadata.update(
        starts={key: ns["INITIAL_STRUCTURES"][key] for key in ns["MULTI_STARTS"]},
        seeds=ns["MULTI_SEEDS"], elapsed_sec=elapsed, execution_wall_sec=wall_sec,
        efficiency_time_basis="sum of per-chain sampling seconds, including warm-up",
    )
    bundle = dict(results=results, rows=rows, draws=draws,
                  diagnostics=ns["diagnose_chains"](draws, elapsed),
                  grouped_diagnostics=grouped, metadata=metadata)
    ns["export_chain_ensemble"](bundle, str(output / "V903_MultiChain"),
                               report_pdf=str(output / "V903_MultiChain_Report.pdf"))
    print("Pooled diagnostics:", flush=True)
    for line in ns["rows_to_text_table"](bundle["diagnostics"],
            ["variable", "rhat_rank", "ess_bulk", "ess_tail", "mcse_mean", "status"]):
        print(line, flush=True)
    write_json(output / "finished.json", {
        "completed_utc": datetime.now(timezone.utc).isoformat(),
        "jobs": len(jobs), "cycles_per_chain": cycles, "warmup_per_chain": equil,
        "wall_sec_this_invocation": wall_sec,
        "all_pooled_diagnostic_checks_passed":
            all(row["status"] == "checks_passed" for row in bundle["diagnostics"]),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--cycles", type=int, default=2000)
    parser.add_argument("--equil", type=int, default=500)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if not 0 <= args.equil < args.cycles or args.workers < 1:
        parser.error("Require 0 <= equil < cycles and workers >= 1.")
    with contextlib.redirect_stdout(io.StringIO()):
        ns = load_cells(args.cycles, args.equil)
    jobs = [{"kind": "single", "name": "single_seed1"}] + [
        {"kind": "multi", "start": start, "seed": seed, "name": f"{start}_seed{seed}"}
        for start in ns["MULTI_STARTS"] for seed in ns["MULTI_SEEDS"]
    ]
    source_files = sorted((ROOT / "_v903cells").glob("cell*.*")) + [Path(__file__)]
    config = {
        "cycles": args.cycles, "equil": args.equil, "jobs": jobs,
        "source_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in source_files},
    }
    output = args.output.resolve()
    if args.resume:
        manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
        if manifest["config"] != config:
            raise RuntimeError("Resume refused: configuration or source files have changed.")
    else:
        output.mkdir(parents=True, exist_ok=False)
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
        write_json(output / "manifest.json", {
            "started_utc": datetime.now(timezone.utc).isoformat(), "base_commit": revision,
            "workers": args.workers, "config": config,
        })
    pending = [job for job in jobs if not (output / job["name"] / "complete.json").exists()]
    print(f"Running {len(pending)} jobs, {args.workers} workers, {args.cycles} cycles each.",
          flush=True)
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_job, job, str(output), args.cycles, args.equil)
                   for job in pending]
        for future in as_completed(futures):
            print(f"FINISHED {future.result()} after {(time.perf_counter()-start)/60:.1f} min",
                  flush=True)
    assemble(jobs, output, args.cycles, args.equil, time.perf_counter() - start)
    print(f"All results written to {output}", flush=True)


if __name__ == "__main__":
    main()
