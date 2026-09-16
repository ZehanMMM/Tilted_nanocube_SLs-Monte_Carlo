# Multi-chain diagnostics and optional magnetic-sweep comparison.
import csv
import json
from pathlib import Path

DIAGNOSTIC_TRAJECTORIES = {
    "energy_kBT": "traj_E",
    "local_beta_deg": "traj_beta",
    "body_tilt_deg": "traj_body_tilt",
    "sl_pca_tilt_deg": "traj_sl_pca_tilt",
    "body_sl_difference_deg": "traj_body_sl_mismatch",
    "body_sl_axis_angle_deg": "traj_body_sl_axis_angle",
    "magnetization": "traj_magnetization",
    "abs_muB": "traj_muB",
    "body_order": "traj_body_order",
    "sl_pca_order": "traj_sl_pca_order",
    "min_gap_nm": "traj_gap_min_nm",
    "rg_nm": "traj_rg_nm",
}


def require_arviz():
    try:
        import arviz as az
    except ImportError as exc:
        raise ImportError("Install requirements.txt for multi-chain diagnostics.") from exc
    return az


def diagnose_chains(draws, elapsed_sec=None):
    """Each variable has shape (chain, post-warm-up draw). Never concatenate chains."""
    az = require_arviz()
    rows = []
    for name, values in draws.items():
        x = np.asarray(values, dtype=float)
        if x.ndim != 2 or not all(x.shape):
            raise ValueError("Diagnostics require nonempty (chain, draw) arrays.")
        chains, samples = x.shape
        row = {
            "variable": name, "chains": chains, "draws": samples,
            "mean": float(np.mean(x)) if np.all(np.isfinite(x)) else np.nan,
            "rhat_rank": np.nan, "ess_bulk": np.nan, "ess_tail": np.nan,
            "ess_mean": np.nan, "mcse_mean": np.nan, "ess_bulk_per_sec": np.nan,
            "mcse_target": DIAGNOSTIC_MCSE_TARGETS.get(name, np.nan),
            "ess_target": max(DIAGNOSTIC_ESS_MIN, DIAGNOSTIC_ESS_PER_CHAIN * chains),
        }
        if not np.all(np.isfinite(x)):
            row["status"] = "nonfinite_data"
        elif chains < 2:
            row["status"] = "insufficient_chains"
        elif samples < max(4, DIAGNOSTIC_MIN_DRAWS):
            row["status"] = "insufficient_draws"
        elif np.any(np.ptp(x, axis=1) == 0):
            # A frozen chain is not evidence of perfect precision.
            row["status"] = "constant_chain"
        else:
            row.update(
                rhat_rank=float(az.rhat(x, method="rank")),
                ess_bulk=float(az.ess(x, method="bulk")),
                ess_tail=float(az.ess(x, method="tail")),
                ess_mean=float(az.ess(x, method="mean")),
                mcse_mean=float(az.mcse(x, method="mean")),
            )
            if elapsed_sec is not None and elapsed_sec > 0:
                row["ess_bulk_per_sec"] = row["ess_bulk"] / elapsed_sec
            checks = [row[k] for k in
                      ["rhat_rank", "ess_bulk", "ess_tail", "ess_mean", "mcse_mean"]]
            reasons = []
            if not np.all(np.isfinite(checks)):
                reasons.append("undefined_diagnostic")
            else:
                if chains < 4:
                    reasons.append("fewer_than_4_chains")
                if row["rhat_rank"] >= DIAGNOSTIC_RHAT_MAX:
                    reasons.append("high_rhat")
                if min(row["ess_bulk"], row["ess_tail"], row["ess_mean"]) < row["ess_target"]:
                    reasons.append("low_ess")
                if np.isfinite(row["mcse_target"]) and row["mcse_mean"] > row["mcse_target"]:
                    reasons.append("high_mcse")
            row["status"] = ",".join(reasons) if reasons else "checks_passed"
        rows.append(row)
    return rows


def sampling_options():
    return dict(
        n_mag_per_cycle=N_MAG_PER_CYCLE,
        trans_step_nm=TRANS_STEP_NM, rot_step_deg=ROT_STEP_DEG,
        dip_step_rad=DIP_STEP_RAD, anis_model=ANIS_MODEL, include_vdw=True,
        move_positions=MOVE_POSITIONS, move_orientations=MOVE_ORIENTATIONS,
        n_cotilt_per_cycle=GLOBAL_COTILT_MOVES_PER_CYCLE,
        cotilt_step_deg=GLOBAL_COTILT_STEP_DEG,
        n_gamma_per_cycle=GLOBAL_GAMMA_MOVES_PER_CYCLE,
        gamma_step_deg=GLOBAL_GAMMA_STEP_DEG,
        show_progress=SHOW_PROGRESS_BARS, progress_backend=PROGRESS_BACKEND,
        progress_every=PROGRESS_EVERY_CYCLES,
    )


def run_chain_ensemble(starts, seeds, n_cycles, n_equil, n_mag_per_cycle=None):
    """Cross initial structures with seeds, all using exactly the same Hamiltonian."""
    az = require_arviz()  # fail before an expensive run if diagnostics are unavailable
    starts, seeds = list(starts), list(seeds)
    if not starts or not seeds or len(set(starts)) != len(starts) or len(set(seeds)) != len(seeds):
        raise ValueError("Provide nonempty, unique initial structures and seeds.")
    if any(s not in INITIAL_STRUCTURES for s in starts):
        raise ValueError("Unknown initial structure.")
    if any(not isinstance(s, (int, np.integer)) or s < 0 for s in seeds):
        raise ValueError("Seeds must be nonnegative integers.")
    if len(starts) > 1 and not (MOVE_POSITIONS and MOVE_ORIENTATIONS):
        raise ValueError("Different structures require both mechanical switches for joint diagnostics.")
    options = sampling_options()
    if n_mag_per_cycle is not None:
        options["n_mag_per_cycle"] = n_mag_per_cycle
    options.update(n_cycles=n_cycles, n_equil=n_equil)
    results, rows = [], []
    t_all = time.perf_counter()
    for structure in starts:
        # Stable within INITIAL_STRUCTURES, even if MULTI_STARTS is reordered.
        start_id = list(INITIAL_STRUCTURES).index(structure)
        for seed in seeds:
            init_seed = int(seed) if MULTI_RESEED_INITIAL_STATE else int(INIT_SEED)
            init_entropy = [init_seed, start_id, 0]
            sample_entropy = [int(seed), start_id, 1]
            state = make_initial_state(np.random.SeedSequence(init_entropy), structure)
            label = f"{structure}, seed {seed}"
            t0 = time.perf_counter()
            r = run_local_collective_mc(
                model, state, A_NM, **options, progress_label=label,
                rng=np.random.default_rng(np.random.SeedSequence(sample_entropy)),
            )
            elapsed = time.perf_counter() - t0
            results.append(r)
            rows.append({
                "chain": len(rows), "start": structure, "seed": int(seed),
                "init_entropy": init_entropy, "sample_entropy": sample_entropy,
                "elapsed_sec": elapsed, "E_mean": r["E_mean"],
                "local_beta": r["beta_mean"], "body_tilt": r["body_tilt_mean"],
                "SL_PCA": r["sl_pca_tilt_mean"],
                "body_SL_difference": r["body_sl_mismatch_mean"],
                "magnetization": r["magnetization_mean"],
                "acc_mech": r["acc_mech"], "acc_cotilt": r["acc_cotilt"],
                "acc_gamma": r["acc_gamma"], "acc_dip": r["acc_dip"],
            })
    elapsed_sec = time.perf_counter() - t_all
    draws = {name: np.stack([r[key][n_equil:] for r in results])
             for name, key in DIAGNOSTIC_TRAJECTORIES.items()}
    diagnostics = diagnose_chains(draws, elapsed_sec)
    # Report both pooled and within-start diagnostics, not averages of chain summaries.
    grouped = {}
    for structure in starts:
        indices = [i for i, row in enumerate(rows) if row["start"] == structure]
        group_sec = sum(rows[i]["elapsed_sec"] for i in indices)
        grouped[structure] = diagnose_chains(
            {name: x[indices] for name, x in draws.items()}, group_sec)
    metadata = {
        "schema_version": 2, "arviz_version": az.__version__,
        "numpy_version": np.__version__, "sampler": "Metropolis-Hastings",
        "proposal_log_ratio": 0.0, "magnetic_update": "random permutation per sweep",
        "n_cycles": n_cycles, "n_equil": n_equil,
        "elapsed_sec": elapsed_sec, "sampling": options,
        "starts": {name: dict(INITIAL_STRUCTURES[name]) for name in starts},
        "seeds": [int(s) for s in seeds], "reseed_initial_state": MULTI_RESEED_INITIAL_STATE,
        "physical_parameters": dict(model.phys), "temperature_K": T_K,
        "L_nm": L_NM, "A_nm": A_NM, "alpha_deg": ALPHA_DEG,
        "phi_init_deg": PHI_INIT_DEG, "dipole_init_mode": DIPOLE_INIT_MODE,
        "diagnostic_thresholds": {
            "min_draws": DIAGNOSTIC_MIN_DRAWS, "rhat_max": DIAGNOSTIC_RHAT_MAX,
            "ess_min": DIAGNOSTIC_ESS_MIN, "mcse_targets": dict(DIAGNOSTIC_MCSE_TARGETS),
            "ess_per_chain": DIAGNOSTIC_ESS_PER_CHAIN,
        },
        "interpretation": "Diagnostics are conditional on the existing finite-cluster model.",
    }
    return dict(results=results, rows=rows, draws=draws, diagnostics=diagnostics,
                grouped_diagnostics=grouped, metadata=metadata)


def write_csv_rows(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def export_chain_ensemble(bundle, prefix, report_pdf=None, show_figures=False,
                          save_data=True):
    """Export saved results using their recorded settings, without rerunning MC."""
    rows, results = bundle["rows"], bundle["results"]
    meta = bundle["metadata"]
    diagnostics = bundle["diagnostics"]
    if save_data:
        Path(prefix).parent.mkdir(parents=True, exist_ok=True)
        write_csv_rows(str(prefix) + "_chains.csv", rows)
        all_diags = [dict(scope="all_starts", **r) for r in diagnostics]
        for group, group_rows in bundle["grouped_diagnostics"].items():
            all_diags.extend(dict(scope=group, **r) for r in group_rows)
        write_csv_rows(str(prefix) + "_diagnostics.csv", all_diags)
        Path(str(prefix) + "_metadata.json").write_text(
            json.dumps(meta, indent=2, allow_nan=False), encoding="utf-8")
        arrays = {key: np.stack([r[key] for r in results])
                  for key in DIAGNOSTIC_TRAJECTORIES.values()}
        arrays.update({f"final_{key}": np.stack([r[key] for r in results])
                       for key in ["idx", "pos", "Q", "mu"]})
        np.savez_compressed(str(prefix) + "_trajectories.npz", **arrays)
    figures = []
    if report_pdf is not None or show_figures:
        for key, ylabel, title in [
            ("traj_body_tilt", "body tilt (deg)", "Coherent cube-body tilt"),
            ("traj_sl_pca_tilt", "SL PCA tilt (deg)", "Superlattice PCA tilt"),
        ]:
            fig, ax = plt.subplots(figsize=(9.0, 5.0))
            for row, r in zip(rows, results):
                ax.plot(np.arange(1, len(r[key]) + 1), r[key], lw=0.7,
                        label=f"{row['start']}, seed {row['seed']}")
            ax.axvline(meta["n_equil"], color="crimson", ls="--", lw=1)
            ax.set(xlabel="MC cycle", ylabel=ylabel, title=title)
            ax.grid(alpha=0.3)
            ax.legend(frameon=False, fontsize=6, ncol=2)
            fig.tight_layout()
            figures.append(fig)
    if report_pdf is not None:
        with PdfPages(report_pdf) as pdf:
            lines = [
                "Independent chains = initial structures x seeds",
                f"Structures: {list(meta['starts'])}", f"Seeds: {meta['seeds']}",
                f"Cycles / warm-up: {meta['n_cycles']} / {meta['n_equil']}",
                f"Magnetic sweeps per cycle: {meta['sampling']['n_mag_per_cycle']}",
                f"Elapsed seconds: {meta['elapsed_sec']:.3f}",
                "Sampler: Metropolis-Hastings with symmetric proposals",
                "", "Diagnostics use post-warm-up (chain, draw) arrays.",
                "Rank R-hat includes split chains and folded rank normalization.",
                "ESS measures independent-sample information, separately for bulk/tails/mean.",
                "MCSE(mean) is sampling uncertainty, in the variable's original units.",
                "Failed or undefined diagnostics do not support equilibrium estimates.",
                "Passing checks is not proof of convergence or physical model validity.",
                "Per-start diagnostics and full settings are available in CSV/JSON.",
            ]
            save_text_page(pdf, "V9.03 Multi-Chain Report", lines, fontsize=9,
                           figsize=(11, 8.5))
            for start in range(0, len(rows), 18):
                save_text_page(pdf, "Per-chain post-warm-up averages",
                               rows_to_text_table(rows[start:start + 18],
                                   ["chain", "start", "seed", "E_mean", "body_tilt", "SL_PCA"]),
                               fontsize=8, figsize=(11, 8.5))
            diagnostic_notes = [
                f"Chains: {len(rows)}. Retained draws per chain: "
                f"{meta['n_cycles'] - meta['n_equil']}.",
                f"Minimum retained draws for estimation: "
                f"{meta['diagnostic_thresholds']['min_draws']}.",
                "N/A means not estimable. It does not indicate a converged or zero-error result.",
                "Diagnostic flags: " + ", ".join(sorted({r["status"] for r in diagnostics})),
                "",
            ]
            save_text_page(pdf, "Pooled rank R-hat, ESS and MCSE",
                           diagnostic_notes + rows_to_text_table(diagnostics,
                               ["variable", "rhat_rank", "ess_bulk", "ess_tail", "ess_mean", "mcse_mean"]),
                           fontsize=7, figsize=(11, 8.5))
            save_text_page(pdf, "Diagnostic flags",
                           rows_to_text_table(diagnostics, ["variable", "status"]),
                           fontsize=8, figsize=(11, 8.5))
            for fig in figures:
                pdf.savefig(fig)
    if show_figures:
        plt.show()
    for fig in figures:
        plt.close(fig)


def compare_magnetic_sweeps(starts, seeds, candidates, n_cycles, n_equil):
    """Pilot only. Candidate chains are never pooled across sweep counts."""
    candidates = list(candidates)
    if not candidates or len(set(candidates)) != len(candidates) or any(
            not isinstance(n, (int, np.integer)) or n < 1 for n in candidates):
        raise ValueError("Sweep candidates must be distinct positive integers.")
    rows, bundles = [], {}
    for sweeps in candidates:
        bundle = run_chain_ensemble(starts, seeds, n_cycles, n_equil, sweeps)
        bundles[sweeps] = bundle
        diags = bundle["diagnostics"]
        # Require every monitored variable to pass before ranking by ESS/second.
        eligible = all(r["status"] == "checks_passed" for r in diags)
        rates = [r["ess_bulk_per_sec"] for r in diags]
        score = min(rates) if np.all(np.isfinite(rates)) else np.nan
        rows.append({"sweeps": sweeps, "elapsed_sec": bundle["metadata"]["elapsed_sec"],
                     "min_bulk_ess_per_sec": score, "eligible": eligible})
    eligible_rows = [r for r in rows if r["eligible"] and np.isfinite(r["min_bulk_ess_per_sec"])]
    recommendation = (max(eligible_rows, key=lambda r: r["min_bulk_ess_per_sec"])["sweeps"]
                      if eligible_rows else None)
    # Never change N_MAG_PER_CYCLE automatically or mix pilot and production draws.
    return dict(rows=rows, recommendation=recommendation, bundles=bundles)
