# Optional multi-seed scan. Set RUN_MULTI_SEED_SCAN=True in Cell 1 to run it.
if not RUN_MULTI_SEED_SCAN:
    print("Multi-seed scan skipped. Set RUN_MULTI_SEED_SCAN=True in Cell 1 to run.")
else:
    multi_results = []
    multi_rows = []

    for seed in MULTI_SEEDS:
        print(f"Running seed {seed} for {MULTI_CYCLES} cycles...")
        state_seed = make_initial_state(seed) if MULTI_RESEED_INITIAL_STATE else state0

        t0 = time.perf_counter()
        r = run_local_collective_mc(
            model, state_seed, A_NM,
            n_cycles=MULTI_CYCLES, n_equil=MULTI_EQUIL,
            n_mag_per_cycle=N_MAG_PER_CYCLE,
            trans_step_nm=TRANS_STEP_NM,
            rot_step_deg=ROT_STEP_DEG,
            dip_step_rad=DIP_STEP_RAD,
            anis_model=ANIS_MODEL,
            include_vdw=True,
            move_positions=MOVE_POSITIONS,
            move_orientations=MOVE_ORIENTATIONS,
            n_cotilt_per_cycle=GLOBAL_COTILT_MOVES_PER_CYCLE,
            cotilt_step_deg=GLOBAL_COTILT_STEP_DEG,
            n_gamma_per_cycle=GLOBAL_GAMMA_MOVES_PER_CYCLE,
            gamma_step_deg=GLOBAL_GAMMA_STEP_DEG,
            show_progress=SHOW_PROGRESS_BARS,
            progress_label=f"multi seed {seed}",
            progress_backend=PROGRESS_BACKEND,
            progress_every=PROGRESS_EVERY_CYCLES,
            rng=np.random.default_rng(seed),
        )
        elapsed = time.perf_counter() - t0

        multi_results.append(r)
        multi_rows.append({
            "seed": seed,
            "elapsed_min": elapsed / 60.0,
            "E_mean": r["E_mean"],
            "local_beta": r["beta_mean"],
            "body_tilt": r["body_tilt_mean"],
            "SL_PCA": r["sl_pca_tilt_mean"],
            "body_SL_mismatch": r["body_sl_mismatch_mean"],
            "body_order": r["body_order_mean"],
            "gap_min": r["gap_min_mean"],
            "gap_p05": r["gap_p05_mean"],
            "muB": r["muB_mean"],
            "acc_mech": r["acc_mech"],
            "acc_cotilt": r["acc_cotilt"],
            "acc_gamma": r["acc_gamma"],
            "acc_dip": r["acc_dip"],
        })

    multi_columns = [
        "seed", "elapsed_min", "E_mean", "local_beta", "body_tilt", "SL_PCA",
        "body_SL_mismatch", "body_order", "gap_min", "gap_p05",
        "muB", "acc_mech", "acc_cotilt", "acc_gamma", "acc_dip",
    ]
    multi_metric_columns = [c for c in multi_columns if c != "seed"]
    multi_agg_rows = []
    for stat, func in [("mean", np.mean), ("std", np.std)]:
        row = {"stat": stat}
        for col in multi_metric_columns:
            vals = np.array([r[col] for r in multi_rows], dtype=float)
            row[col] = func(vals, ddof=1) if stat == "std" and len(vals) > 1 else func(vals)
        multi_agg_rows.append(row)

    try:
        import pandas as pd
        df_multi = pd.DataFrame(multi_rows)
        df_multi_summary = df_multi.drop(columns=["seed"]).agg(["mean", "std"])
        try:
            display(df_multi)
            display(df_multi_summary)
        except NameError:
            print(df_multi.to_string(index=False))
            print(df_multi_summary.to_string())
    except ImportError:
        df_multi = multi_rows
        for row in multi_rows:
            print(row)

    fig_body, ax_body = plt.subplots(figsize=(7.0, 4.2))
    for row, r in zip(multi_rows, multi_results):
        ax_body.plot(r["traj_body_tilt"], lw=0.8, label=f"seed {row['seed']}")
    ax_body.axvline(MULTI_EQUIL, color="crimson", ls="--", lw=1)
    ax_body.set(xlabel="MC cycle", ylabel="body tilt vs lab-x (deg)",
                title="multi-seed coherent cube-body tilt")
    ax_body.grid(alpha=0.3)
    ax_body.legend(frameon=False, fontsize=8, ncol=2)
    plt.tight_layout()

    fig_sl, ax_sl = plt.subplots(figsize=(7.0, 4.2))
    for row, r in zip(multi_rows, multi_results):
        ax_sl.plot(r["traj_sl_pca_tilt"], lw=0.8, label=f"seed {row['seed']}")
    ax_sl.axvline(MULTI_EQUIL, color="crimson", ls="--", lw=1)
    ax_sl.set(xlabel="MC cycle", ylabel="SL PCA tilt vs lab-x (deg)",
              title="multi-seed superlattice PCA tilt")
    ax_sl.grid(alpha=0.3)
    ax_sl.legend(frameon=False, fontsize=8, ncol=2)
    plt.tight_layout()

    if SAVE_REPORT_PDFS:
        with PdfPages(MULTI_SEED_REPORT_PDF) as pdf:
            config_lines = [
                "Multi-seed configuration",
                "",
                *base_report_lines(),
                "",
                f"MULTI_SEEDS                = {MULTI_SEEDS}",
                f"MULTI_CYCLES, MULTI_EQUIL  = {MULTI_CYCLES}, {MULTI_EQUIL}",
                f"MULTI_RESEED_INITIAL_STATE = {MULTI_RESEED_INITIAL_STATE}",
                "",
                "Each per-seed row is the equilibrated mean over cycles "
                "MULTI_EQUIL..MULTI_CYCLES-1.",
            ]
            save_text_page(pdf, "V9.03 Multi-Seed Configuration", config_lines)

            per_seed_groups = [
                ("Tilt and energy",
                 ["seed", "elapsed_min", "E_mean", "local_beta", "body_tilt",
                  "SL_PCA", "body_SL_mismatch"]),
                ("Order and gaps",
                 ["seed", "body_order", "gap_min", "gap_p05", "muB"]),
                ("Acceptance",
                 ["seed", "acc_mech", "acc_cotilt", "acc_gamma", "acc_dip"]),
            ]
            for label, cols in per_seed_groups:
                save_text_page(
                    pdf,
                    f"V9.03 Multi-Seed Per-Seed Results: {label}",
                    rows_to_text_table(multi_rows, cols),
                    fontsize=8,
                    figsize=(11.0, 8.5),
                )

            agg_groups = [
                ("Tilt and energy",
                 ["stat", "elapsed_min", "E_mean", "local_beta", "body_tilt",
                  "SL_PCA", "body_SL_mismatch"]),
                ("Order and gaps",
                 ["stat", "body_order", "gap_min", "gap_p05", "muB"]),
                ("Acceptance",
                 ["stat", "acc_mech", "acc_cotilt", "acc_gamma", "acc_dip"]),
            ]
            for label, cols in agg_groups:
                save_text_page(
                    pdf,
                    f"V9.03 Multi-Seed Mean and Std: {label}",
                    rows_to_text_table(multi_agg_rows, cols),
                    fontsize=8,
                    figsize=(11.0, 8.5),
                )

            pdf.savefig(fig_body, bbox_inches="tight")
            pdf.savefig(fig_sl, bbox_inches="tight")
        print(f"Wrote PDF report: {MULTI_SEED_REPORT_PDF}")

    if SHOW_FIGURES:
        plt.show()
    else:
        plt.close(fig_body)
        plt.close(fig_sl)
