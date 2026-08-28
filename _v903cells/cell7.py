# Export already-computed multi-seed results without rerunning MC.
# Use this if a long multi-seed scan has already finished in the current kernel.
if "multi_rows" not in globals() or "multi_results" not in globals() or not multi_rows:
    print("No completed multi-seed results found in memory. Run the multi-seed cell first.")
else:
    export_columns = [
        "seed", "elapsed_min", "E_mean", "local_beta", "body_tilt", "SL_PCA",
        "body_SL_mismatch", "body_order", "gap_min", "gap_p05",
        "muB", "acc_mech", "acc_cotilt", "acc_gamma", "acc_dip",
    ]
    export_metric_columns = [c for c in export_columns if c != "seed"]
    export_agg_rows = []
    for stat, func in [("mean", np.mean), ("std", np.std)]:
        row = {"stat": stat}
        for col in export_metric_columns:
            vals = np.array([r.get(col, np.nan) for r in multi_rows], dtype=float)
            row[col] = func(vals, ddof=1) if stat == "std" and len(vals) > 1 else func(vals)
        export_agg_rows.append(row)

    export_fig_body, export_ax_body = plt.subplots(figsize=(7.0, 4.2))
    for row, r in zip(multi_rows, multi_results):
        export_ax_body.plot(r["traj_body_tilt"], lw=0.8, label=f"seed {row['seed']}")
    export_ax_body.axvline(MULTI_EQUIL, color="crimson", ls="--", lw=1)
    export_ax_body.set(xlabel="MC cycle", ylabel="body tilt vs lab-x (deg)",
                       title="multi-seed coherent cube-body tilt")
    export_ax_body.grid(alpha=0.3)
    export_ax_body.legend(frameon=False, fontsize=8, ncol=2)
    plt.tight_layout()

    export_fig_sl, export_ax_sl = plt.subplots(figsize=(7.0, 4.2))
    for row, r in zip(multi_rows, multi_results):
        export_ax_sl.plot(r["traj_sl_pca_tilt"], lw=0.8, label=f"seed {row['seed']}")
    export_ax_sl.axvline(MULTI_EQUIL, color="crimson", ls="--", lw=1)
    export_ax_sl.set(xlabel="MC cycle", ylabel="SL PCA tilt vs lab-x (deg)",
                     title="multi-seed superlattice PCA tilt")
    export_ax_sl.grid(alpha=0.3)
    export_ax_sl.legend(frameon=False, fontsize=8, ncol=2)
    plt.tight_layout()

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
                rows_to_text_table(export_agg_rows, cols),
                fontsize=8,
                figsize=(11.0, 8.5),
            )

        pdf.savefig(export_fig_body, bbox_inches="tight")
        pdf.savefig(export_fig_sl, bbox_inches="tight")

    print(f"Wrote PDF report from current in-memory results: {MULTI_SEED_REPORT_PDF}")
    if SHOW_FIGURES:
        plt.show()
    else:
        plt.close(export_fig_body)
        plt.close(export_fig_sl)
