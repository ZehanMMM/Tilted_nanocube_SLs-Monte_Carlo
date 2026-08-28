# Full V9.03 N27 run.
t0 = time.perf_counter()
res = run_local_collective_mc(model, state0, A_NM,
                              n_cycles=N_CYCLES, n_equil=N_EQUIL,
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
                              progress_label=f"full seed {FULL_RUN_SEED}",
                              progress_backend=PROGRESS_BACKEND,
                              progress_every=PROGRESS_EVERY_CYCLES,
                              rng=np.random.default_rng(FULL_RUN_SEED))
elapsed = time.perf_counter() - t0

print("V9.03 N27 run summary")
print(f"  init/full seed = {INIT_SEED} / {FULL_RUN_SEED}")
print(f"  elapsed = {elapsed:.1f} s = {elapsed/60:.2f} min")
print(f"  <E>     = {res['E_mean']:+.3f} kBT")
print(f"  local cube-body <beta> = {res['beta_mean']:.3f} deg")
print(f"  coherent cube-body tilt = {res['body_tilt_mean']:.3f} deg "
      f"(body order={res['body_order_mean']:.3f})")
print(f"  coherent body-SL mismatch = {res['body_sl_mismatch_mean']:.3f} deg")
print(f"  SL tilt PCA      = {res['sl_pca_tilt_mean']:.3f} deg "
      f"(order={res['sl_pca_order_mean']:.3f})")
print(f"  surface gap min/p05 = {res['gap_min_mean']:.3f} / {res['gap_p05_mean']:.3f} nm")
print(f"  <mu.B>  = {res['muB_mean']:.3f}")
print(f"  acceptance mech/cotilt/gamma/dip = {res['acc_mech']:.3f} / "
      f"{res['acc_cotilt']:.3f} / {res['acc_gamma']:.3f} / {res['acc_dip']:.3f}")

full_report_rows = [
    {
        "stat": "equil_mean",
        "E_kBT": res["E_mean"],
        "local_beta": res["beta_mean"],
        "body_tilt": res["body_tilt_mean"],
        "SL_PCA": res["sl_pca_tilt_mean"],
        "body_SL_mismatch": res["body_sl_mismatch_mean"],
        "body_order": res["body_order_mean"],
        "gap_min": res["gap_min_mean"],
        "gap_p05": res["gap_p05_mean"],
        "muB": res["muB_mean"],
    },
    {
        "stat": "final",
        "E_kBT": res["traj_E"][-1],
        "local_beta": res["traj_beta"][-1],
        "body_tilt": res["traj_body_tilt"][-1],
        "SL_PCA": res["traj_sl_pca_tilt"][-1],
        "body_SL_mismatch": res["traj_body_sl_mismatch"][-1],
        "body_order": res["traj_body_order"][-1],
        "gap_min": res["traj_gap_min_nm"][-1],
        "gap_p05": res["traj_gap_p05_nm"][-1],
        "muB": res["traj_muB"][-1],
    },
]

fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.6))

ax = axes[0, 0]
ax.plot(res["traj_E"], lw=0.7, color="black")
ax.axvline(N_EQUIL, color="crimson", ls="--", lw=1)
ax.set(xlabel="MC cycle", ylabel="E (kBT)", title="N27 energy trajectory")
ax.grid(alpha=0.3)

ax = axes[0, 1]
ax.plot(res["traj_body_tilt"], lw=0.9, color="teal", label="coherent cube-body [111]")
ax.plot(res["traj_sl_pca_tilt"], lw=0.9, color="darkorange", label="SL PCA")
ax.plot(res["traj_beta"], lw=0.6, color="teal", alpha=0.35, ls=":", label="mean local beta")
ax.axvline(N_EQUIL, color="crimson", ls="--", lw=1)
ax.set(xlabel="MC cycle", ylabel="tilt vs lab-x (deg)", title="coherent body vs SL tilt")
ax.legend(frameon=False, fontsize=8)
ax.grid(alpha=0.3)

ax = axes[1, 0]
ax.plot(res["traj_muB"], lw=0.7, color="purple")
ax.axvline(N_EQUIL, color="crimson", ls="--", lw=1)
ax.set(xlabel="MC cycle", ylabel="<|mu.Bhat|>", title="dipole unlocking")
ax.set_ylim(0.0, 1.02)
ax.grid(alpha=0.3)

ax = axes[1, 1]
ax.plot(res["traj_body_order"], lw=0.9, color="teal", label="body vector order")
ax.plot(res["traj_sl_pca_order"], lw=0.9, color="darkorange", label="SL PCA order")
ax.axvline(N_EQUIL, color="crimson", ls="--", lw=1)
ax.set(xlabel="MC cycle", ylabel="axis order", title="body and SL PCA order parameters")
ax.set_ylim(0.0, 1.0)
ax.legend(frameon=False, fontsize=8)
ax.grid(alpha=0.3)

fig.suptitle(f"V9.03 N27 collective MC, B=500 G, cubic-first; {RUN_LABEL}", y=1.02)
plt.tight_layout()
plt.savefig(FULL_RUN_FIG_PDF, bbox_inches="tight", dpi=150)

if SAVE_REPORT_PDFS:
    with PdfPages(FULL_RUN_REPORT_PDF) as pdf:
        lines = [
            "Full-run configuration",
            "",
            *base_report_lines(),
            "",
            f"INIT_SEED, FULL_RUN_SEED = {INIT_SEED}, {FULL_RUN_SEED}",
            f"N_CYCLES, N_EQUIL        = {N_CYCLES}, {N_EQUIL}",
            f"elapsed_min              = {elapsed / 60.0:.6g}",
            f"acc_mech/cotilt/gamma/dip = {res['acc_mech']:.6g} / "
            f"{res['acc_cotilt']:.6g} / {res['acc_gamma']:.6g} / {res['acc_dip']:.6g}",
            "",
            "Values marked equil_mean are averaged over cycles N_EQUIL..N_CYCLES-1.",
            "Values marked final are the last recorded cycle.",
        ]
        save_text_page(pdf, "V9.03 Full Run Report", lines)
        save_text_page(
            pdf,
            "V9.03 Full Run Tilt Summary",
            rows_to_text_table(
                full_report_rows,
                ["stat", "E_kBT", "local_beta", "body_tilt", "SL_PCA",
                 "body_SL_mismatch"],
            ),
            fontsize=8,
            figsize=(11.0, 8.5),
        )
        save_text_page(
            pdf,
            "V9.03 Full Run Order, Gap, and Acceptance",
            rows_to_text_table(
                full_report_rows,
                ["stat", "body_order", "gap_min", "gap_p05", "muB"],
            ),
            fontsize=8,
            figsize=(11.0, 8.5),
        )
        pdf.savefig(fig, bbox_inches="tight")
    print(f"Wrote PDF report: {FULL_RUN_REPORT_PDF}")

if SHOW_FIGURES:
    plt.show()
else:
    plt.close(fig)
