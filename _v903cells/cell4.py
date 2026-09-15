# Optional pilot. Results are separate from production chains.
if RUN_MAGNETIC_PILOT:
    magnetic_pilot = compare_magnetic_sweeps(
        MULTI_STARTS, MULTI_SEEDS, MAGNETIC_SWEEP_CANDIDATES, PILOT_CYCLES, PILOT_EQUIL)
    for line in rows_to_text_table(magnetic_pilot["rows"],
            ["sweeps", "elapsed_sec", "min_bulk_ess_per_sec", "eligible"]):
        print(line)
    print("Recommended sweeps:", magnetic_pilot["recommendation"])
    if magnetic_pilot["recommendation"] is None:
        print("No eligible candidate. Improve mixing or extend the pilot before selecting.")
    if SAVE_CHAIN_DATA:
        write_csv_rows(MULTI_DATA_PREFIX + "_pilot.csv", magnetic_pilot["rows"])
        for sweeps, bundle in magnetic_pilot["bundles"].items():
            export_chain_ensemble(bundle, MULTI_DATA_PREFIX + f"_pilot_{sweeps}")

# Short benchmark: estimate runtime before the full 2000-cycle run.

t0 = time.perf_counter()
bench = run_local_collective_mc(model, state0, A_NM,
                                n_cycles=BENCH_CYCLES, n_equil=BENCH_EQUIL,
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
                                show_progress=SHOW_BENCH_PROGRESS,
                                progress_label=f"bench seed {BENCH_SEED}",
                                progress_backend=PROGRESS_BACKEND,
                                progress_every=PROGRESS_EVERY_CYCLES,
                                rng=np.random.default_rng(BENCH_SEED))
bench_sec = time.perf_counter() - t0
sec_per_cycle = bench_sec / BENCH_CYCLES
est_full_min = sec_per_cycle * N_CYCLES / 60.0

print("Benchmark")
print(f"  seed = {BENCH_SEED}")
print(f"  {BENCH_CYCLES} cycles took {bench_sec:.2f} s")
print(f"  {sec_per_cycle:.3f} s/cycle")
print(f"  estimated {N_CYCLES} cycles: {est_full_min:.1f} min")
print(f"  bench <E>={bench['E_mean']:+.2f} kBT, local <beta>={bench['beta_mean']:.2f} deg, "
      f"<|mu.B|>={bench['muB_mean']:.3f}")
print(f"  bench coherent body tilt={bench['body_tilt_mean']:.2f} deg "
      f"(order={bench['body_order_mean']:.3f}), "
      f"body-SL tilt difference={bench['body_sl_mismatch_mean']:.2f} deg")
print(f"  bench SL tilt PCA={bench['sl_pca_tilt_mean']:.2f} deg "
      f"(order={bench['sl_pca_order_mean']:.3f})")
print(f"  bench gap min/p05={bench['gap_min_mean']:.3f} / {bench['gap_p05_mean']:.3f} nm")
print(f"  acceptance mech/cotilt/gamma/dip = {bench['acc_mech']:.3f} / "
      f"{bench['acc_cotilt']:.3f} / {bench['acc_gamma']:.3f} / {bench['acc_dip']:.3f}")
