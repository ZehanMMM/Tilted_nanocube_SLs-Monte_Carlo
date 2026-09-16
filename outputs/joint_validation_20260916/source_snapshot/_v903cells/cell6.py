# Optional cross of initial structures and independent seeds.
if not RUN_MULTI_SEED_SCAN:
    print("Multi-chain scan skipped. Set RUN_MULTI_SEED_SCAN=True in Cell 1 to run.")
else:
    multi_bundle = run_chain_ensemble(MULTI_STARTS, MULTI_SEEDS, MULTI_CYCLES, MULTI_EQUIL)
    multi_results = multi_bundle["results"]
    multi_rows = multi_bundle["rows"]
    multi_diagnostics = multi_bundle["diagnostics"]
    print("Post-warm-up diagnostics. Passing checks is not proof of convergence.")
    for line in rows_to_text_table(multi_diagnostics,
            ["variable", "rhat_rank", "ess_bulk", "ess_tail", "ess_mean", "mcse_mean", "status"]):
        print(line)
    export_chain_ensemble(
        multi_bundle, MULTI_DATA_PREFIX,
        report_pdf=MULTI_SEED_REPORT_PDF if SAVE_REPORT_PDFS else None,
        show_figures=SHOW_FIGURES, save_data=SAVE_CHAIN_DATA,
    )
