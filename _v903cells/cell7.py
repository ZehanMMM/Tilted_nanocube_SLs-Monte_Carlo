# Export the saved multi-chain bundle with its original run settings.
if "multi_bundle" not in globals():
    print("No multi-chain bundle found. Old multi_rows alone lack diagnostic provenance.")
else:
    export_chain_ensemble(
        multi_bundle, MULTI_DATA_PREFIX, report_pdf=MULTI_SEED_REPORT_PDF,
        show_figures=SHOW_FIGURES, save_data=SAVE_CHAIN_DATA,
    )
    print(f"Exported existing results to {MULTI_SEED_REPORT_PDF} without rerunning MC.")
