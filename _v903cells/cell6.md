## Notes

The default run starts from `16.0: {'d': 21.0, 'alpha': 74.2}` and
`beta_init = 0 deg`. To switch among experimental sizes, edit
`DEFAULT_SIZE_NM`. To test tilted candidates suggested by V9.02, set
`RUN_PRESET` to `19p5_tilted_candidate` or `29p6_tilted_candidate`.

For a direct straight-vs-tilted comparison, duplicate Cell 5 with
`DEFAULT_BETA_INIT_DEG = 0` / `20` or switch `RUN_PRESET`, then compare
equilibrated `<E>`, local cube-body `<beta>`, coherent cube-body tilt/order,
SL PCA tilt, body-SL mismatch, and `<|mu.B|>`.

`GLOBAL_GAMMA_MOVES_PER_CYCLE` controls the V6.20-style coherent cube twist:
all cube bodies rotate about the coherent cube `[111]` axis while the cube
centers remain fixed. This samples face-alignment/gamma relaxation without
creating local cube-SL tilt mismatch.

If the 2000-cycle trajectory is still drifting, increase to `10000` cycles and
use the final 50-70% of the trajectory for averages.

All random seeds are set in the first code cell. The default is
`DEFAULT_SEED = 1`. To run several trajectories, set
`RUN_MULTI_SEED_SCAN = True` and edit `MULTI_SEEDS`, for example
`[1, 2, 3]`.

PDF reports are saved by default. The single full run writes
`V903_FullRun_Report.pdf`; the optional multi-seed scan writes
`V903_MultiSeed_Report.pdf`. Set `SAVE_REPORT_PDFS = False` to disable this.
Set `SHOW_FIGURES = False` to save reports without displaying figures.

If a long multi-seed scan has already finished and the kernel still has
`multi_rows` / `multi_results` in memory, run the final export-only cell to
write `V903_MultiSeed_Report.pdf` without rerunning MC.
