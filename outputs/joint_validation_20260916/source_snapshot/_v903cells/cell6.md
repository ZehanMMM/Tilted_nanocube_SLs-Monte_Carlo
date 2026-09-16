## Sampling and diagnostics

Set all controls in the first code cell. The interaction potentials and their
physical parameters are unchanged in this sampling update.

- `N_MAG_PER_CYCLE = 5` means five complete random-order magnetic sweeps.
- `RUN_MAGNETIC_PILOT = True` compares `MAGNETIC_SWEEP_CANDIDATES` using a
  separate pilot. Only candidates passing all diagnostic checks are ranked by
  the smallest bulk ESS per second across monitored observables. A short pilot
  may return no recommendation. Change the production count explicitly afterward.
- `RUN_MULTI_SEED_SCAN = True` runs every combination of `MULTI_STARTS` and
  `MULTI_SEEDS`. Defaults produce 12 independent chains. This is opt-in work.
- Initial conditions include aligned compact, co-tilted compact and expanded
  co-tilted clusters. Expansion changes initial positions only, not `A_NM` or
  the vdW cutoff. All chains use one Hamiltonian.
- Each chain uses separate initialization and sampling streams derived from
  `SeedSequence([seed, structure_id, stream_id])`. The actual sequences are saved.

Rank-normalized split R-hat compares between-chain and within-chain variation
using ranked and folded draws. Values below 1.01 are a diagnostic target, not
proof of equilibrium. ESS estimates independent-sample information. Bulk ESS
covers the distribution center, tail ESS covers the 5% and 95% tails, and mean
ESS is used for mean precision. MCSE(mean) estimates sampling uncertainty of the
mean in the observable's original units. It is not the physical spread of the
observable or uncertainty of the interaction model.

Diagnostics use unthinned post-warm-up arrays with shape `(chain, draw)`.
They are computed both across all starts and within each start. Fewer than
100 retained draws, constant chains, nonfinite values, high R-hat and low ESS
are flagged. At least four chains are required to pass the diagnostic gate.
The default ESS target is 400 and tilt MCSE targets are 0.5 degrees.
No automatic extension or stopping rule is implied by these thresholds.

Body tilt and SL PCA tilt remain on separate multi-chain figures. Raw trajectories,
final states, chain identities, run settings and diagnostic results are exported
as NPZ, JSON and CSV. Use the final cell to re-export the saved `multi_bundle`
without MC. Its recorded settings are retained even if current controls change.
The 0-based diagnostic slice `N_EQUIL:` removes exactly `N_EQUIL` cycles.

The single-run and multi-chain PDFs are controlled by `SAVE_REPORT_PDFS`.
The separate Fig1 PDF follows the same switch. `SAVE_CHAIN_DATA` controls raw
multi-chain data. `SHOW_FIGURES` controls display only. Existing published PDFs
and the SI are historical results and have not been recomputed.

The current model has no explicit finite container. Diagnostics assess sampling
under the existing cluster assumptions and do not establish unrestricted
assembly equilibrium. Entropy and interaction models will be reviewed separately.
