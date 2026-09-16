# Sampling optimization protocol

This experiment keeps the 16 nm, 27-cube Hamiltonian and physical parameters unchanged.
The previous 2000-cycle experiment remains in `../sampling_2000cycles_20260915/`.
The pilot search was expanded sequentially after poor mixing at the smaller
proposal angles. This is a tuning study, followed by independent validation.

## Magnetic pilot

- Frozen geometries: compact/expanded spacing crossed with 0/40-degree tilt.
- Proposal maximum angles: 0.30, 0.80 and 1.60 rad.
- Complete magnetic sweeps per cycle: 1, 5 and 10.
- Seeds 11-14, with random, field, opposite-field and body-easy-axis magnetic starts.
- Each chain: 2000 cycles, first 500 discarded, no thinning.
- Four processes, one BLAS thread per process. Sampling seconds include warm-up.
- Initial and sampling random streams are separate. Initial dipoles are shared
  across candidate settings for each geometry/seed label.
- Diagnose each geometry separately. Never pool different frozen geometries.
- Monitor energy, mean local beta, signed/absolute magnetic alignment and all
  27 individual field projections and 27 individual local beta angles.
- Each candidate has 232 monitored geometry/variable combinations.
- Eligibility: rank split/folded R-hat < 1.01, bulk/tail/mean ESS >= 400,
  mean local-beta MCSE <= 0.25 degrees and signed-magnetization MCSE <= 0.002.
- Among eligible settings, rank by the minimum bulk ESS per sampling second.
  This selects only among the tested candidates, not a universal optimum.

## Independent validation

After selection, freeze the chosen magnetic settings.

- Conditional validation: the same four geometries, seeds 21-24, 4000 cycles,
  1000 warm-up. Recheck all magnetic variables, including individual particles.
- Joint validation: the four initial geometries, seeds 31-34, 4000 cycles,
  1000 warm-up. Sample centers, bodies and dipoles.
- Joint proposals include one uniform center-scaling move per cycle with
  maximum log-scale 0.002, using the `3*(N-1)*epsilon` volume Jacobian, and
  one rigid body/center co-tilt move with maximum angle 3 degrees.
- The body/SL MCSE target remains 0.5 degrees. ESS threshold is the greater
  of 400 and 100 times the chain count. Report pooled and within-start checks.
- Keep every completed run even if diagnostics fail. Do not tune on validation
  draws or alter the warm-up slice to make checks pass.

Conditional magnetic targets are normalizable. The full position space remains
unbounded, so joint diagnostics characterize finite-run mixing and cannot prove
unrestricted assembly equilibrium. The scale Jacobian preserves the original
Cartesian sampling measure and is not an added configurational-entropy energy.

## Provenance and checks

Per-chain checkpoints contain all trajectories, final state, final RNG state,
configuration, elapsed time, acceptance rates and numerical checks. Each run
stores exact source snapshots and SHA256 hashes. Rebuild the notebook and run
the regression suite before releasing changes. Preserve historical results.
