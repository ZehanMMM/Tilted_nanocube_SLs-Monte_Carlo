# V9.03 Nanocube Collective Monte Carlo Simulation

This repository contains the V9.03 Monte Carlo model used to study the coupled reorientation of a compact 3 x 3 x 3 magnetite nanocube superlattice in a 500 G magnetic field.

## Model

The simulated system contains 27 cubes with a default inorganic edge length of 16 nm. The initial structure is compact and rhombohedrally tilted. Each Monte Carlo cycle attempts collective cube translation, cube-body rotation, and magnetic dipole rotation according to the switches in the first notebook cell.

The energy includes Zeeman, cubic magnetocrystalline anisotropy, dipole-dipole, van der Waals, and steric terms. The superlattice tilt is obtained from principal component analysis of the cube-center coordinates. Per-cube body tilt is measured between body [111] and the magnetic field. The reported coherent body tilt uses the mean body-[111] direction and treats its sign as equivalent. Body order records the length of that mean direction. Monte Carlo cycles are sampling steps and are not interpreted as physical time.

Default settings include:

- Magnetic field: 500 G (0.05 T)
- Temperature: 298.15 K
- Particle size: 16 nm
- Cluster: 27 cubes in a full 3 x 3 x 3 arrangement
- Cycles: 2000
- Equilibration cycles: 500
- Random seed: 1
- Effective initial surface gap: 1.8 nm
- Corner-rounding parameter: 1.5 nm
- Magnetic updates: 10 complete random-order sweeps per cycle
- Magnetic proposal maximum angle: 1.60 rad, selected by conditional pilots
- Collective center scaling: one move per cycle, maximum log-scale 0.002
- Rigid body/center co-tilt: one move per cycle, maximum angle 3 degrees
- Optional multi-chain design: 5 initial structures x 4 seeds

## Repository Structure

- `V9.03_N27_500G.ipynb`: executable notebook
- `_v903cells/`: editable notebook source cells and rebuild script
- `tests/`: short sampler, diagnostic and notebook regression tests
- `docs/MCMC_DIAGNOSTICS_ZH.md`: Chinese explanation of the experiment and diagnostics
- `outputs/`: 16 nm reports, calibration and validation checkpoints, and diagnostics
- `supporting_information/`: concise Supporting Information text in Word format

The notebook is generated from `_v903cells`. Edit the latest source cells and run `_v903cells/build.py` to rebuild it. The original historical versions are not included or modified.

## Running the Simulation

Create a Python environment and install the dependencies:

```bash
python -m pip install -r requirements.txt
```

Open the notebook:

```bash
jupyter notebook V9.03_N27_500G.ipynb
```

Simulation controls are grouped near the beginning of `_v903cells/cell1.py` and in the first notebook setup cell. Set `DEFAULT_SEED`, `N_CYCLES`, `N_EQUIL`, and the move switches there before running the notebook from top to bottom. Running all cells executes the benchmark and the full single-chain run. Multi-chain and magnetic-sweep pilot experiments remain disabled by default.

## Sampling Experiment

The sampler uses Metropolis-Hastings with acceptance probability
`min(1, exp(-delta_U/kBT + log_q_reverse - log_q_forward))`.
Current translation, body rotation, dipole rotation and global rotation proposals
are symmetric, so the log proposal ratio is zero. Acceptance is evaluated in
log space. Increasing the magnetic sweep count changes sampling efficiency,
not the interaction model.

`N_MAG_PER_CYCLE = 10` means ten complete random permutations of the
particles per cycle. Each particle receives one dipole proposal per sweep.
Magnetic proposals evaluate only Zeeman, anisotropy and dipolar terms.
The unchanged vdW and steric terms cancel for these proposals.

Set `RUN_MULTI_SEED_SCAN = True` to cross `MULTI_STARTS` with `MULTI_SEEDS`.
Defaults use seeds 1, 2, 3 and 4 with these starts:

- `compact_aligned`: the compact cluster at 0 degrees
- `compact_tilted`: the same compact geometry co-rotated by 20 degrees
- `expanded_tilted`: an 8% expansion, co-rotated by 40 degrees
- `compact_tilted40`: the compact cluster co-rotated by 40 degrees
- `expanded_aligned`: an 8% expansion at zero tilt

The third start changes interparticle distances, not just global orientation.
These are dispersed starts within the existing N27 cluster family, not a survey
of all assembly topologies. Bodies and the SL axis are initially aligned.
All chains use the same size, field, temperature, energy parameters and `A_NM`
cutoff reference. Different structures require both mechanical switches to
remain enabled for the pooled comparison.

Initialization and sampling have separate reproducible NumPy SeedSequence
streams. The structure identifier is its position in `INITIAL_STRUCTURES`.
The actual entropy tuples are recorded in the chain CSV, together with the seed
and start label. No pilot samples are included in production diagnostics.

Set `RUN_MAGNETIC_PILOT = True` to compare 1, 5 and 10 sweeps in joint
sampling using `PILOT_CYCLES` and `PILOT_EQUIL`. Fixed-geometry magnetic
calibration uses the separate script described below. A candidate must pass the diagnostic gates
for every monitored observable before it can be recommended. Eligible candidates
are ranked by their worst bulk ESS per elapsed second, including warm-up time.
If none qualifies, the recommendation is absent. The pilot never changes
`N_MAG_PER_CYCLE` automatically. The current ten-sweep default was selected
by the separate fixed-geometry calibration, not this optional joint pilot.

## Interpreting Diagnostics

ArviZ 0.22.0 computes diagnostics on unthinned post-warm-up arrays shaped
`(chain, draw)`, both pooled across starts and separately within each start.
Chain means are not used as a substitute for the underlying samples.

- **Rank-normalized split R-hat** compares variation within and across split
  chains after rank normalization, including a folded diagnostic sensitive to
  scale differences. The default target is below 1.01.
- **Bulk, tail and mean ESS** estimate how much independent-sample information
  remains after autocorrelation. Tail ESS uses the 5% and 95% quantiles.
  Mean ESS is distinct from bulk ESS. The minimum is the larger of 400 and
  100 times the number of chains.
- **MCSE of the mean** estimates Monte Carlo uncertainty in the mean, in the
  original units. The default body-tilt and SL-tilt targets are 0.5 degrees.
  MCSE is neither the physical standard deviation nor model uncertainty.

Fewer than two chains or 100 retained draws, constant chains and nonfinite
observations produce unavailable diagnostics and explicit flags. Passing the
gate also requires at least four chains. High R-hat, insufficient ESS and
excessive MCSE are flagged. `checks_passed` is not proof of convergence.
The fixed 2000/500 cycle settings do not imply equilibrium.

Multi-chain export writes `V903_MultiChain_chains.csv`,
`V903_MultiChain_diagnostics.csv`, `V903_MultiChain_metadata.json` and
`V903_MultiChain_trajectories.npz`. The NPZ retains every cycle, including
warm-up, and final states. Metadata records the warm-up slice and sampler/model
settings. Body and SL tilts remain on separate PDF figures. The final notebook
cell re-exports `multi_bundle` with its recorded settings without rerunning MC.

Methods: [Vehtari et al. (2021)](https://arxiv.org/abs/1903.08008),
[ArviZ R-hat](https://python.arviz.org/en/v0.22.0/api/generated/arviz.rhat.html),
[ESS](https://python.arviz.org/en/v0.22.0/api/generated/arviz.ess.html), and
[MCSE](https://python.arviz.org/en/v0.22.0/api/generated/arviz.mcse.html).

## Short Validation

```bash
python _v903cells/build.py
python -m unittest discover -s tests -v
```

The tests use short N27 chains and synthetic diagnostic fixtures. They do not
run the 2000-cycle experiment or the full multi-chain/pilot configurations.

## Full Experiment with Checkpoints

Run the default single chain and all 20 multi-start chains in separate processes:

```bash
python scripts/run_sampling_experiment.py --output outputs/my_full_run --cycles 2000 --equil 500 --workers 4
```

Use a new output directory to preserve previous results. Each completed chain
is saved immediately as NPZ and JSON, together with a progress log and energy
consistency checks. A manifest records the source hashes, base commit and chain
configuration. Resume an interrupted run with the same command plus `--resume`.
Resume requires identical source files and cycle settings and skips completed
chains. The final report pools only the multi-start chains, excluding the
separate single-chain reproduction.

All chains use the sampling parameters in the first source cell. This command
does not run the magnetic-sweep optimization pilot. Reported ESS/second uses the sum of per-chain
sampling times, including warm-up. Elapsed wall time is recorded separately.

The historical [2000-cycle experiment](outputs/sampling_2000cycles_20260915/RESULTS_ZH.md)
contains one default single-chain reproduction and 12 multi-start chains, each
with 500 warm-up cycles. All trajectories and estimated diagnostics are finite.
None of the 12 pooled observables passes the predefined diagnostic checks.
Body and SL tilt rank R-hat values are 1.830 and 1.807, with bulk ESS of 17.4
and 17.7. Expanded starts retain substantial drift after warm-up. These outputs
document incomplete mixing and must not be interpreted as equilibrium estimates.
That run used five sweeps, a 0.30 rad magnetic step and no collective scaling.

## Magnetic Calibration and Longer Validation

The additional compact 40-degree and expanded zero-degree starts separate
initial tilt from initial spacing. Their structure IDs are appended to preserve
the existing random streams. ESS checks now require at least the larger of
400 and 100 times the number of chains. Saved historical results retain their
original settings and thresholds.

```bash
python scripts/run_sampling_optimization.py --mode conditional --output outputs/my_magnetic_calibration --cycles 2000 --equil 500 --sweeps 1 5 10 --dip-step-rad 1.6 --workers 4
python scripts/run_sampling_optimization.py --mode joint --output outputs/my_joint_validation --cycles 4000 --equil 1000 --sweeps 10 --seeds 31 32 33 34 --workers 4
```

Conditional calibration freezes each of four geometries and uses seeds 11-14
with random, field-aligned, opposite-field and body-easy-axis magnetic starts.
It checks individual particles as well as averages. Different frozen geometries
are diagnosed separately because they define different conditional distributions.
Only candidates passing all checks are ranked by the worst bulk ESS per second.
A conditional recommendation is not a proven optimum for joint sampling.

The completed pilots tested 0.30, 0.80 and 1.60 rad, each crossed with 1, 5 and
10 sweeps. Only 1.60 rad with 10 sweeps passed all 232 geometry/variable checks
in the pilot. Its maximum rank R-hat was 1.0079 and minimum bulk ESS was 583.3.
The [experiment protocol](outputs/sampling_optimization_20260916/PROTOCOL.md)
separates tuning from independent conditional and joint validation.
The independent conditional validation used seeds 21-24 and 4000/1000
cycles. All 232 checks passed, with maximum rank R-hat 1.0030 and minimum
bulk ESS 1129.0. These checks apply to the four tested frozen geometries.

The independent joint validation used four structures, seeds 31-34 and the same
4000/1000 cycle lengths. All 60 pooled and within-start diagnostic rows failed
at least one predefined check. Pooled body tilt had rank R-hat 1.1915, bulk ESS
59.2 and MCSE 0.663 degrees. Radius of gyration rose throughout the retained
windows, with rank R-hat 1.8963. See the
[full optimization and validation report](outputs/sampling_optimization_20260916/README.md).

Joint center scaling includes the `3*(N-1)*log(scale)` MH volume correction.
The fixed central position leaves 78 free position coordinates. This proposal
preserves the existing Cartesian target measure without adding an entropy term
to the energy. The front-cell position switch disables it. Explicit CLI proposal
options override the front-cell defaults and are recorded in each run manifest.

Both modes save per-chain checkpoints, source hashes, diagnostics, acceptance
rates and early/late window means. Add `--resume` to resume the same configuration.
The joint command retains the existing unconfined position space. Its output
characterizes finite-run mixing and does not establish unrestricted equilibrium.
The [Chinese diagnostic guide](docs/MCMC_DIAGNOSTICS_ZH.md) explains why large
R-hat does not prove that body tilt lacks a stationary distribution and discusses
the separate normalization problem for the full unbounded position space.

## Rebuilding the Notebook

From the repository root, run:

```bash
python _v903cells/build.py
```

## Notes

Ligand chains are not represented explicitly. Their separation and repulsion are incorporated through the effective gap and steric interaction. Nearest-neighbor bond tilt was removed from the latest output because PCA provides the retained global superlattice orientation measure.

Local beta now denotes the directed dipole-to-body-[111] angle, from 0 to 180
degrees. The previous implementation reported a body-to-field angle under that
name. Body-SL mismatch now denotes signed body tilt minus SL PCA tilt. The
actual angle between those axes is retained as a separate trajectory.
Signed magnetization is recorded alongside the original absolute alignment.
The central position remains fixed, but the central body can receive local
rotations. Both mechanical switches now cover global proposals as well.

The magnetic, vdW, steric and configurational-entropy models have not been
reparameterized. The current cluster model has no explicit finite container.
Diagnostics do not establish unrestricted assembly equilibrium or justify
interpreting MC cycles as physical time.

The older PDFs directly under `outputs/` and the Supporting Information document
describe the historical 16 nm calculation. They have not been regenerated for
this sampling update, and their old angle labels must be interpreted accordingly.
The reports under `outputs/sampling_2000cycles_20260915/` use the current angle
definitions and include the completed single-chain and multi-chain experiment.
