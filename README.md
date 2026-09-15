# V9.03 Nanocube Collective Monte Carlo Simulation

This repository contains the V9.03 Monte Carlo model used to study the coupled reorientation of a compact 3 x 3 x 3 magnetite nanocube superlattice in a 500 G magnetic field.

## Model

The simulated system contains 27 cubes with a default inorganic edge length of 16 nm. The initial structure is compact and rhombohedrally tilted. Each Monte Carlo cycle attempts collective cube translation, cube-body rotation, and magnetic dipole rotation according to the switches in the first notebook cell.

The energy includes Zeeman, cubic magnetocrystalline anisotropy, dipole-dipole, van der Waals, and steric terms. The superlattice tilt is obtained from principal component analysis of the cube-center coordinates. Cube-body tilt is the angle between each cube's body [111] axis and the laboratory magnetic-field direction. Monte Carlo cycles are sampling steps and are not interpreted as physical time.

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
- Magnetic updates: 5 complete random-order sweeps per cycle, provisional
- Optional multi-chain design: 3 initial structures x 4 seeds

## Repository Structure

- `V9.03_N27_500G.ipynb`: executable notebook
- `_v903cells/`: editable notebook source cells and rebuild script
- `tests/`: short sampler, diagnostic and notebook regression tests
- `docs/MCMC_DIAGNOSTICS_ZH.md`: Chinese explanation of the experiment and diagnostics
- `outputs/`: representative 16 nm full-run and multi-seed PDF reports
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

`N_MAG_PER_CYCLE = 5` now means five complete random permutations of the
particles per cycle. Each particle receives one dipole proposal per sweep.
Magnetic proposals evaluate only Zeeman, anisotropy and dipolar terms.
The unchanged vdW and steric terms cancel for these proposals.

Set `RUN_MULTI_SEED_SCAN = True` to cross `MULTI_STARTS` with `MULTI_SEEDS`.
Defaults use seeds 1, 2, 3 and 4 with these starts:

- `compact_aligned`: the compact cluster at 0 degrees
- `compact_tilted`: the same compact geometry co-rotated by 20 degrees
- `expanded_tilted`: an 8% expansion, co-rotated by 40 degrees

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

Set `RUN_MAGNETIC_PILOT = True` to compare 1, 5 and 10 sweeps using
`PILOT_CYCLES` and `PILOT_EQUIL`. A candidate must pass the diagnostic gates
for every monitored observable before it can be recommended. Eligible candidates
are ranked by their worst bulk ESS per elapsed second, including warm-up time.
If none qualifies, the recommendation is absent. The pilot never changes
`N_MAG_PER_CYCLE` automatically. Five sweeps is a provisional default, not a
measured optimum.

## Interpreting Diagnostics

ArviZ 0.22.0 computes diagnostics on unthinned post-warm-up arrays shaped
`(chain, draw)`, both pooled across starts and separately within each start.
Chain means are not used as a substitute for the underlying samples.

- **Rank-normalized split R-hat** compares variation within and across split
  chains after rank normalization, including a folded diagnostic sensitive to
  scale differences. The default target is below 1.01.
- **Bulk, tail and mean ESS** estimate how much independent-sample information
  remains after autocorrelation. Tail ESS uses the 5% and 95% quantiles.
  Mean ESS is distinct from bulk ESS. The default minimum is 400.
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

The PDFs under `outputs/` and the Supporting Information document describe the
historical 16 nm calculation. They have not been regenerated for this sampling
update, and their old angle labels must be interpreted accordingly.
