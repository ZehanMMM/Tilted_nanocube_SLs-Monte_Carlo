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

## Repository Structure

- `V9.03_N27_500G.ipynb`: executable notebook
- `_v903cells/`: editable notebook source cells and rebuild script
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

Simulation controls are grouped near the beginning of `_v903cells/cell1.py` and in the first notebook setup cell. Set `DEFAULT_SEED`, `N_CYCLES`, `N_EQUIL`, and the move switches there before running the notebook from top to bottom. Multi-seed calculations repeat the simulation independently for each listed seed and report statistics over the post-equilibration samples.

## Rebuilding the Notebook

From the repository root, run:

```bash
python _v903cells/build.py
```

## Notes

Ligand chains are not represented explicitly. Their separation and repulsion are incorporated through the effective gap and steric interaction. Nearest-neighbor bond tilt was removed from the latest output because PCA provides the retained global superlattice orientation measure.

The files in this repository correspond to the 16 nm V9.03 calculation. See the Supporting Information document for model equations, parameter provenance, and interpretation of the four full-run figures.

