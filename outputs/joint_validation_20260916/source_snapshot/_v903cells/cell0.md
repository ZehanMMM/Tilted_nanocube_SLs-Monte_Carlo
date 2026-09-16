# V9.03 -- N27 collective MC at B = 500 G

Local-delta Metropolis-Hastings sampling of coupled mechanical and magnetic
degrees of freedom. MC cycles are sampling steps, not physical time.

Default settings:

- field: `B = 500 G = 0.05 T`;
- anisotropy: Singh-like `cubic_first_raw`;
- cluster: full `3x3x3` rhombohedral block, `N = 27`;
- MC: local-delta collective MC with independent cube translations, cube-body
  rotations, dipole rotations, plus a small global co-tilt proposal that rotates
  cube centers and cube bodies together;
- V6.20-style coherent gamma twist: an optional global move rotates all cube
  bodies about the coherent cube `[111]` axis while leaving cube centers fixed;
- run length: `2000 cycles`, with first `500` treated as equilibration;
- default starting ansatz: `16.0: {'d': 21.0, 'alpha': 74.2}` with
  `beta_init = 0 deg`. Edit `DEFAULT_SIZE_NM` for experimental-size runs, or
  set `RUN_PRESET` to `19p5_tilted_candidate` / `29p6_tilted_candidate` for
  the tilted candidates.
- default random seed: `DEFAULT_SEED = 1`, editable in the first code cell.

Body tilt measures body [111] relative to the field. Local beta measures the
magnetic dipole relative to its own body [111], using a directed 0-180 degree
angle. SL tilt uses PCA only. Body-SL mismatch is the signed difference of the
two tilts. The true axis separation is recorded separately.

Each cycle now includes five complete random-order magnetic sweeps. This is a
provisional setting. The optional 1/5/10-sweep pilot compares diagnostic quality
and effective samples per second without changing production settings.

The optional multi-chain experiment crosses three initial structures with four
seeds. It reports rank-normalized split R-hat, bulk/tail/mean ESS and MCSE of the
mean. Short or frozen chains are flagged rather than declared converged.
