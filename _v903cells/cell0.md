# V9.03 -- N27 collective MC at B = 500 G

Final benchmark branch after V9.02 indicated that dipole unlocking can make
nonzero beta competitive.

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

The notebook now separates local cube rocking from coherent cube-body tilt. It
records the mean local cube `[111]` angle, the coherent cube-body `[111]` axis
tilt/order, and an independent superlattice/chain tilt from cube-center
geometry. The reported SL diagnostic is the PCA axis of all centers. The goal
is timing and a first N27 stability check, not final statistics.
