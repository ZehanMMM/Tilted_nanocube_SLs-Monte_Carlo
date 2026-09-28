# Production v1 (archived): the energy hole

16 chains x 3000 cycles from the four lattice starts, target WITHOUT the
D0 = 0.165 nm sharp-core minimum.  Chain `compact_tilted40_s1` dropped in
steps to -462 kBT of vdW from cycle ~1300 on; at cycle 3000 a single pair
(16, 18) carried -458 kBT with sharp cores 1.5e-5 nm apart and the rounded
steric gap at 1.86 nm.

- `anomaly_pair_state.npz`: the cycle-3000 configuration of that chain and
  the six most attractive pairs.
- The value is the true sharp-cube integral: `v904.vdw.octree_energy`
  (value-blind volume refinement, no shared code) gives -10.97 kBT at 0.01 nm
  and -1.500 kBT at 0.05 nm for this geometry, identical to the surface form.
- `corner_search.py`: 600 random steric-wall contacts; with D0 the strongest
  attraction is -1.25 kBT, so the constraint removes only a sliver of
  configuration space.

These chains are NOT samples of the corrected target.  Their states (the
hole chain's at cycle 1300) seed production v2 via `scripts/make_v2_sources.py`.
