# Run configurations

| file | purpose | chains x cycles |
|---|---|---|
| pilot_relax.json | pilot 1: relax each lattice start; samples discarded | 4 x 300 |
| (scripts/pilot_tuning.py) | pilot 2: step ladder scored by MSJD; writes tuned_moves.json | 20 x 120 per round |
| production.json (v1, archived) | tuned steps frozen; 4 starts x 4 seeds; target WITHOUT the D0 core minimum; one chain fell into the sharp-corner energy hole | 16 x 3000 |
| production.json (v2, current) | continues v1 (as relaxation) under the corrected target, fresh seeds 11-14 | 16 x 4000 |
| stab_halfstep.json | all tuned mechanical/collective steps halved: only efficiency may change | 4 x 1500 |
| stab_bond36.json | cluster definition 30 -> 36 nm: is the constraint active? | 4 x 2000 |
| stab_inherited_vdw.json | V9.03 4^3 vdW, same ensemble and steps: what the vdW update changes | 4 x 2000 |
| (scripts/pair_binding.py) | one pair in the same ensemble, bond 30 and 36 nm | 2 x 20000 |

The stability runs continue v2 production chains compact_aligned and
expanded_tilted, seeds 11 and 12 (`init_from`, `init_seed_offset` = 10), with
their own sampling seeds 21 and 22.

Seeds never repeat between pilot, production and stability runs, so no
tuning sample can leak into a reported estimate.  The quadrature accuracy is
checked without new sampling by `scripts/quadrature_reweighting.py`.
