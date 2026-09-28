# Stability: stab_bond36
4 chains x 2000 cycles; fixed warm-up 1000 cycles; 1000 retained draws per chain, no thinning.
## Pooled diagnostics (all starts)
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -555.8 | 8.928 | 1.031 | 155.7 | 452.9 | 0.7159 | n/a | 38.1 | 2.994 | high_rhat,low_ess |
| vdw_kBT | -3.094 | 0.8694 | 1.047 | 94.25 | 193.4 | 0.08895 | n/a | 93.88 | 2.901 | high_rhat,low_ess |
| local_beta_deg | 89.32 | 8.202 | 1.021 | 158.5 | 470.1 | 0.6509 | n/a | 28.75 | 1.43 | high_rhat,low_ess,high_mcse |
| body_tilt_deg | 57.66 | 21.82 | 1.008 | 485 | 659.2 | 1.028 | 0.9254 | 11.5 | 1.639 | high_mcse |
| body_order | 0.1804 | 0.07477 | 1.017 | 285.3 | 521.8 | 0.004514 | 0.004086 | 14.73 | 2.58 | high_rhat,low_ess |
| sl_pca_tilt_deg | 13.25 | 7.222 | 1.073 | 34.87 | 110.1 | 1.277 | n/a | 397.8 | 1.962 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9308 | 0.01326 | 1.004 | 903.7 | 2360 | 0.0004369 | 0.0004245 | 9.189 | 1.448 | passed |
| min_gap_nm | 1.964 | 0.1779 | 1.005 | 1049 | 1357 | 0.006588 | 0.006032 | 8.974 | 3.778 | passed |
| rg_nm | 83.56 | 4.882 | 1.752 | 6.139 | 25.62 | 1.973 | n/a | 459.8 | 5.252 | high_rhat,low_ess |
| n_bonds | 31.41 | 2.492 | 1.068 | 81.27 | 92.8 | 0.2825 | n/a | 209.4 | 2.429 | high_rhat,low_ess |

## Per-start diagnostics

### compact_aligned
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -556.5 | 8.371 | 1.036 | 93.43 | 363.3 | 0.8705 | 0.6601 | 18.35 | 2.199 | fewer_than_4_chains,high_rhat,low_ess |
| vdw_kBT | -3.107 | 0.7272 | 1.017 | 74.23 | 175.4 | 0.08437 | n/a | 31.72 | 1.404 | fewer_than_4_chains,high_rhat,low_ess |
| local_beta_deg | 89.64 | 8.344 | 1.026 | 83.55 | 237.3 | 0.9139 | n/a | 28.75 | 0.9211 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| body_tilt_deg | 59.37 | 20.44 | 1.005 | 282.8 | 528.7 | 1.259 | 1.343 | 8.685 | 1.639 | fewer_than_4_chains,low_ess,high_mcse |
| body_order | 0.1886 | 0.07809 | 1.009 | 123 | 339.7 | 0.007131 | 0.006867 | 14.3 | 1.538 | fewer_than_4_chains,low_ess |
| sl_pca_tilt_deg | 12.69 | 7.689 | 1.147 | 14.41 | 56.49 | 2.432 | n/a | 397.8 | 1.962 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| magnetization | 0.9309 | 0.01311 | 1.001 | 579.8 | 1054 | 0.0005456 | 0.0005444 | 3.814 | 1.448 | fewer_than_4_chains |
| min_gap_nm | 1.952 | 0.1603 | 1.005 | 634.6 | 898.8 | 0.006904 | 0.006805 | 3.946 | 3.778 | fewer_than_4_chains |
| rg_nm | 80.42 | 3.291 | 1.623 | 3.412 | 22.74 | 1.814 | n/a | 459.8 | 5.252 | fewer_than_4_chains,high_rhat,low_ess |
| n_bonds | 31.66 | 2.696 | 1.1 | 33.11 | 41.29 | 0.483 | n/a | 209.4 | 1.119 | fewer_than_4_chains,high_rhat,low_ess |

### expanded_tilted
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -555.1 | 9.399 | 1.029 | 68.43 | 253 | 1.136 | n/a | 38.1 | 2.994 | fewer_than_4_chains,high_rhat,low_ess |
| vdw_kBT | -3.081 | 0.9913 | 1.077 | 29.05 | 70.4 | 0.1782 | n/a | 93.88 | 2.901 | fewer_than_4_chains,high_rhat,low_ess |
| local_beta_deg | 89 | 8.045 | 1.019 | 64.61 | 214.6 | 1 | n/a | 28.27 | 1.43 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| body_tilt_deg | 55.96 | 22.99 | 1.009 | 219.2 | 370 | 1.605 | 1.477 | 11.5 | 0.5991 | fewer_than_4_chains,low_ess,high_mcse |
| body_order | 0.1722 | 0.07036 | 1.021 | 168.2 | 255.4 | 0.00551 | 0.004714 | 14.73 | 2.58 | fewer_than_4_chains,high_rhat,low_ess |
| sl_pca_tilt_deg | 13.8 | 6.677 | 1.039 | 19.05 | 71.17 | 1.545 | n/a | 149.3 | 1.657 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| magnetization | 0.9307 | 0.0134 | 1.008 | 357.2 | 1255 | 0.0007012 | 0.0006466 | 9.189 | 0.7102 | fewer_than_4_chains,low_ess |
| min_gap_nm | 1.977 | 0.193 | 1.004 | 516.4 | 689.9 | 0.01055 | 0.00991 | 8.974 | 3.49 | fewer_than_4_chains |
| rg_nm | 86.69 | 4.143 | 1.219 | 7.498 | 24.21 | 1.595 | n/a | 177.8 | 4.22 | fewer_than_4_chains,high_rhat,low_ess |
| n_bonds | 31.16 | 2.243 | 1.041 | 50.66 | 87.7 | 0.3149 | n/a | 35.32 | 2.429 | fewer_than_4_chains,high_rhat,low_ess |

## Warm-up check: MSER-5 truncation point per chain (cycles)
| chain | mser_energy_kBT | mser_rg_nm | mser_body_tilt_deg | mser_vdw_kBT |
|---|---|---|---|---|
| compact_aligned_s21 | 0 | 915 | 0 | 195 |
| compact_aligned_s22 | 5 | 180 | 5 | 100 |
| expanded_tilted_s21 | 55 | 810 | 0 | 335 |
| expanded_tilted_s22 | 70 | 225 | 0 | 40 |

## Stationarity: first vs second half of the kept window (chain-averaged)
| start | variable | first_half | second_half | change |
|---|---|---|---|---|
| compact_aligned | energy_kBT | -556.8 | -556.3 | 0.4392 |
| compact_aligned | vdw_kBT | -3.15 | -3.064 | 0.08641 |
| compact_aligned | body_tilt_deg | 61.43 | 57.31 | -4.121 |
| compact_aligned | sl_pca_tilt_deg | 14.07 | 11.32 | -2.757 |
| compact_aligned | rg_nm | 79.89 | 80.95 | 1.057 |
| compact_aligned | magnetization | 0.9303 | 0.9315 | 0.001131 |
| compact_aligned | local_beta_deg | 88.95 | 90.33 | 1.371 |
| compact_aligned | body_order | 0.1909 | 0.1863 | -0.004645 |
| compact_aligned | n_bonds | 32.06 | 31.26 | -0.801 |
| compact_aligned | min_gap_nm | 1.938 | 1.965 | 0.0265 |
| expanded_tilted | energy_kBT | -556.6 | -553.6 | 3.029 |
| expanded_tilted | vdw_kBT | -2.968 | -3.194 | -0.2257 |
| expanded_tilted | body_tilt_deg | 55.89 | 56.03 | 0.1354 |
| expanded_tilted | sl_pca_tilt_deg | 12.87 | 14.74 | 1.865 |
| expanded_tilted | rg_nm | 86.19 | 87.2 | 1.009 |
| expanded_tilted | magnetization | 0.932 | 0.9293 | -0.002696 |
| expanded_tilted | local_beta_deg | 89.55 | 88.44 | -1.111 |
| expanded_tilted | body_order | 0.1592 | 0.1852 | 0.02599 |
| expanded_tilted | n_bonds | 31.06 | 31.26 | 0.198 |
| expanded_tilted | min_gap_nm | 1.98 | 1.973 | -0.006684 |

## Acceptance, constraint rejections, energy drift
| chain | minutes | acc_trans | acc_rot | acc_dip | acc_cotilt | acc_gamma | acc_scale | blocked_trans | blocked_rot | blocked_dip | blocked_cotilt | blocked_gamma | blocked_scale | drift_kBT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| compact_aligned_s21 | 25.3 | 0.637 | 0.6335 | 0.273 | 0.3935 | 0.708 | 0.208 | 0.1523 | 0.01757 | 0 | 0 | 0 | 0.2695 | 9.357e-14 |
| compact_aligned_s22 | 25.65 | 0.6402 | 0.6357 | 0.2733 | 0.4185 | 0.745 | 0.2165 | 0.15 | 0.01798 | 0 | 0 | 0 | 0.266 | 9.357e-14 |
| expanded_tilted_s21 | 25.73 | 0.6446 | 0.6334 | 0.2729 | 0.397 | 0.724 | 0.2285 | 0.1482 | 0.01726 | 0 | 0 | 0 | 0.236 | 9.357e-14 |
| expanded_tilted_s22 | 25.57 | 0.6382 | 0.6341 | 0.2712 | 0.404 | 0.7325 | 0.221 | 0.1559 | 0.01554 | 0 | 0 | 0.001 | 0.282 | 9.357e-14 |

Figures: traces.png, distributions_by_start.png, autocorrelation.png
