# Stability: stab_halfstep
4 chains x 1500 cycles; fixed warm-up 300 cycles; 1200 retained draws per chain, no thinning.
## Pooled diagnostics (all starts)
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -569.2 | 9.232 | 1.103 | 27.84 | 115.5 | 1.745 | n/a | 197.3 | 1.31 | high_rhat,low_ess |
| vdw_kBT | -4.64 | 0.9997 | 1.066 | 61.51 | 139.5 | 0.1214 | n/a | 86.14 | 1.61 | high_rhat,low_ess |
| local_beta_deg | 91.01 | 8.737 | 1.104 | 29.27 | 83.02 | 1.614 | n/a | 126 | 1.835 | high_rhat,low_ess,high_mcse |
| body_tilt_deg | 55.1 | 22.27 | 1.022 | 88.68 | 323.4 | 2.444 | n/a | 181.2 | 0.8258 | high_rhat,low_ess,high_mcse |
| body_order | 0.182 | 0.07222 | 1.029 | 120.5 | 154.2 | 0.006729 | n/a | 63.08 | 2.852 | high_rhat,low_ess |
| sl_pca_tilt_deg | 10.74 | 6.167 | 1.219 | 14.46 | 21.09 | 1.779 | n/a | 351.5 | 3.77 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9308 | 0.01293 | 1.005 | 739.1 | 1981 | 0.0004854 | 0.0003956 | 6.618 | 3.233 | passed |
| min_gap_nm | 1.904 | 0.1079 | 1.006 | 1551 | 2739 | 0.002765 | 0.002609 | 3.381 | 0.8867 | passed |
| rg_nm | 67.55 | 4.11 | 2.321 | 4.998 | 27.19 | 1.9 | n/a | 943.3 | 7.676 | high_rhat,low_ess |
| n_bonds | 29.14 | 1.975 | 1.13 | 23.81 | 68.8 | 0.4233 | n/a | 449.2 | 3.905 | high_rhat,low_ess |

## Per-start diagnostics

### compact_aligned
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -566.6 | 8.769 | 1.044 | 34.31 | 194.8 | 1.484 | n/a | 197.3 | 0.6096 | fewer_than_4_chains,high_rhat,low_ess |
| vdw_kBT | -4.443 | 0.9178 | 1.04 | 39.11 | 91.28 | 0.1485 | n/a | 86.14 | 0.3696 | fewer_than_4_chains,high_rhat,low_ess |
| local_beta_deg | 88.45 | 8.349 | 1.056 | 35.91 | 78.85 | 1.4 | n/a | 110.1 | 1.835 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| body_tilt_deg | 55 | 22.57 | 1.009 | 52.16 | 188.5 | 3.163 | n/a | 43.92 | 0.1372 | fewer_than_4_chains,low_ess,high_mcse |
| body_order | 0.174 | 0.07009 | 1.009 | 75.88 | 151.3 | 0.008151 | n/a | 33.53 | 1.202 | fewer_than_4_chains,low_ess |
| sl_pca_tilt_deg | 12.5 | 6.809 | 1.307 | 5.546 | 26.92 | 3.001 | n/a | 351.5 | 3.77 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| magnetization | 0.9302 | 0.01294 | 1.001 | 548.2 | 1113 | 0.0005525 | 0.0005554 | 4.756 | 1.289 | fewer_than_4_chains |
| min_gap_nm | 1.91 | 0.1132 | 1.008 | 877.3 | 1075 | 0.004078 | 0.003851 | 3.381 | 0.8867 | fewer_than_4_chains |
| rg_nm | 64.49 | 2.373 | 1.774 | 2.999 | 20.01 | 1.404 | n/a | 660.9 | 7.676 | fewer_than_4_chains,high_rhat,low_ess |
| n_bonds | 28.72 | 1.778 | 1.055 | 31.72 | 82.72 | 0.3014 | n/a | 175 | 3.905 | fewer_than_4_chains,high_rhat,low_ess |

### expanded_tilted
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -571.8 | 8.954 | 1.079 | 19.73 | 176.7 | 2.033 | n/a | 149.4 | 1.31 | fewer_than_4_chains,high_rhat,low_ess |
| vdw_kBT | -4.837 | 1.038 | 1.063 | 45.57 | 80.57 | 0.1608 | n/a | 65.5 | 1.61 | fewer_than_4_chains,high_rhat,low_ess |
| local_beta_deg | 93.57 | 8.359 | 1.061 | 21.67 | 46.71 | 1.803 | n/a | 126 | 1.521 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| body_tilt_deg | 55.2 | 21.95 | 1.045 | 43.01 | 75.64 | 3.596 | n/a | 181.2 | 0.8258 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| body_order | 0.1899 | 0.07343 | 1.042 | 50.73 | 84.08 | 0.0106 | n/a | 63.08 | 2.852 | fewer_than_4_chains,high_rhat,low_ess |
| sl_pca_tilt_deg | 8.973 | 4.847 | 1.081 | 24.95 | 64.04 | 1.052 | n/a | 106.5 | 1.443 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| magnetization | 0.9315 | 0.0129 | 1.008 | 296.3 | 948.9 | 0.0007677 | 0.0005634 | 6.618 | 3.233 | fewer_than_4_chains,low_ess |
| min_gap_nm | 1.899 | 0.1019 | 1.003 | 716 | 1633 | 0.003612 | 0.003522 | 3.008 | 0.3172 | fewer_than_4_chains |
| rg_nm | 70.6 | 3.082 | 1.701 | 3.297 | 36.73 | 1.829 | n/a | 943.3 | 1.965 | fewer_than_4_chains,high_rhat,low_ess |
| n_bonds | 29.55 | 2.073 | 1.188 | 9.012 | 23.95 | 0.7326 | n/a | 449.2 | 2.889 | fewer_than_4_chains,high_rhat,low_ess |

## Warm-up check: MSER-5 truncation point per chain (cycles)
| chain | mser_energy_kBT | mser_rg_nm | mser_body_tilt_deg | mser_vdw_kBT |
|---|---|---|---|---|
| compact_aligned_s21 | 0 | 455 | 0 | 0 |
| compact_aligned_s22 | 0 | 10 | 0 | 0 |
| expanded_tilted_s21 | 45 | 0 | 0 | 175 |
| expanded_tilted_s22 | 0 | 0 | 0 | 0 |

## Stationarity: first vs second half of the kept window (chain-averaged)
| start | variable | first_half | second_half | change |
|---|---|---|---|---|
| compact_aligned | energy_kBT | -565 | -568.3 | -3.362 |
| compact_aligned | vdw_kBT | -4.472 | -4.413 | 0.05972 |
| compact_aligned | body_tilt_deg | 53.83 | 56.18 | 2.351 |
| compact_aligned | sl_pca_tilt_deg | 15 | 9.992 | -5.009 |
| compact_aligned | rg_nm | 64.16 | 64.82 | 0.661 |
| compact_aligned | magnetization | 0.9297 | 0.9307 | 0.0009737 |
| compact_aligned | local_beta_deg | 88.7 | 88.21 | -0.4916 |
| compact_aligned | body_order | 0.1799 | 0.1681 | -0.01188 |
| compact_aligned | n_bonds | 29.02 | 28.42 | -0.5958 |
| compact_aligned | min_gap_nm | 1.912 | 1.909 | -0.002979 |
| expanded_tilted | energy_kBT | -574.1 | -569.5 | 4.524 |
| expanded_tilted | vdw_kBT | -5.044 | -4.63 | 0.4141 |
| expanded_tilted | body_tilt_deg | 59.2 | 51.2 | -7.993 |
| expanded_tilted | sl_pca_tilt_deg | 9.37 | 8.577 | -0.7929 |
| expanded_tilted | rg_nm | 69.74 | 71.46 | 1.714 |
| expanded_tilted | magnetization | 0.931 | 0.9319 | 0.0009072 |
| expanded_tilted | local_beta_deg | 92.57 | 94.56 | 1.992 |
| expanded_tilted | body_order | 0.172 | 0.2078 | 0.0358 |
| expanded_tilted | n_bonds | 30.02 | 29.09 | -0.9283 |
| expanded_tilted | min_gap_nm | 1.897 | 1.901 | 0.004064 |

## Acceptance, constraint rejections, energy drift
| chain | minutes | acc_trans | acc_rot | acc_dip | acc_cotilt | acc_gamma | acc_scale | blocked_trans | blocked_rot | blocked_dip | blocked_cotilt | blocked_gamma | blocked_scale | drift_kBT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| compact_aligned_s21 | 23.91 | 0.672 | 0.6679 | 0.2694 | 0.6347 | 0.806 | 0.2627 | 0.1207 | 0.0202 | 0 | 0 | 0 | 0.2433 | 9.357e-14 |
| compact_aligned_s22 | 24.01 | 0.6731 | 0.672 | 0.2694 | 0.6373 | 0.8247 | 0.2807 | 0.1211 | 0.01963 | 0 | 0 | 0.0006667 | 0.2333 | 9.357e-14 |
| expanded_tilted_s21 | 25.14 | 0.6645 | 0.6603 | 0.2703 | 0.628 | 0.794 | 0.252 | 0.1163 | 0.02583 | 0 | 0 | 0.0006667 | 0.2147 | 9.357e-14 |
| expanded_tilted_s22 | 24.24 | 0.6672 | 0.671 | 0.2666 | 0.6227 | 0.814 | 0.264 | 0.1165 | 0.02096 | 0 | 0 | 0.001333 | 0.2053 | 0 |

Figures: traces.png, distributions_by_start.png, autocorrelation.png
