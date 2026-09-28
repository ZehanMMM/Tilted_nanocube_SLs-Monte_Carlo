# Stability: V9.03 inherited 4^3 vdW
4 chains x 2000 cycles; fixed warm-up 1000 cycles; 1000 retained draws per chain, no thinning.
## Pooled diagnostics (all starts)
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -571.4 | 9.597 | 1.115 | 23.89 | 67.92 | 1.963 | n/a | 215.5 | 4.152 | high_rhat,low_ess |
| vdw_kBT | -3.392 | 0.5474 | 1.096 | 45.05 | 188.5 | 0.07789 | n/a | 268.3 | 2.935 | high_rhat,low_ess |
| local_beta_deg | 89.44 | 8.077 | 1.041 | 134.3 | 401.9 | 0.6919 | n/a | 69.55 | 4.747 | high_rhat,low_ess,high_mcse |
| body_tilt_deg | 56.45 | 21.87 | 1.007 | 328 | 1068 | 1.245 | 1.156 | 13.49 | 0.6059 | low_ess,high_mcse |
| body_order | 0.1729 | 0.07425 | 1.011 | 310.5 | 393.5 | 0.004368 | 0.004077 | 16.44 | 1.72 | high_rhat,low_ess |
| sl_pca_tilt_deg | 9.481 | 4.717 | 1.06 | 54.51 | 85.72 | 0.6505 | n/a | 123.9 | 1.866 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9324 | 0.01317 | 1.016 | 348 | 1726 | 0.0006785 | 0.0004414 | 6.769 | 2.92 | high_rhat,low_ess |
| min_gap_nm | 1.91 | 0.114 | 1.003 | 1156 | 1920 | 0.003498 | 0.00326 | 4.803 | 1.977 | passed |
| rg_nm | 75.54 | 8.941 | 2.422 | 4.911 | 13.46 | 4.351 | n/a | 886 | 10.71 | high_rhat,low_ess |
| n_bonds | 29.71 | 2.083 | 1.105 | 30.2 | 76.53 | 0.3792 | n/a | 30.15 | 2.568 | high_rhat,low_ess |

## Per-start diagnostics

### compact_aligned
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -568.5 | 9.576 | 1.101 | 16.52 | 134.3 | 2.347 | n/a | 215.5 | 4.152 | fewer_than_4_chains,high_rhat,low_ess |
| vdw_kBT | -3.393 | 0.5511 | 1.069 | 39.51 | 137 | 0.08864 | n/a | 135.2 | 1.674 | fewer_than_4_chains,high_rhat,low_ess |
| local_beta_deg | 89.96 | 8.089 | 1.048 | 77.58 | 230 | 0.9231 | n/a | 31.44 | 4.747 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| body_tilt_deg | 57.64 | 21.41 | 1.009 | 169.7 | 464.6 | 1.63 | 1.602 | 13.49 | 0.6059 | fewer_than_4_chains,low_ess,high_mcse |
| body_order | 0.1794 | 0.079 | 1.008 | 145.7 | 176.3 | 0.006873 | 0.006772 | 16.44 | 1.72 | fewer_than_4_chains,low_ess |
| sl_pca_tilt_deg | 10.53 | 4.906 | 1.055 | 42.79 | 88.56 | 0.7863 | n/a | 123.9 | 1.866 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| magnetization | 0.9311 | 0.01334 | 1.021 | 171.8 | 588.3 | 0.001014 | 0.0006618 | 6.365 | 2.503 | fewer_than_4_chains,high_rhat,low_ess |
| min_gap_nm | 1.91 | 0.1158 | 1.001 | 661.3 | 998.3 | 0.004773 | 0.004598 | 3.924 | 1.188 | fewer_than_4_chains |
| rg_nm | 67.97 | 2.492 | 1.105 | 10.57 | 35.21 | 0.7806 | n/a | 321.7 | 2.58 | fewer_than_4_chains,high_rhat,low_ess |
| n_bonds | 29.89 | 2.119 | 1.105 | 20.24 | 43.81 | 0.4611 | n/a | 28.16 | 1.715 | fewer_than_4_chains,high_rhat,low_ess |

### expanded_tilted
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -574.3 | 8.722 | 1.055 | 51.78 | 61.04 | 1.22 | n/a | 48.68 | 1.716 | fewer_than_4_chains,high_rhat,low_ess |
| vdw_kBT | -3.391 | 0.5437 | 1.155 | 11.2 | 64.36 | 0.1574 | n/a | 268.3 | 2.935 | fewer_than_4_chains,high_rhat,low_ess |
| local_beta_deg | 88.92 | 8.032 | 1.043 | 44.85 | 142.8 | 1.196 | n/a | 69.55 | 2.213 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| body_tilt_deg | 55.27 | 22.25 | 1.004 | 174.4 | 687.3 | 1.783 | 1.614 | 12.29 | 0.3512 | fewer_than_4_chains,low_ess,high_mcse |
| body_order | 0.1664 | 0.06855 | 1.013 | 174.7 | 274.2 | 0.005233 | 0.004701 | 12.08 | 0.5606 | fewer_than_4_chains,high_rhat,low_ess |
| sl_pca_tilt_deg | 8.435 | 4.271 | 1.011 | 43.79 | 113.4 | 0.6589 | n/a | 46.86 | 0.5127 | fewer_than_4_chains,high_rhat,low_ess,high_mcse |
| magnetization | 0.9338 | 0.01287 | 1.004 | 355.2 | 1087 | 0.000676 | 0.0006112 | 6.769 | 2.92 | fewer_than_4_chains,low_ess |
| min_gap_nm | 1.91 | 0.1121 | 1.007 | 616.3 | 1017 | 0.004989 | 0.004623 | 4.803 | 1.977 | fewer_than_4_chains |
| rg_nm | 83.11 | 6.261 | 2.484 | 2.466 | 21.44 | 4.287 | n/a | 886 | 10.71 | fewer_than_4_chains,high_rhat,low_ess |
| n_bonds | 29.53 | 2.031 | 1.127 | 13.37 | 51.14 | 0.542 | n/a | 30.15 | 2.568 | fewer_than_4_chains,high_rhat,low_ess |

## Warm-up check: MSER-5 truncation point per chain (cycles)
| chain | mser_energy_kBT | mser_rg_nm | mser_body_tilt_deg | mser_vdw_kBT |
|---|---|---|---|---|
| compact_aligned_s21 | 440 | 215 | 0 | 0 |
| compact_aligned_s22 | 0 | 0 | 30 | 0 |
| expanded_tilted_s21 | 0 | 995 | 0 | 140 |
| expanded_tilted_s22 | 80 | 200 | 0 | 0 |

## Stationarity: first vs second half of the kept window (chain-averaged)
| start | variable | first_half | second_half | change |
|---|---|---|---|---|
| compact_aligned | energy_kBT | -565.7 | -571.3 | -5.583 |
| compact_aligned | vdw_kBT | -3.264 | -3.522 | -0.2585 |
| compact_aligned | body_tilt_deg | 59.79 | 55.49 | -4.3 |
| compact_aligned | sl_pca_tilt_deg | 11.1 | 9.95 | -1.153 |
| compact_aligned | rg_nm | 67.55 | 68.39 | 0.8392 |
| compact_aligned | magnetization | 0.9306 | 0.9316 | 0.0009422 |
| compact_aligned | local_beta_deg | 88.2 | 91.72 | 3.527 |
| compact_aligned | body_order | 0.1845 | 0.1744 | -0.01014 |
| compact_aligned | n_bonds | 29.7 | 30.08 | 0.389 |
| compact_aligned | min_gap_nm | 1.914 | 1.907 | -0.007106 |
| expanded_tilted | energy_kBT | -573.1 | -575.4 | -2.364 |
| expanded_tilted | vdw_kBT | -3.304 | -3.478 | -0.1744 |
| expanded_tilted | body_tilt_deg | 55.34 | 55.19 | -0.1462 |
| expanded_tilted | sl_pca_tilt_deg | 8.597 | 8.273 | -0.3236 |
| expanded_tilted | rg_nm | 80.99 | 85.22 | 4.224 |
| expanded_tilted | magnetization | 0.9333 | 0.9342 | 0.0008312 |
| expanded_tilted | local_beta_deg | 90.48 | 87.36 | -3.116 |
| expanded_tilted | body_order | 0.1697 | 0.1631 | -0.006588 |
| expanded_tilted | n_bonds | 29.61 | 29.45 | -0.151 |
| expanded_tilted | min_gap_nm | 1.919 | 1.901 | -0.01834 |

## Acceptance, constraint rejections, energy drift
| chain | minutes | acc_trans | acc_rot | acc_dip | acc_cotilt | acc_gamma | acc_scale | blocked_trans | blocked_rot | blocked_dip | blocked_cotilt | blocked_gamma | blocked_scale | drift_kBT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| compact_aligned_s21 | 4.795 | 0.5406 | 0.6056 | 0.2679 | 0.4045 | 0.674 | 0.1585 | 0.2628 | 0.02594 | 0 | 0 | 0.001 | 0.326 | 9.357e-14 |
| compact_aligned_s22 | 4.89 | 0.5412 | 0.6114 | 0.269 | 0.399 | 0.6635 | 0.154 | 0.2636 | 0.02252 | 0 | 0 | 0.001 | 0.3205 | 9.357e-14 |
| expanded_tilted_s21 | 4.628 | 0.5365 | 0.6063 | 0.2663 | 0.4 | 0.6615 | 0.14 | 0.2668 | 0.02394 | 0 | 0 | 0.001 | 0.325 | 9.357e-14 |
| expanded_tilted_s22 | 4.911 | 0.5451 | 0.6077 | 0.265 | 0.4065 | 0.699 | 0.1555 | 0.256 | 0.02209 | 0 | 0 | 0.0005 | 0.3125 | 9.357e-14 |

Figures: traces.png, distributions_by_start.png, autocorrelation.png
