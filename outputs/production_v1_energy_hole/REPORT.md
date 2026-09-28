# V9.04 production: 16 chains, converged vdW, Stillinger cluster (bond 30 nm)
16 chains x 3000 cycles; fixed warm-up 1000 cycles; 2000 retained draws per chain, no thinning.
## Pooled diagnostics (all starts)
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -572.8 | 46.51 | 1.24 | 47.06 | 19.44 | 9.116 | n/a | 1078 | 5.686 | high_rhat,low_ess |
| vdw_kBT | -12.43 | 45.75 | 1.223 | 53.28 | 19.21 | 9.119 | n/a | 1096 | 3.771 | high_rhat,low_ess |
| local_beta_deg | 90.05 | 7.89 | 1.014 | 1092 | 2473 | 0.2389 | 0.2132 | 38.17 | 2.35 | high_rhat,low_ess |
| body_tilt_deg | 57.65 | 21.51 | 1.007 | 2910 | 6440 | 0.413 | 0.3941 | 20.12 | 2.232 | passed |
| body_order | 0.1753 | 0.07332 | 1.01 | 2246 | 3747 | 0.001581 | 0.001575 | 26.46 | 3.31 | passed |
| sl_pca_tilt_deg | 14.98 | 11.76 | 1.289 | 42.8 | 23.57 | 2.311 | n/a | 1400 | 5.555 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9288 | 0.01377 | 1.016 | 792.5 | 9317 | 0.0004922 | n/a | 87.32 | 3.549 | high_rhat,low_ess |
| min_gap_nm | 1.9 | 0.105 | 1.025 | 3767 | 1.166e+04 | 0.002415 | n/a | 41.11 | 2.939 | high_rhat |
| rg_nm | 61.75 | 4.428 | 1.855 | 23.01 | 27.59 | 0.9319 | n/a | 1649 | 10.21 | high_rhat,low_ess |
| n_bonds | 29.69 | 2.229 | 1.099 | 112.7 | 290.5 | 0.2092 | n/a | 1050 | 4.248 | high_rhat,low_ess |

## Per-start diagnostics

### compact_aligned
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -564.5 | 9.498 | 1.089 | 32.33 | 427.7 | 1.669 | n/a | 480.8 | 5.686 | high_rhat,low_ess |
| vdw_kBT | -4.794 | 1.127 | 1.093 | 36.04 | 36.31 | 0.2229 | n/a | 300.6 | 2.928 | high_rhat,low_ess |
| local_beta_deg | 89.58 | 7.573 | 1.005 | 359.7 | 693.3 | 0.3997 | 0.4051 | 25.92 | 1.291 | low_ess,high_mcse |
| body_tilt_deg | 59.16 | 21.02 | 1.002 | 885.7 | 1616 | 0.737 | 0.715 | 11.92 | 2.232 | high_mcse |
| body_order | 0.1695 | 0.07181 | 1.011 | 701.5 | 1084 | 0.00277 | 0.002882 | 26.46 | 1.717 | high_rhat |
| sl_pca_tilt_deg | 20.96 | 18.48 | 1.727 | 6.289 | 25.17 | 8.059 | n/a | 1400 | 4.013 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9286 | 0.01361 | 1.009 | 1109 | 3751 | 0.0004039 | 0.0003156 | 6.458 | 1.121 | passed |
| min_gap_nm | 1.905 | 0.1095 | 1.006 | 1903 | 3793 | 0.002533 | 0.002257 | 3.919 | 1.8 | passed |
| rg_nm | 61.56 | 3.707 | 1.601 | 6.736 | 29.01 | 1.428 | n/a | 1295 | 10.21 | high_rhat,low_ess |
| n_bonds | 29.26 | 2.041 | 1.1 | 31.99 | 154.7 | 0.3639 | n/a | 99.98 | 2.601 | high_rhat,low_ess |

### compact_tilted40
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -594.7 | 87.94 | 1.517 | 7.206 | 10.82 | 34.7 | n/a | 1078 | 1.962 | high_rhat,low_ess |
| vdw_kBT | -35.2 | 87.61 | 1.495 | 7.547 | 12.73 | 34.7 | n/a | 1096 | 2.4 | high_rhat,low_ess |
| local_beta_deg | 90.89 | 7.829 | 1.02 | 263.5 | 699.9 | 0.4829 | 0.4517 | 37.34 | 1.039 | high_rhat,low_ess,high_mcse |
| body_tilt_deg | 57.09 | 21.76 | 1.01 | 658.6 | 1550 | 0.8887 | 0.8315 | 20.12 | 2.006 | high_mcse |
| body_order | 0.1749 | 0.07263 | 1.009 | 489.8 | 869.8 | 0.003336 | 0.003178 | 21.75 | 3.31 | passed |
| sl_pca_tilt_deg | 12.98 | 7.074 | 1.137 | 21.96 | 55.73 | 1.675 | n/a | 1111 | 2.322 | high_rhat,low_ess,high_mcse |
| magnetization | 0.928 | 0.01403 | 1.017 | 905 | 2769 | 0.0005053 | n/a | 87.32 | 3.549 | high_rhat |
| min_gap_nm | 1.888 | 0.09609 | 1.15 | 276.3 | 1619 | 0.008029 | n/a | 41.11 | 2.939 | high_rhat,low_ess |
| rg_nm | 59.96 | 3.514 | 1.865 | 5.794 | 28.09 | 1.469 | n/a | 1502 | 9.846 | high_rhat,low_ess |
| n_bonds | 29.76 | 2.195 | 1.048 | 95.12 | 224.8 | 0.2216 | n/a | 74.4 | 0.4698 | high_rhat,low_ess |

### expanded_aligned
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -565.3 | 9.563 | 1.107 | 28.89 | 252.4 | 1.783 | n/a | 471.1 | 2.929 | high_rhat,low_ess |
| vdw_kBT | -4.877 | 1.063 | 1.102 | 60.82 | 119.4 | 0.1272 | n/a | 623.5 | 3.046 | high_rhat,low_ess |
| local_beta_deg | 89.27 | 8.194 | 1.006 | 275.4 | 739.7 | 0.4939 | 0.4675 | 28.89 | 0.2982 | low_ess,high_mcse |
| body_tilt_deg | 56.74 | 21.43 | 1.004 | 755.6 | 1690 | 0.7965 | 0.7387 | 12.26 | 0.8861 | high_mcse |
| body_order | 0.1814 | 0.07461 | 1.002 | 521.5 | 1228 | 0.00333 | 0.003379 | 19.08 | 0.5111 | passed |
| sl_pca_tilt_deg | 12.6 | 7.262 | 1.195 | 15.69 | 21.57 | 2.072 | n/a | 675.9 | 2.459 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9288 | 0.0137 | 1.017 | 378 | 2535 | 0.0006953 | 0.0004286 | 18.21 | 0.9832 | high_rhat,low_ess |
| min_gap_nm | 1.903 | 0.1038 | 1.002 | 2073 | 4090 | 0.002173 | 0.002158 | 4.334 | 2.089 | passed |
| rg_nm | 62.08 | 4.865 | 1.873 | 5.75 | 11.37 | 2.096 | n/a | 1649 | 3.565 | high_rhat,low_ess |
| n_bonds | 30.02 | 2.236 | 1.125 | 23.89 | 207.4 | 0.4768 | n/a | 1050 | 4.248 | high_rhat,low_ess |

### expanded_tilted
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -566.9 | 9.923 | 1.154 | 17.21 | 86.57 | 2.427 | n/a | 168.9 | 2.512 | high_rhat,low_ess |
| vdw_kBT | -4.835 | 1.217 | 1.134 | 22.81 | 34.59 | 0.2805 | n/a | 770.6 | 3.771 | high_rhat,low_ess |
| local_beta_deg | 90.46 | 7.847 | 1.015 | 248.7 | 601.6 | 0.5025 | 0.4482 | 38.17 | 2.35 | high_rhat,low_ess,high_mcse |
| body_tilt_deg | 57.6 | 21.75 | 1.011 | 612.4 | 1701 | 0.914 | 0.8175 | 17.44 | 2.035 | high_rhat,high_mcse |
| body_order | 0.1754 | 0.0737 | 1.016 | 561.2 | 724.6 | 0.003211 | 0.003116 | 19.62 | 1.496 | high_rhat |
| sl_pca_tilt_deg | 13.39 | 7.788 | 1.086 | 32.26 | 61.01 | 1.471 | n/a | 636 | 5.555 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9298 | 0.01366 | 1.024 | 151.3 | 1998 | 0.00112 | 0.0003621 | 11.69 | 3.129 | high_rhat,low_ess |
| min_gap_nm | 1.904 | 0.1089 | 1.006 | 1461 | 3972 | 0.002744 | 0.002604 | 16.23 | 2.771 | passed |
| rg_nm | 63.38 | 4.759 | 1.874 | 5.695 | 11.25 | 2.032 | n/a | 1553 | 4.976 | high_rhat,low_ess |
| n_bonds | 29.72 | 2.363 | 1.125 | 22.66 | 258.7 | 0.4804 | n/a | 520 | 2.117 | high_rhat,low_ess |

## Warm-up check: MSER-5 truncation point per chain (cycles)
| chain | mser_energy_kBT | mser_rg_nm | mser_body_tilt_deg | mser_vdw_kBT |
|---|---|---|---|---|
| compact_aligned_s1 | 340 | 1405 | 325 | 425 |
| compact_aligned_s2 | 610 | 1495 | 285 | 295 |
| compact_aligned_s3 | 240 | 540 | 280 | 290 |
| compact_aligned_s4 | 490 | 1485 | 235 | 270 |
| compact_tilted40_s1 | 0 | 1495 | 200 | 0 |
| compact_tilted40_s2 | 270 | 890 | 240 | 375 |
| compact_tilted40_s3 | 555 | 1480 | 235 | 415 |
| compact_tilted40_s4 | 265 | 960 | 300 | 1295 |
| expanded_aligned_s1 | 930 | 1080 | 95 | 325 |
| expanded_aligned_s2 | 320 | 415 | 160 | 445 |
| expanded_aligned_s3 | 80 | 860 | 130 | 130 |
| expanded_aligned_s4 | 70 | 835 | 125 | 200 |
| expanded_tilted_s1 | 1250 | 1495 | 140 | 440 |
| expanded_tilted_s2 | 125 | 1130 | 95 | 420 |
| expanded_tilted_s3 | 470 | 1190 | 140 | 115 |
| expanded_tilted_s4 | 315 | 975 | 85 | 100 |

## Stationarity: first vs second half of the kept window (chain-averaged)
| start | variable | first_half | second_half | change |
|---|---|---|---|---|
| compact_aligned | energy_kBT | -562.9 | -566.1 | -3.219 |
| compact_aligned | vdw_kBT | -4.699 | -4.889 | -0.1903 |
| compact_aligned | body_tilt_deg | 58.66 | 59.67 | 1.008 |
| compact_aligned | sl_pca_tilt_deg | 19.93 | 21.98 | 2.05 |
| compact_aligned | rg_nm | 59.7 | 63.43 | 3.725 |
| compact_aligned | magnetization | 0.9286 | 0.9286 | -4.676e-05 |
| compact_aligned | local_beta_deg | 89.52 | 89.65 | 0.1366 |
| compact_aligned | body_order | 0.171 | 0.168 | -0.00294 |
| compact_aligned | n_bonds | 29.11 | 29.4 | 0.2845 |
| compact_aligned | min_gap_nm | 1.906 | 1.905 | -0.001008 |
| compact_tilted40 | energy_kBT | -572.8 | -616.5 | -43.69 |
| compact_tilted40 | vdw_kBT | -13.27 | -57.13 | -43.86 |
| compact_tilted40 | body_tilt_deg | 55.86 | 58.32 | 2.465 |
| compact_tilted40 | sl_pca_tilt_deg | 10.98 | 14.97 | 3.981 |
| compact_tilted40 | rg_nm | 58.77 | 61.14 | 2.376 |
| compact_tilted40 | magnetization | 0.9267 | 0.9292 | 0.002473 |
| compact_tilted40 | local_beta_deg | 91.53 | 90.25 | -1.274 |
| compact_tilted40 | body_order | 0.1729 | 0.1769 | 0.00392 |
| compact_tilted40 | n_bonds | 29.7 | 29.81 | 0.1065 |
| compact_tilted40 | min_gap_nm | 1.889 | 1.887 | -0.001955 |
| expanded_aligned | energy_kBT | -564.4 | -566.2 | -1.781 |
| expanded_aligned | vdw_kBT | -4.786 | -4.967 | -0.1806 |
| expanded_aligned | body_tilt_deg | 56.57 | 56.91 | 0.3391 |
| expanded_aligned | sl_pca_tilt_deg | 9.925 | 15.28 | 5.354 |
| expanded_aligned | rg_nm | 59.82 | 64.35 | 4.522 |
| expanded_aligned | magnetization | 0.9286 | 0.929 | 0.0003442 |
| expanded_aligned | local_beta_deg | 89.04 | 89.51 | 0.4722 |
| expanded_aligned | body_order | 0.1821 | 0.1807 | -0.00142 |
| expanded_aligned | n_bonds | 30.37 | 29.67 | -0.693 |
| expanded_aligned | min_gap_nm | 1.902 | 1.904 | 0.002499 |
| expanded_tilted | energy_kBT | -565.2 | -568.6 | -3.409 |
| expanded_tilted | vdw_kBT | -4.62 | -5.051 | -0.4308 |
| expanded_tilted | body_tilt_deg | 56.75 | 58.44 | 1.687 |
| expanded_tilted | sl_pca_tilt_deg | 13.65 | 13.13 | -0.5231 |
| expanded_tilted | rg_nm | 61.5 | 65.26 | 3.755 |
| expanded_tilted | magnetization | 0.93 | 0.9296 | -0.0004867 |
| expanded_tilted | local_beta_deg | 89.96 | 90.95 | 0.9847 |
| expanded_tilted | body_order | 0.1677 | 0.1831 | 0.01539 |
| expanded_tilted | n_bonds | 29.52 | 29.93 | 0.4032 |
| expanded_tilted | min_gap_nm | 1.908 | 1.901 | -0.007012 |

## Acceptance, constraint rejections, energy drift
| chain | minutes | acc_trans | acc_rot | acc_dip | acc_cotilt | acc_gamma | acc_scale | blocked_trans | blocked_rot | blocked_dip | blocked_cotilt | blocked_gamma | blocked_scale | drift_kBT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| compact_aligned_s1 | 54.48 | 0.5138 | 0.5787 | 0.2694 | 0.391 | 0.639 | 0.1337 | 0.2649 | 0.04026 | 0 | 0 | 0 | 0.291 | 9.357e-14 |
| compact_aligned_s2 | 56.65 | 0.5226 | 0.5818 | 0.2683 | 0.3927 | 0.629 | 0.1373 | 0.2498 | 0.03483 | 0 | 0 | 0.0003333 | 0.2767 | 9.357e-14 |
| compact_aligned_s3 | 56.15 | 0.5118 | 0.5699 | 0.2713 | 0.3927 | 0.6047 | 0.128 | 0.2574 | 0.04325 | 0 | 0 | 0.0003333 | 0.307 | 9.357e-14 |
| compact_aligned_s4 | 55.35 | 0.5282 | 0.5852 | 0.267 | 0.414 | 0.646 | 0.1407 | 0.2558 | 0.03253 | 0 | 0 | 0.001667 | 0.312 | 9.357e-14 |
| compact_tilted40_s1 | 55.87 | 0.5001 | 0.5546 | 0.2708 | 0.4003 | 0.284 | 0.06567 | 0.2598 | 0.03609 | 0 | 0 | 0.139 | 0.5643 | 9.357e-14 |
| compact_tilted40_s2 | 56.48 | 0.5212 | 0.5807 | 0.2695 | 0.4027 | 0.645 | 0.1337 | 0.25 | 0.03473 | 0 | 0 | 0 | 0.297 | 9.357e-14 |
| compact_tilted40_s3 | 57.23 | 0.5199 | 0.5719 | 0.2684 | 0.377 | 0.5917 | 0.1187 | 0.2416 | 0.03953 | 0 | 0 | 0 | 0.2713 | 9.357e-14 |
| compact_tilted40_s4 | 55.09 | 0.5243 | 0.5816 | 0.2681 | 0.3823 | 0.6463 | 0.1353 | 0.2563 | 0.03479 | 0 | 0 | 0 | 0.315 | 9.357e-14 |
| expanded_aligned_s1 | 54.36 | 0.5377 | 0.5985 | 0.2719 | 0.399 | 0.6697 | 0.1507 | 0.2522 | 0.02575 | 0 | 0 | 0 | 0.325 | 9.357e-14 |
| expanded_aligned_s2 | 55.9 | 0.5351 | 0.5947 | 0.269 | 0.3977 | 0.6457 | 0.1273 | 0.2474 | 0.02725 | 0 | 0 | 0 | 0.3077 | 9.357e-14 |
| expanded_aligned_s3 | 55.6 | 0.5364 | 0.599 | 0.2704 | 0.3907 | 0.675 | 0.1403 | 0.2504 | 0.02281 | 0 | 0 | 0.0003333 | 0.3133 | 9.357e-14 |
| expanded_aligned_s4 | 53.68 | 0.5351 | 0.601 | 0.2689 | 0.4087 | 0.653 | 0.1403 | 0.2519 | 0.02414 | 0 | 0 | 0 | 0.328 | 9.357e-14 |
| expanded_tilted_s1 | 54.49 | 0.5396 | 0.5989 | 0.2686 | 0.392 | 0.653 | 0.1367 | 0.2538 | 0.0243 | 0 | 0 | 0.0006667 | 0.303 | 9.357e-14 |
| expanded_tilted_s2 | 56.37 | 0.5324 | 0.5875 | 0.269 | 0.3877 | 0.6157 | 0.12 | 0.2428 | 0.02963 | 0 | 0 | 0.002333 | 0.2977 | 9.357e-14 |
| expanded_tilted_s3 | 54.02 | 0.5271 | 0.5876 | 0.2685 | 0.407 | 0.6253 | 0.13 | 0.2515 | 0.02748 | 0 | 0 | 0.0003333 | 0.3057 | 9.357e-14 |
| expanded_tilted_s4 | 52.77 | 0.5439 | 0.6098 | 0.2723 | 0.3907 | 0.671 | 0.1567 | 0.2546 | 0.02115 | 0 | 0 | 0 | 0.3113 | 9.357e-14 |

Figures: traces.png, distributions_by_start.png, autocorrelation.png
