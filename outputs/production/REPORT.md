# V9.04 production v2: 16 chains x 4000 cycles, converged vdW, D0 core minimum, Stillinger cluster (bond 30 nm)
16 chains x 4000 cycles; fixed warm-up 1000 cycles; 3000 retained draws per chain, no thinning.
## Pooled diagnostics (all starts)
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -570.6 | 10.15 | 1.131 | 79.86 | 182.6 | 1.134 | n/a | 884.9 | 4.849 | high_rhat,low_ess |
| vdw_kBT | -4.885 | 1.091 | 1.041 | 360.9 | 913.6 | 0.06013 | n/a | 235.7 | 3.092 | high_rhat,low_ess |
| local_beta_deg | 90.3 | 8.029 | 1.005 | 1684 | 3689 | 0.1957 | 0.1817 | 44.54 | 2.897 | passed |
| body_tilt_deg | 57.33 | 21.49 | 1.002 | 4401 | 9343 | 0.3365 | 0.3337 | 19.49 | 2.303 | passed |
| body_order | 0.176 | 0.07476 | 1.004 | 2879 | 4711 | 0.001425 | 0.001389 | 24.65 | 2.784 | passed |
| sl_pca_tilt_deg | 12.62 | 9.274 | 1.217 | 53.26 | 30.46 | 1.626 | n/a | 1890 | 2.77 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9306 | 0.01344 | 1.012 | 1397 | 1.56e+04 | 0.0003616 | 0.0001956 | 47.95 | 4.523 | high_rhat,low_ess |
| min_gap_nm | 1.898 | 0.1007 | 1.002 | 1.14e+04 | 2.524e+04 | 0.0009678 | 0.0008396 | 4.694 | 2.514 | passed |
| rg_nm | 70.98 | 8.645 | 2.449 | 19.37 | 21.42 | 2.035 | n/a | 2637 | 10.55 | high_rhat,low_ess |
| n_bonds | 29.46 | 2.17 | 1.07 | 160.4 | 373 | 0.176 | n/a | 889.9 | 6.156 | high_rhat,low_ess |

## Per-start diagnostics

### compact_aligned
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -569 | 9.767 | 1.119 | 23.92 | 176.9 | 2.003 | n/a | 660.2 | 4.551 | high_rhat,low_ess |
| vdw_kBT | -4.807 | 1.039 | 1.067 | 54.59 | 218.5 | 0.1538 | n/a | 130.8 | 3.092 | high_rhat,low_ess |
| local_beta_deg | 90.48 | 8.144 | 1.013 | 318.4 | 740 | 0.4563 | 0.4083 | 44.54 | 2.897 | high_rhat,low_ess,high_mcse |
| body_tilt_deg | 57.08 | 21.27 | 1.004 | 1027 | 3196 | 0.6827 | 0.684 | 19.49 | 2.203 | high_mcse |
| body_order | 0.1746 | 0.07389 | 1.003 | 698.4 | 1006 | 0.002868 | 0.002798 | 19.67 | 1.292 | passed |
| sl_pca_tilt_deg | 17.9 | 13.89 | 1.377 | 8.941 | 38.29 | 5.307 | n/a | 812.1 | 1.259 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9303 | 0.01344 | 1.009 | 837 | 4018 | 0.0004636 | 0.0003248 | 13.35 | 0.8323 | passed |
| min_gap_nm | 1.899 | 0.1004 | 1.002 | 3361 | 6020 | 0.001759 | 0.001587 | 4.576 | 2.104 | passed |
| rg_nm | 70.28 | 5.306 | 2.181 | 5.143 | 11.38 | 2.337 | n/a | 2360 | 6.956 | high_rhat,low_ess |
| n_bonds | 29.07 | 2.053 | 1.093 | 31.05 | 327.1 | 0.3769 | n/a | 787.4 | 6.156 | high_rhat,low_ess |

### compact_tilted40
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -570.9 | 9.91 | 1.026 | 122.5 | 506.2 | 0.8948 | n/a | 142.5 | 1.827 | high_rhat,low_ess |
| vdw_kBT | -4.929 | 1.157 | 1.031 | 90.24 | 322.8 | 0.1228 | n/a | 99.58 | 2.231 | high_rhat,low_ess |
| local_beta_deg | 90.34 | 8.019 | 1.003 | 465.8 | 1113 | 0.3719 | 0.3649 | 31.2 | 1.152 | high_mcse |
| body_tilt_deg | 57.56 | 21.55 | 1.003 | 1032 | 2435 | 0.7088 | 0.622 | 12.26 | 1.342 | high_mcse |
| body_order | 0.1782 | 0.0752 | 1.005 | 812.4 | 1367 | 0.002691 | 0.002691 | 21.2 | 2.784 | passed |
| sl_pca_tilt_deg | 11.09 | 6.334 | 1.153 | 19.11 | 42.86 | 1.564 | n/a | 1890 | 2.197 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9304 | 0.01341 | 1.005 | 1140 | 5509 | 0.0003998 | 0.00029 | 9.544 | 2.582 | passed |
| min_gap_nm | 1.899 | 0.1023 | 1.001 | 3387 | 7269 | 0.001969 | 0.001798 | 3.691 | 1.851 | passed |
| rg_nm | 68.73 | 5.377 | 1.861 | 5.824 | 12 | 2.309 | n/a | 2581 | 6.26 | high_rhat,low_ess |
| n_bonds | 29.57 | 2.286 | 1.09 | 30.8 | 110 | 0.4083 | n/a | 889.9 | 1.478 | high_rhat,low_ess |

### expanded_aligned
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -574.2 | 10.25 | 1.19 | 14.62 | 93.11 | 2.68 | n/a | 118.6 | 4.849 | high_rhat,low_ess |
| vdw_kBT | -4.949 | 1.112 | 1.024 | 144.8 | 214.7 | 0.09716 | n/a | 101.4 | 2.696 | high_rhat,low_ess |
| local_beta_deg | 89.74 | 7.801 | 1.003 | 478.9 | 1022 | 0.3567 | 0.3256 | 30.83 | 1.91 | high_mcse |
| body_tilt_deg | 57.09 | 21.39 | 1.001 | 1231 | 2323 | 0.6321 | 0.6107 | 12.06 | 2.303 | high_mcse |
| body_order | 0.1712 | 0.07421 | 1.002 | 792 | 1329 | 0.002676 | 0.002506 | 17.54 | 0.9895 | passed |
| sl_pca_tilt_deg | 10.43 | 5.977 | 1.201 | 16.43 | 55.71 | 1.595 | n/a | 686.8 | 2.77 | high_rhat,low_ess,high_mcse |
| magnetization | 0.932 | 0.01327 | 1.02 | 245 | 4085 | 0.0008351 | 0.0002775 | 7.387 | 4.523 | high_rhat,low_ess |
| min_gap_nm | 1.896 | 0.09903 | 1.002 | 3121 | 6969 | 0.001823 | 0.001621 | 4.694 | 2.073 | passed |
| rg_nm | 78.32 | 11.78 | 2.363 | 4.912 | 11.2 | 5.617 | n/a | 2637 | 10.55 | high_rhat,low_ess |
| n_bonds | 29.52 | 2.071 | 1.044 | 127.8 | 345.4 | 0.1824 | n/a | 649 | 3.063 | high_rhat,low_ess |

### expanded_tilted
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -568.4 | 9.651 | 1.125 | 20.72 | 53.73 | 2.115 | n/a | 884.9 | 3.896 | high_rhat,low_ess |
| vdw_kBT | -4.855 | 1.047 | 1.057 | 88.97 | 192 | 0.1197 | n/a | 235.7 | 1.976 | high_rhat,low_ess |
| local_beta_deg | 90.66 | 8.117 | 1.003 | 429.2 | 926.5 | 0.3923 | 0.3853 | 37.13 | 1.449 | high_mcse |
| body_tilt_deg | 57.58 | 21.75 | 1.001 | 1033 | 1631 | 0.7061 | 0.6586 | 15.7 | 1.886 | high_mcse |
| body_order | 0.1799 | 0.07542 | 1.004 | 605.1 | 1207 | 0.003166 | 0.003133 | 24.65 | 1.07 | passed |
| sl_pca_tilt_deg | 11.05 | 6.137 | 1.036 | 83.01 | 159.3 | 0.7148 | n/a | 234.1 | 2.731 | high_rhat,low_ess,high_mcse |
| magnetization | 0.9298 | 0.01354 | 1.01 | 1166 | 5372 | 0.0003959 | 0.0003939 | 47.95 | 3.608 | passed |
| min_gap_nm | 1.898 | 0.1012 | 1.003 | 3487 | 5965 | 0.001812 | 0.001703 | 3.713 | 2.514 | passed |
| rg_nm | 66.59 | 4.936 | 2.095 | 5.265 | 12.28 | 2.201 | n/a | 2524 | 3.791 | high_rhat,low_ess |
| n_bonds | 29.66 | 2.214 | 1.054 | 56.7 | 171 | 0.3211 | n/a | 718.8 | 1.085 | high_rhat,low_ess |

## Warm-up check: MSER-5 truncation point per chain (cycles)
| chain | mser_energy_kBT | mser_rg_nm | mser_body_tilt_deg | mser_vdw_kBT |
|---|---|---|---|---|
| compact_aligned_s11 | 0 | 435 | 0 | 0 |
| compact_aligned_s12 | 5 | 610 | 0 | 0 |
| compact_aligned_s13 | 105 | 0 | 0 | 65 |
| compact_aligned_s14 | 0 | 0 | 0 | 0 |
| compact_tilted40_s11 | 0 | 0 | 0 | 60 |
| compact_tilted40_s12 | 125 | 1660 | 0 | 1210 |
| compact_tilted40_s13 | 75 | 75 | 0 | 0 |
| compact_tilted40_s14 | 0 | 0 | 0 | 0 |
| expanded_aligned_s11 | 0 | 1995 | 0 | 1845 |
| expanded_aligned_s12 | 0 | 0 | 0 | 0 |
| expanded_aligned_s13 | 315 | 305 | 0 | 0 |
| expanded_aligned_s14 | 590 | 1995 | 0 | 0 |
| expanded_tilted_s11 | 0 | 45 | 0 | 0 |
| expanded_tilted_s12 | 0 | 0 | 0 | 830 |
| expanded_tilted_s13 | 0 | 0 | 0 | 0 |
| expanded_tilted_s14 | 185 | 0 | 15 | 165 |

## Stationarity: first vs second half of the kept window (chain-averaged)
| start | variable | first_half | second_half | change |
|---|---|---|---|---|
| compact_aligned | energy_kBT | -567.4 | -570.5 | -3.071 |
| compact_aligned | vdw_kBT | -4.778 | -4.835 | -0.05699 |
| compact_aligned | body_tilt_deg | 56.53 | 57.63 | 1.103 |
| compact_aligned | sl_pca_tilt_deg | 19.26 | 16.54 | -2.716 |
| compact_aligned | rg_nm | 68.42 | 72.15 | 3.731 |
| compact_aligned | magnetization | 0.9302 | 0.9305 | 0.0003548 |
| compact_aligned | local_beta_deg | 90.66 | 90.31 | -0.3507 |
| compact_aligned | body_order | 0.1772 | 0.1721 | -0.005056 |
| compact_aligned | n_bonds | 28.81 | 29.34 | 0.5257 |
| compact_aligned | min_gap_nm | 1.898 | 1.899 | 0.001259 |
| compact_tilted40 | energy_kBT | -570.7 | -571.1 | -0.3733 |
| compact_tilted40 | vdw_kBT | -4.81 | -5.049 | -0.2388 |
| compact_tilted40 | body_tilt_deg | 57.85 | 57.26 | -0.5894 |
| compact_tilted40 | sl_pca_tilt_deg | 11.96 | 10.22 | -1.744 |
| compact_tilted40 | rg_nm | 65.52 | 71.93 | 6.409 |
| compact_tilted40 | magnetization | 0.9305 | 0.9302 | -0.0002454 |
| compact_tilted40 | local_beta_deg | 90.18 | 90.49 | 0.3045 |
| compact_tilted40 | body_order | 0.1806 | 0.1758 | -0.004781 |
| compact_tilted40 | n_bonds | 29.69 | 29.45 | -0.2382 |
| compact_tilted40 | min_gap_nm | 1.899 | 1.899 | -0.0007833 |
| expanded_aligned | energy_kBT | -573.4 | -575.1 | -1.673 |
| expanded_aligned | vdw_kBT | -5.058 | -4.841 | 0.2168 |
| expanded_aligned | body_tilt_deg | 57.63 | 56.56 | -1.07 |
| expanded_aligned | sl_pca_tilt_deg | 12.04 | 8.819 | -3.224 |
| expanded_aligned | rg_nm | 74.03 | 82.62 | 8.595 |
| expanded_aligned | magnetization | 0.9319 | 0.9322 | 0.0002908 |
| expanded_aligned | local_beta_deg | 89.7 | 89.77 | 0.07543 |
| expanded_aligned | body_order | 0.1714 | 0.171 | -0.0004346 |
| expanded_aligned | n_bonds | 29.84 | 29.21 | -0.6317 |
| expanded_aligned | min_gap_nm | 1.895 | 1.897 | 0.002124 |
| expanded_tilted | energy_kBT | -565.7 | -571.2 | -5.444 |
| expanded_tilted | vdw_kBT | -4.752 | -4.958 | -0.2056 |
| expanded_tilted | body_tilt_deg | 57.01 | 58.15 | 1.142 |
| expanded_tilted | sl_pca_tilt_deg | 11.58 | 10.52 | -1.062 |
| expanded_tilted | rg_nm | 63.74 | 69.44 | 5.7 |
| expanded_tilted | magnetization | 0.9291 | 0.9305 | 0.001461 |
| expanded_tilted | local_beta_deg | 90.35 | 90.97 | 0.6167 |
| expanded_tilted | body_order | 0.1797 | 0.1802 | 0.0005053 |
| expanded_tilted | n_bonds | 29.67 | 29.65 | -0.01433 |
| expanded_tilted | min_gap_nm | 1.901 | 1.895 | -0.005921 |

## Acceptance, constraint rejections, energy drift
| chain | minutes | acc_trans | acc_rot | acc_dip | acc_cotilt | acc_gamma | acc_scale | blocked_trans | blocked_rot | blocked_dip | blocked_cotilt | blocked_gamma | blocked_scale | drift_kBT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| compact_aligned_s11 | 83.22 | 0.5278 | 0.6007 | 0.2706 | 0.3902 | 0.6412 | 0.1348 | 0.2765 | 0.02874 | 0 | 0 | 0.00025 | 0.374 | 9.357e-14 |
| compact_aligned_s12 | 85.86 | 0.5325 | 0.5973 | 0.2675 | 0.4148 | 0.656 | 0.13 | 0.2681 | 0.02807 | 0 | 0 | 0.001 | 0.3448 | 9.357e-14 |
| compact_aligned_s13 | 84.48 | 0.5307 | 0.6018 | 0.2678 | 0.4105 | 0.6685 | 0.136 | 0.2711 | 0.0246 | 0 | 0 | 0.0005 | 0.3345 | 9.357e-14 |
| compact_aligned_s14 | 87.18 | 0.5276 | 0.5919 | 0.2689 | 0.4005 | 0.628 | 0.1258 | 0.2614 | 0.03275 | 0 | 0 | 0 | 0.3185 | 9.357e-14 |
| compact_tilted40_s11 | 86.66 | 0.5318 | 0.5983 | 0.2691 | 0.3925 | 0.639 | 0.128 | 0.2576 | 0.03097 | 0 | 0 | 0.00075 | 0.3377 | 9.357e-14 |
| compact_tilted40_s12 | 84.95 | 0.5355 | 0.6039 | 0.2676 | 0.3947 | 0.671 | 0.1368 | 0.2663 | 0.02544 | 0 | 0 | 0 | 0.3247 | 9.357e-14 |
| compact_tilted40_s13 | 87.07 | 0.5363 | 0.5995 | 0.2686 | 0.3822 | 0.6485 | 0.138 | 0.2606 | 0.02773 | 0 | 0 | 0.001 | 0.3242 | 9.357e-14 |
| compact_tilted40_s14 | 85.94 | 0.5281 | 0.5989 | 0.2686 | 0.4005 | 0.6495 | 0.1222 | 0.2641 | 0.02779 | 0 | 0 | 0.0005 | 0.328 | 9.357e-14 |
| expanded_aligned_s11 | 87.08 | 0.531 | 0.5968 | 0.2672 | 0.404 | 0.637 | 0.1258 | 0.2643 | 0.0286 | 0 | 0 | 0.0005 | 0.332 | 1.871e-13 |
| expanded_aligned_s12 | 86.87 | 0.5333 | 0.6033 | 0.2685 | 0.4 | 0.6432 | 0.1323 | 0.2637 | 0.02811 | 0 | 0 | 0.00125 | 0.3295 | 1.871e-13 |
| expanded_aligned_s13 | 87.43 | 0.5371 | 0.5971 | 0.2666 | 0.3812 | 0.6597 | 0.137 | 0.2562 | 0.02724 | 0 | 0 | 0 | 0.3098 | 9.357e-14 |
| expanded_aligned_s14 | 87.71 | 0.5269 | 0.5941 | 0.2656 | 0.3965 | 0.6288 | 0.1305 | 0.2602 | 0.03345 | 0 | 0 | 0.0005 | 0.3145 | 1.871e-13 |
| expanded_tilted_s11 | 86.67 | 0.5299 | 0.5956 | 0.2682 | 0.3972 | 0.6435 | 0.119 | 0.2615 | 0.02817 | 0 | 0 | 0.00125 | 0.3292 | 9.357e-14 |
| expanded_tilted_s12 | 87.04 | 0.5373 | 0.6034 | 0.2693 | 0.395 | 0.6548 | 0.129 | 0.2607 | 0.02807 | 0 | 0 | 0.0015 | 0.3445 | 9.357e-14 |
| expanded_tilted_s13 | 85.91 | 0.537 | 0.6012 | 0.27 | 0.4027 | 0.6595 | 0.1398 | 0.2604 | 0.02832 | 0 | 0 | 0.00025 | 0.3468 | 9.357e-14 |
| expanded_tilted_s14 | 87.08 | 0.534 | 0.6005 | 0.2689 | 0.3955 | 0.6432 | 0.132 | 0.2591 | 0.02893 | 0 | 0 | 0.0005 | 0.3335 | 9.357e-14 |

Figures: traces.png, distributions_by_start.png, autocorrelation.png
