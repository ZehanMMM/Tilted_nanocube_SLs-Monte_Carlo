# Rigid-lattice ensemble at 298 K, continuation: 16 chains x 1500 cycles (all of the first 2000 cycles treated as warm-up)
16 chains x 1500 cycles; fixed warm-up 0 cycles; 1500 retained draws per chain, no thinning.
## Pooled diagnostics (all starts)
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -660.6 | 5.961 | 1.005 | 5789 | 1.54e+04 | 0.07794 | 0.06553 | 5.73 | 2.924 | passed |
| vdw_kBT | -48.06 | 1.007 | 1.664 | 25.85 | 26.33 | 0.1969 | n/a | 1293 | 9.593 | high_rhat,low_ess |
| dipole_kBT | -85.05 | 2.771 | 1 | 1.68e+04 | 2.122e+04 | 0.02136 | 0.02133 | 1.859 | 2.969 | passed |
| local_beta_deg | 13.61 | 1.631 | 1.01 | 1770 | 1.685e+04 | 0.03924 | 0.01607 | 7.037 | 3.159 | high_rhat |
| body_tilt_deg | 4.805 | 2.438 | 1.022 | 732.2 | 1248 | 0.09167 | n/a | 314.4 | 2.927 | high_rhat,low_ess |
| magnetization | 0.9612 | 0.009638 | 1.007 | 4151 | 1.453e+04 | 0.000145 | 0.0001061 | 7.015 | 3.54 | passed |
| sl_axis_tilt_deg | 4.814 | 2.447 | 1.025 | 708.4 | 987.7 | 0.09408 | n/a | 330.1 | 4.063 | high_rhat,low_ess |
| body_rotation_deg | 16.5 | 1.381 | 2.015 | 21.6 | 50.96 | 0.2971 | n/a | 1205 | 9.802 | high_rhat,low_ess |
| body_lattice_angle_deg | 5.667 | 1.604 | 2.263 | 20.21 | 29.88 | 0.3731 | n/a | 1290 | 9.93 | high_rhat,low_ess |
| anisotropy_kBT | -159.5 | 3.349 | 1.008 | 4916 | 1.721e+04 | 0.04554 | 0.02991 | 3.199 | 3.093 | passed |

## Per-start diagnostics

### lat_tilt00
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -661 | 5.962 | 1.001 | 2339 | 3681 | 0.1234 | 0.1232 | 2.83 | 2.311 | passed |
| vdw_kBT | -48.09 | 1.18 | 1.846 | 5.837 | 13.14 | 0.4974 | n/a | 1050 | 6.82 | high_rhat,low_ess |
| dipole_kBT | -85.08 | 2.725 | 1.001 | 3977 | 5123 | 0.04317 | 0.04208 | 1.809 | 0.884 | passed |
| local_beta_deg | 13.52 | 1.646 | 1.017 | 825.3 | 3934 | 0.0533 | 0.03121 | 3.187 | 3.159 | high_rhat |
| body_tilt_deg | 4.597 | 2.292 | 1.025 | 203.4 | 398.4 | 0.1642 | n/a | 37.63 | 1.822 | high_rhat,low_ess |
| magnetization | 0.9617 | 0.009716 | 1.014 | 1365 | 2887 | 0.0002629 | 0.0002098 | 2.833 | 2.205 | high_rhat |
| sl_axis_tilt_deg | 4.597 | 2.279 | 1.021 | 208.3 | 397.4 | 0.162 | n/a | 49.29 | 1.871 | high_rhat,low_ess |
| body_rotation_deg | 15.94 | 1.616 | 2.237 | 5.096 | 15.52 | 0.7285 | n/a | 562.7 | 1.254 | high_rhat,low_ess |
| body_lattice_angle_deg | 5.12 | 2.046 | 2.421 | 4.884 | 15.61 | 0.977 | n/a | 1228 | 4.951 | high_rhat,low_ess |
| anisotropy_kBT | -159.7 | 3.359 | 1.014 | 2024 | 4368 | 0.0779 | 0.05523 | 2.302 | 2.44 | high_rhat |

### lat_tilt20
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -660.5 | 5.821 | 1.006 | 1940 | 4163 | 0.1324 | 0.1275 | 2.952 | 2.924 | passed |
| vdw_kBT | -48.55 | 0.8894 | 1.3 | 10.35 | 52.48 | 0.2834 | n/a | 393.3 | 9.593 | high_rhat,low_ess |
| dipole_kBT | -85.05 | 2.797 | 1.001 | 4242 | 5228 | 0.04296 | 0.04005 | 1.754 | 2.066 | passed |
| local_beta_deg | 13.81 | 1.596 | 1.003 | 2747 | 4399 | 0.03053 | 0.03203 | 2.518 | 2.744 | passed |
| body_tilt_deg | 4.724 | 2.442 | 1.023 | 152.4 | 289.9 | 0.2007 | n/a | 50.77 | 2.035 | high_rhat,low_ess |
| magnetization | 0.9604 | 0.00935 | 1.001 | 2453 | 4459 | 0.0001884 | 0.0001937 | 2.866 | 2.026 | passed |
| sl_axis_tilt_deg | 4.741 | 2.46 | 1.029 | 147.5 | 237 | 0.2093 | n/a | 55.63 | 2.477 | high_rhat,low_ess |
| body_rotation_deg | 17.07 | 1.506 | 2.128 | 5.264 | 20.01 | 0.6919 | n/a | 1193 | 6.996 | high_rhat,low_ess |
| body_lattice_angle_deg | 6.848 | 0.7845 | 1.429 | 8.213 | 35.64 | 0.2776 | n/a | 606.1 | 8.767 | high_rhat,low_ess |
| anisotropy_kBT | -159.2 | 3.32 | 1.002 | 3263 | 4241 | 0.05819 | 0.06049 | 2.263 | 2.903 | passed |

### lat_tilt40
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -660.8 | 6.063 | 1.006 | 1884 | 3872 | 0.14 | 0.1351 | 5.73 | 2.248 | passed |
| vdw_kBT | -48.06 | 0.7392 | 1.384 | 9.197 | 14.01 | 0.2379 | n/a | 878.1 | 3.519 | high_rhat,low_ess |
| dipole_kBT | -85.07 | 2.803 | 1 | 4186 | 5054 | 0.04327 | 0.04661 | 1.859 | 2.969 | passed |
| local_beta_deg | 13.58 | 1.651 | 1.008 | 2105 | 3904 | 0.03588 | 0.03345 | 7.037 | 2.857 | passed |
| body_tilt_deg | 4.85 | 2.445 | 1.016 | 231.7 | 374.9 | 0.1602 | n/a | 32.28 | 2.367 | high_rhat,low_ess |
| magnetization | 0.9612 | 0.00976 | 1.006 | 1411 | 3540 | 0.0002543 | 0.0002253 | 7.015 | 3.54 | passed |
| sl_axis_tilt_deg | 4.858 | 2.438 | 1.02 | 235.5 | 369.9 | 0.1588 | 0.1756 | 27.82 | 2.588 | high_rhat,low_ess |
| body_rotation_deg | 16.47 | 1.132 | 1.689 | 6.331 | 28.13 | 0.444 | n/a | 1205 | 9.802 | high_rhat,low_ess |
| body_lattice_angle_deg | 5.735 | 1.066 | 1.817 | 5.89 | 16.75 | 0.4515 | n/a | 1177 | 7.692 | high_rhat,low_ess |
| anisotropy_kBT | -159.6 | 3.359 | 1.005 | 2628 | 3975 | 0.06543 | 0.06111 | 3.199 | 3.093 | passed |

### lat_tilt60
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -660.2 | 5.966 | 1.005 | 1801 | 4102 | 0.1376 | 0.138 | 4.854 | 2.63 | passed |
| vdw_kBT | -47.53 | 0.8973 | 1.661 | 6.506 | 15.49 | 0.3515 | n/a | 1293 | 7.43 | high_rhat,low_ess |
| dipole_kBT | -84.99 | 2.756 | 1 | 4569 | 5240 | 0.04064 | 0.0416 | 1.401 | 1.474 | passed |
| local_beta_deg | 13.53 | 1.612 | 1.005 | 2307 | 4085 | 0.03354 | 0.03182 | 2.703 | 2.646 | passed |
| body_tilt_deg | 5.052 | 2.544 | 1.043 | 105.8 | 284.9 | 0.2578 | n/a | 314.4 | 2.927 | high_rhat,low_ess |
| magnetization | 0.9614 | 0.009678 | 1.004 | 1494 | 4295 | 0.000243 | 0.0002191 | 3.383 | 1.976 | passed |
| sl_axis_tilt_deg | 5.06 | 2.579 | 1.049 | 89.35 | 276 | 0.283 | n/a | 330.1 | 4.063 | high_rhat,low_ess |
| body_rotation_deg | 16.53 | 0.9107 | 1.271 | 11.65 | 31.05 | 0.2635 | n/a | 803.8 | 4.61 | high_rhat,low_ess |
| body_lattice_angle_deg | 4.964 | 1.471 | 2.457 | 4.854 | 15.31 | 0.6846 | n/a | 1290 | 9.93 | high_rhat,low_ess |
| anisotropy_kBT | -159.7 | 3.331 | 1.005 | 2576 | 4453 | 0.06558 | 0.06217 | 2.495 | 2.399 | passed |

## Warm-up check: MSER-5 truncation point per chain (cycles)
| chain | mser_energy_kBT | mser_body_tilt_deg | mser_vdw_kBT | mser_sl_axis_tilt_deg |
|---|---|---|---|---|
| lat_tilt00_s51 | 0 | 35 | 70 | 35 |
| lat_tilt00_s52 | 10 | 0 | 745 | 0 |
| lat_tilt00_s53 | 0 | 0 | 0 | 0 |
| lat_tilt00_s54 | 10 | 25 | 0 | 5 |
| lat_tilt20_s51 | 35 | 0 | 120 | 0 |
| lat_tilt20_s52 | 0 | 0 | 240 | 0 |
| lat_tilt20_s53 | 10 | 0 | 405 | 0 |
| lat_tilt20_s54 | 0 | 5 | 250 | 0 |
| lat_tilt40_s51 | 0 | 0 | 0 | 0 |
| lat_tilt40_s52 | 45 | 0 | 55 | 0 |
| lat_tilt40_s53 | 0 | 0 | 745 | 0 |
| lat_tilt40_s54 | 0 | 25 | 0 | 25 |
| lat_tilt60_s51 | 0 | 5 | 0 | 5 |
| lat_tilt60_s52 | 0 | 0 | 0 | 0 |
| lat_tilt60_s53 | 0 | 730 | 680 | 730 |
| lat_tilt60_s54 | 0 | 0 | 520 | 0 |

## Stationarity: first vs second half of the kept window (chain-averaged)
| start | variable | first_half | second_half | change |
|---|---|---|---|---|
| lat_tilt00 | energy_kBT | -661.1 | -661 | 0.09515 |
| lat_tilt00 | sl_axis_tilt_deg | 4.399 | 4.795 | 0.3967 |
| lat_tilt00 | body_rotation_deg | 15.67 | 16.21 | 0.5457 |
| lat_tilt00 | body_tilt_deg | 4.393 | 4.801 | 0.408 |
| lat_tilt00 | body_lattice_angle_deg | 4.5 | 5.74 | 1.241 |
| lat_tilt00 | magnetization | 0.9619 | 0.9614 | -0.0004988 |
| lat_tilt00 | local_beta_deg | 13.46 | 13.58 | 0.1273 |
| lat_tilt00 | dipole_kBT | -85.15 | -85.02 | 0.1299 |
| lat_tilt00 | anisotropy_kBT | -159.8 | -159.5 | 0.2757 |
| lat_tilt00 | vdw_kBT | -47.84 | -48.34 | -0.4982 |
| lat_tilt20 | energy_kBT | -660.4 | -660.5 | -0.07407 |
| lat_tilt20 | sl_axis_tilt_deg | 4.743 | 4.739 | -0.004151 |
| lat_tilt20 | body_rotation_deg | 16.88 | 17.26 | 0.3833 |
| lat_tilt20 | body_tilt_deg | 4.762 | 4.685 | -0.07714 |
| lat_tilt20 | body_lattice_angle_deg | 6.584 | 7.112 | 0.5281 |
| lat_tilt20 | magnetization | 0.9605 | 0.9603 | -0.0002253 |
| lat_tilt20 | local_beta_deg | 13.78 | 13.84 | 0.06218 |
| lat_tilt20 | dipole_kBT | -85.09 | -85.02 | 0.07124 |
| lat_tilt20 | anisotropy_kBT | -159.2 | -159.1 | 0.1114 |
| lat_tilt20 | vdw_kBT | -48.38 | -48.72 | -0.3404 |
| lat_tilt40 | energy_kBT | -660.6 | -661 | -0.4386 |
| lat_tilt40 | sl_axis_tilt_deg | 4.825 | 4.89 | 0.06505 |
| lat_tilt40 | body_rotation_deg | 16.41 | 16.54 | 0.1283 |
| lat_tilt40 | body_tilt_deg | 4.836 | 4.864 | 0.02767 |
| lat_tilt40 | body_lattice_angle_deg | 5.313 | 6.157 | 0.8444 |
| lat_tilt40 | magnetization | 0.9613 | 0.9611 | -0.0001869 |
| lat_tilt40 | local_beta_deg | 13.59 | 13.56 | -0.02874 |
| lat_tilt40 | dipole_kBT | -85.1 | -85.04 | 0.05921 |
| lat_tilt40 | anisotropy_kBT | -159.6 | -159.6 | -0.003129 |
| lat_tilt40 | vdw_kBT | -47.77 | -48.34 | -0.5639 |
| lat_tilt60 | energy_kBT | -660.1 | -660.4 | -0.3346 |
| lat_tilt60 | sl_axis_tilt_deg | 5.204 | 4.917 | -0.287 |
| lat_tilt60 | body_rotation_deg | 16.64 | 16.42 | -0.22 |
| lat_tilt60 | body_tilt_deg | 5.31 | 4.793 | -0.5178 |
| lat_tilt60 | body_lattice_angle_deg | 4.228 | 5.701 | 1.473 |
| lat_tilt60 | magnetization | 0.9617 | 0.961 | -0.0006838 |
| lat_tilt60 | local_beta_deg | 13.43 | 13.62 | 0.1967 |
| lat_tilt60 | dipole_kBT | -85.04 | -84.95 | 0.0845 |
| lat_tilt60 | anisotropy_kBT | -159.9 | -159.5 | 0.3966 |
| lat_tilt60 | vdw_kBT | -47 | -48.07 | -1.078 |

## Acceptance, constraint rejections, energy drift
| chain | minutes | acc_trans | acc_rot | acc_dip | acc_cotilt | acc_gamma | acc_scale | acc_latrot | blocked_trans | blocked_rot | blocked_dip | blocked_cotilt | blocked_gamma | blocked_scale | blocked_latrot | drift_kBT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lat_tilt00_s51 | 36.49 | n/a | 0.4365 | 0.2384 | 0.194 | 0.054 | n/a | n/a | 0 | 4.938e-05 | 0 | 0 | 0 | 0 | 0 | 1.303e-10 |
| lat_tilt00_s52 | 35.68 | n/a | 0.3054 | 0.2372 | 0.1947 | 0.08933 | n/a | n/a | 0 | 0.0005432 | 0 | 0 | 0 | 0 | 0 | 2.489e-10 |
| lat_tilt00_s53 | 37.27 | n/a | 0.4985 | 0.2394 | 0.2053 | 0.042 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.895e-08 |
| lat_tilt00_s54 | 37.3 | n/a | 0.4266 | 0.2362 | 0.196 | 0.064 | n/a | n/a | 0 | 0.0001481 | 0 | 0 | 0 | 0 | 0 | 2.555e-11 |
| lat_tilt20_s51 | 37.05 | n/a | 0.4881 | 0.2381 | 0.184 | 0.05667 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.043e-10 |
| lat_tilt20_s52 | 36.93 | n/a | 0.4844 | 0.2395 | 0.196 | 0.04933 | n/a | n/a | 0 | 9.877e-05 | 0 | 0 | 0 | 0 | 0 | 4.215e-09 |
| lat_tilt20_s53 | 36.82 | n/a | 0.5199 | 0.2375 | 0.1853 | 0.05 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.183e-08 |
| lat_tilt20_s54 | 36.44 | n/a | 0.444 | 0.2371 | 0.19 | 0.05267 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.158e-09 |
| lat_tilt40_s51 | 36.82 | n/a | 0.4516 | 0.2383 | 0.216 | 0.06733 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.178e-10 |
| lat_tilt40_s52 | 36.93 | n/a | 0.4247 | 0.2364 | 0.1927 | 0.05467 | n/a | n/a | 0 | 0.0002222 | 0 | 0 | 0 | 0 | 0 | 4.241e-09 |
| lat_tilt40_s53 | 36.5 | n/a | 0.4038 | 0.2369 | 0.2007 | 0.04533 | n/a | n/a | 0 | 0.0006667 | 0 | 0 | 0 | 0 | 0 | 1.88e-09 |
| lat_tilt40_s54 | 36.03 | n/a | 0.4697 | 0.2382 | 0.2013 | 0.04933 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.22e-08 |
| lat_tilt60_s51 | 37.23 | n/a | 0.395 | 0.2375 | 0.2027 | 0.06 | n/a | n/a | 0 | 0.0007407 | 0 | 0 | 0.009333 | 0 | 0 | 2.3e-09 |
| lat_tilt60_s52 | 36.21 | n/a | 0.4486 | 0.2373 | 0.1933 | 0.06133 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 9.022e-09 |
| lat_tilt60_s53 | 36.93 | n/a | 0.3592 | 0.237 | 0.1807 | 0.07467 | n/a | n/a | 0 | 7.407e-05 | 0 | 0 | 0 | 0 | 0 | 1.026e-08 |
| lat_tilt60_s54 | 36.73 | n/a | 0.4731 | 0.2376 | 0.1987 | 0.056 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4.9e-09 |

Figures: traces.png, distributions_by_start.png, autocorrelation.png
