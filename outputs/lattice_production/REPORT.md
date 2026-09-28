# Rigid-lattice ensemble at 298 K: 16 chains x 2000 cycles (4 SL tilts x 4 seeds)
16 chains x 2000 cycles; fixed warm-up 500 cycles; 1500 retained draws per chain, no thinning.
## Pooled diagnostics (all starts)
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -660.3 | 5.943 | 1.006 | 6405 | 1.586e+04 | 0.074 | 0.06281 | 4.006 | 3.024 | passed |
| vdw_kBT | -46.25 | 0.8264 | 1.661 | 25.65 | 48.8 | 0.1672 | n/a | 1267 | 12.19 | high_rhat,low_ess |
| dipole_kBT | -85.13 | 2.734 | 1 | 1.761e+04 | 2.141e+04 | 0.02057 | 0.02033 | 1.555 | 2.696 | passed |
| local_beta_deg | 13.23 | 1.564 | 1.005 | 9752 | 1.721e+04 | 0.01587 | 0.01481 | 3.944 | 2.909 | passed |
| body_tilt_deg | 4.7 | 2.416 | 1.028 | 594.1 | 1021 | 0.09915 | n/a | 345.9 | 4.966 | high_rhat,low_ess |
| magnetization | 0.963 | 0.009372 | 1.004 | 6581 | 1.551e+04 | 0.0001136 | 0.0001022 | 4.757 | 2.662 | passed |
| sl_axis_tilt_deg | 4.706 | 2.424 | 1.029 | 575.7 | 1012 | 0.1015 | n/a | 337.5 | 5.006 | high_rhat,low_ess |
| body_rotation_deg | 14.46 | 2.739 | 1.839 | 23.31 | 49.54 | 0.577 | n/a | 1214 | 14.48 | high_rhat,low_ess |
| body_lattice_angle_deg | 2.804 | 1.367 | 2.173 | 20.64 | 30.95 | 0.3031 | n/a | 1396 | 8.1 | high_rhat,low_ess |
| anisotropy_kBT | -160.2 | 3.229 | 1.004 | 1.149e+04 | 1.83e+04 | 0.0301 | 0.02852 | 2.843 | 2.921 | passed |

## Per-start diagnostics

### lat_tilt00
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -660.7 | 5.97 | 1.005 | 2357 | 3659 | 0.123 | 0.1222 | 3.337 | 3.024 | passed |
| vdw_kBT | -46.41 | 0.8625 | 1.675 | 6.389 | 11.25 | 0.3625 | n/a | 1191 | 9.752 | high_rhat,low_ess |
| dipole_kBT | -85.11 | 2.709 | 1 | 4098 | 4977 | 0.04241 | 0.04394 | 1.522 | 2.119 | passed |
| local_beta_deg | 13.17 | 1.561 | 1.003 | 2793 | 3797 | 0.02974 | 0.02813 | 3.944 | 1.229 | passed |
| body_tilt_deg | 4.524 | 2.364 | 1.022 | 221.9 | 313.7 | 0.1618 | n/a | 36.38 | 1.853 | high_rhat,low_ess |
| magnetization | 0.9634 | 0.009249 | 1.002 | 2642 | 4134 | 0.0001801 | 0.0001832 | 3.753 | 1.204 | passed |
| sl_axis_tilt_deg | 4.524 | 2.365 | 1.021 | 223.4 | 313.1 | 0.1613 | n/a | 35.12 | 1.893 | high_rhat,low_ess |
| body_rotation_deg | 13.92 | 2.303 | 1.859 | 5.797 | 10.79 | 0.992 | n/a | 1210 | 11.54 | high_rhat,low_ess |
| body_lattice_angle_deg | 2.521 | 1.33 | 2.44 | 4.871 | 10.84 | 0.6035 | n/a | 1182 | 5.665 | high_rhat,low_ess |
| anisotropy_kBT | -160.3 | 3.24 | 1.002 | 3062 | 4009 | 0.05872 | 0.0544 | 2.792 | 0.7205 | passed |

### lat_tilt20
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -660 | 5.993 | 1.006 | 1722 | 4199 | 0.1441 | 0.1391 | 4.006 | 2.688 | passed |
| vdw_kBT | -46.28 | 0.8916 | 1.704 | 6.241 | 26.01 | 0.3623 | n/a | 1267 | 9.456 | high_rhat,low_ess |
| dipole_kBT | -85.1 | 2.745 | 1 | 4430 | 4983 | 0.04107 | 0.03849 | 1.525 | 1.132 | passed |
| local_beta_deg | 13.3 | 1.576 | 1.005 | 2483 | 3580 | 0.032 | 0.03203 | 3.225 | 2.909 | passed |
| body_tilt_deg | 4.878 | 2.46 | 1.033 | 137 | 129.5 | 0.2142 | n/a | 73.93 | 1.667 | high_rhat,low_ess |
| magnetization | 0.9626 | 0.009483 | 1.004 | 1578 | 2826 | 0.0002392 | 0.0002302 | 4.757 | 2.662 | passed |
| sl_axis_tilt_deg | 4.887 | 2.475 | 1.037 | 128.4 | 205.4 | 0.2222 | n/a | 82.96 | 1.503 | high_rhat,low_ess |
| body_rotation_deg | 14.46 | 3.065 | 1.817 | 5.896 | 20.01 | 1.33 | n/a | 1214 | 12.18 | high_rhat,low_ess |
| body_lattice_angle_deg | 3.157 | 1.576 | 1.89 | 5.686 | 12.46 | 0.6762 | n/a | 1289 | 5.497 | high_rhat,low_ess |
| anisotropy_kBT | -160.1 | 3.256 | 1.004 | 3003 | 4343 | 0.05936 | 0.06073 | 2.843 | 1.854 | passed |

### lat_tilt40
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -660.1 | 5.84 | 1.004 | 2452 | 4290 | 0.1163 | 0.1164 | 3.083 | 0.6659 | passed |
| vdw_kBT | -46.22 | 0.8417 | 1.573 | 6.855 | 11.34 | 0.3368 | n/a | 1170 | 5.316 | high_rhat,low_ess |
| dipole_kBT | -85.16 | 2.756 | 1.002 | 4379 | 5475 | 0.04153 | 0.04001 | 1.497 | 2.696 | passed |
| local_beta_deg | 13.26 | 1.548 | 1.004 | 2690 | 4405 | 0.02973 | 0.0281 | 2.659 | 2.337 | passed |
| body_tilt_deg | 4.812 | 2.407 | 1.047 | 114.5 | 297.9 | 0.2203 | n/a | 345.9 | 1.717 | high_rhat,low_ess |
| magnetization | 0.9627 | 0.009343 | 1.003 | 2160 | 3848 | 0.0002018 | 0.000194 | 3.363 | 1.11 | passed |
| sl_axis_tilt_deg | 4.815 | 2.403 | 1.045 | 118.7 | 276.4 | 0.2168 | n/a | 337.5 | 1.661 | high_rhat,low_ess |
| body_rotation_deg | 14.32 | 3.009 | 2.03 | 5.414 | 23.34 | 1.36 | n/a | 1036 | 5.573 | high_rhat,low_ess |
| body_lattice_angle_deg | 2.962 | 1.304 | 2.186 | 5.128 | 11.33 | 0.5914 | n/a | 1396 | 8.1 | high_rhat,low_ess |
| anisotropy_kBT | -160.2 | 3.179 | 1.003 | 3041 | 4173 | 0.05739 | 0.05277 | 2.188 | 2.921 | passed |

### lat_tilt60
| variable | mean | sd | rhat | ess_bulk | ess_tail | mcse_mean | mcse_batch | tau_int_max | geweke_max_abs | status |
|---|---|---|---|---|---|---|---|---|---|---|
| energy_kBT | -660.3 | 5.946 | 1.005 | 1989 | 3365 | 0.1359 | 0.1237 | 2.944 | 0.6868 | passed |
| vdw_kBT | -46.1 | 0.6587 | 1.702 | 6.267 | 11.61 | 0.2681 | n/a | 1154 | 12.19 | high_rhat,low_ess |
| dipole_kBT | -85.14 | 2.726 | 1 | 4445 | 5396 | 0.04105 | 0.03996 | 1.555 | 1.693 | passed |
| local_beta_deg | 13.16 | 1.567 | 1.007 | 2334 | 3868 | 0.03253 | 0.03004 | 2.671 | 2.321 | passed |
| body_tilt_deg | 4.586 | 2.412 | 1.021 | 144 | 305.1 | 0.2 | n/a | 85.67 | 4.966 | high_rhat,low_ess |
| magnetization | 0.9632 | 0.00939 | 1.006 | 1552 | 4117 | 0.0002329 | 0.0002075 | 3.705 | 2.505 | passed |
| sl_axis_tilt_deg | 4.597 | 2.435 | 1.026 | 132.4 | 276.1 | 0.2095 | n/a | 106 | 5.006 | high_rhat,low_ess |
| body_rotation_deg | 15.13 | 2.346 | 1.423 | 8.401 | 13.38 | 0.8174 | n/a | 731.7 | 14.48 | high_rhat,low_ess |
| body_lattice_angle_deg | 2.578 | 1.115 | 2.276 | 5.026 | 13.04 | 0.5082 | n/a | 1278 | 5.776 | high_rhat,low_ess |
| anisotropy_kBT | -160.3 | 3.234 | 1.006 | 2575 | 4362 | 0.0639 | 0.05982 | 2.501 | 2.404 | passed |

## Warm-up check: MSER-5 truncation point per chain (cycles)
| chain | mser_energy_kBT | mser_body_tilt_deg | mser_vdw_kBT | mser_sl_axis_tilt_deg |
|---|---|---|---|---|
| lat_tilt00_s31 | 5 | 160 | 775 | 160 |
| lat_tilt00_s32 | 10 | 25 | 620 | 25 |
| lat_tilt00_s33 | 5 | 5 | 0 | 5 |
| lat_tilt00_s34 | 5 | 25 | 995 | 25 |
| lat_tilt20_s31 | 25 | 25 | 145 | 25 |
| lat_tilt20_s32 | 5 | 10 | 0 | 10 |
| lat_tilt20_s33 | 55 | 580 | 995 | 580 |
| lat_tilt20_s34 | 10 | 145 | 5 | 145 |
| lat_tilt40_s31 | 20 | 35 | 940 | 35 |
| lat_tilt40_s32 | 25 | 20 | 50 | 20 |
| lat_tilt40_s33 | 145 | 65 | 195 | 65 |
| lat_tilt40_s34 | 35 | 40 | 0 | 45 |
| lat_tilt60_s31 | 65 | 70 | 925 | 70 |
| lat_tilt60_s32 | 50 | 50 | 995 | 50 |
| lat_tilt60_s33 | 55 | 70 | 565 | 70 |
| lat_tilt60_s34 | 45 | 55 | 0 | 55 |

## Stationarity: first vs second half of the kept window (chain-averaged)
| start | variable | first_half | second_half | change |
|---|---|---|---|---|
| lat_tilt00 | energy_kBT | -660.3 | -661.1 | -0.7647 |
| lat_tilt00 | sl_axis_tilt_deg | 4.43 | 4.619 | 0.1888 |
| lat_tilt00 | body_rotation_deg | 12.9 | 14.95 | 2.057 |
| lat_tilt00 | body_tilt_deg | 4.411 | 4.638 | 0.2267 |
| lat_tilt00 | body_lattice_angle_deg | 1.97 | 3.072 | 1.101 |
| lat_tilt00 | magnetization | 0.9636 | 0.9631 | -0.0005757 |
| lat_tilt00 | local_beta_deg | 13.12 | 13.22 | 0.09547 |
| lat_tilt00 | dipole_kBT | -85.08 | -85.15 | -0.07182 |
| lat_tilt00 | anisotropy_kBT | -160.4 | -160.3 | 0.1015 |
| lat_tilt00 | vdw_kBT | -45.9 | -46.92 | -1.019 |
| lat_tilt20 | energy_kBT | -659.6 | -660.4 | -0.7248 |
| lat_tilt20 | sl_axis_tilt_deg | 4.775 | 4.999 | 0.2245 |
| lat_tilt20 | body_rotation_deg | 13.17 | 15.74 | 2.565 |
| lat_tilt20 | body_tilt_deg | 4.763 | 4.994 | 0.2312 |
| lat_tilt20 | body_lattice_angle_deg | 1.965 | 4.349 | 2.384 |
| lat_tilt20 | magnetization | 0.963 | 0.9621 | -0.000881 |
| lat_tilt20 | local_beta_deg | 13.21 | 13.39 | 0.1814 |
| lat_tilt20 | dipole_kBT | -85.07 | -85.13 | -0.05896 |
| lat_tilt20 | anisotropy_kBT | -160.3 | -159.9 | 0.3218 |
| lat_tilt20 | vdw_kBT | -45.62 | -46.95 | -1.325 |
| lat_tilt40 | energy_kBT | -659.8 | -660.4 | -0.6374 |
| lat_tilt40 | sl_axis_tilt_deg | 4.957 | 4.673 | -0.2834 |
| lat_tilt40 | body_rotation_deg | 13.6 | 15.04 | 1.443 |
| lat_tilt40 | body_tilt_deg | 4.961 | 4.663 | -0.2978 |
| lat_tilt40 | body_lattice_angle_deg | 2.11 | 3.814 | 1.704 |
| lat_tilt40 | magnetization | 0.963 | 0.9624 | -0.0005448 |
| lat_tilt40 | local_beta_deg | 13.18 | 13.35 | 0.1657 |
| lat_tilt40 | dipole_kBT | -85.16 | -85.15 | 0.001542 |
| lat_tilt40 | anisotropy_kBT | -160.3 | -160 | 0.3317 |
| lat_tilt40 | vdw_kBT | -45.63 | -46.81 | -1.179 |
| lat_tilt60 | energy_kBT | -660.3 | -660.4 | -0.1168 |
| lat_tilt60 | sl_axis_tilt_deg | 4.399 | 4.796 | 0.3972 |
| lat_tilt60 | body_rotation_deg | 14.74 | 15.52 | 0.7846 |
| lat_tilt60 | body_tilt_deg | 4.368 | 4.805 | 0.4364 |
| lat_tilt60 | body_lattice_angle_deg | 1.98 | 3.175 | 1.195 |
| lat_tilt60 | magnetization | 0.9636 | 0.9628 | -0.0008695 |
| lat_tilt60 | local_beta_deg | 13.07 | 13.26 | 0.1931 |
| lat_tilt60 | dipole_kBT | -85.11 | -85.16 | -0.05284 |
| lat_tilt60 | anisotropy_kBT | -160.5 | -160.2 | 0.3592 |
| lat_tilt60 | vdw_kBT | -45.72 | -46.47 | -0.7507 |

## Acceptance, constraint rejections, energy drift
| chain | minutes | acc_trans | acc_rot | acc_dip | acc_cotilt | acc_gamma | acc_scale | acc_latrot | blocked_trans | blocked_rot | blocked_dip | blocked_cotilt | blocked_gamma | blocked_scale | blocked_latrot | drift_kBT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lat_tilt00_s31 | 80.81 | n/a | 0.3218 | 0.2366 | 0.1915 | 0.1255 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.398e-08 |
| lat_tilt00_s32 | 75.06 | n/a | 0.2748 | 0.236 | 0.193 | 0.145 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4.865e-09 |
| lat_tilt00_s33 | 81.91 | n/a | 0.3398 | 0.2368 | 0.186 | 0.1095 | n/a | n/a | 0 | 9.259e-05 | 0 | 0 | 0.001 | 0 | 0 | 3.396e-10 |
| lat_tilt00_s34 | 78.37 | n/a | 0.2914 | 0.2352 | 0.1975 | 0.115 | n/a | n/a | 0 | 3.704e-05 | 0 | 0 | 0 | 0 | 0 | 1.079e-08 |
| lat_tilt20_s31 | 75.72 | n/a | 0.3063 | 0.2369 | 0.1895 | 0.1365 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 6.686e-09 |
| lat_tilt20_s32 | 74.19 | n/a | 0.3467 | 0.2385 | 0.1935 | 0.1185 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 8.735e-08 |
| lat_tilt20_s33 | 81.37 | n/a | 0.318 | 0.2368 | 0.197 | 0.1315 | n/a | n/a | 0 | 0.0002963 | 0 | 0 | 0 | 0 | 0 | 3.389e-08 |
| lat_tilt20_s34 | 80.35 | n/a | 0.3116 | 0.2357 | 0.1845 | 0.138 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.07e-09 |
| lat_tilt40_s31 | 79.8 | n/a | 0.298 | 0.2364 | 0.199 | 0.121 | n/a | n/a | 0 | 1.852e-05 | 0 | 0 | 0 | 0 | 0 | 3.597e-10 |
| lat_tilt40_s32 | 80.13 | n/a | 0.3099 | 0.2358 | 0.1925 | 0.1225 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.567e-08 |
| lat_tilt40_s33 | 78.23 | n/a | 0.2982 | 0.237 | 0.1995 | 0.111 | n/a | n/a | 0 | 1.852e-05 | 0 | 0 | 0 | 0 | 0 | 6.602e-08 |
| lat_tilt40_s34 | 78.22 | n/a | 0.3475 | 0.2361 | 0.1965 | 0.106 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 2.574e-09 |
| lat_tilt60_s31 | 81.75 | n/a | 0.3052 | 0.2374 | 0.197 | 0.1335 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 1.178e-07 |
| lat_tilt60_s32 | 79.92 | n/a | 0.3198 | 0.2366 | 0.189 | 0.1155 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 4.825e-09 |
| lat_tilt60_s33 | 80.33 | n/a | 0.2761 | 0.2364 | 0.1895 | 0.173 | n/a | n/a | 0 | 0.0001296 | 0 | 0 | 0 | 0 | 0 | 2.127e-10 |
| lat_tilt60_s34 | 80.28 | n/a | 0.326 | 0.2371 | 0.1875 | 0.104 | n/a | n/a | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 3.573e-08 |

Figures: traces.png, distributions_by_start.png, autocorrelation.png
