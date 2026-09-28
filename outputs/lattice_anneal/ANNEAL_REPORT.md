# Simulated annealing: best state of every restart (sorted by energy)

| run | best_cycle | energy_kBT | sl_axis_tilt_deg | body_tilt_deg | body_rotation_deg | body_rotation_max_deg | body_lattice_angle_deg | magnetization | local_beta_deg | zeeman_kBT | anisotropy_kBT | dipole_kBT | vdw_kBT |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lat_tilt40_s41 | 1500 | -700.9 | 0.2445 | 0.3363 | 19.8 | 33.16 | 1.571 | 0.9965 | 3.959 | -381.5 | -176.8 | -87.35 | -55.23 |
| lat_tilt60_s42 | 1500 | -700.3 | 0.05833 | 0.08097 | 14.98 | 32.18 | 0.7283 | 0.9968 | 3.795 | -381.6 | -176.9 | -87.34 | -54.5 |
| lat_tilt20_s41 | 1500 | -700 | 0.315 | 0.3893 | 15.21 | 32.22 | 2.053 | 0.9964 | 3.841 | -381.5 | -176.9 | -87.37 | -54.33 |
| lat_tilt40_s42 | 1500 | -700 | 0.2703 | 0.1507 | 17.88 | 32.52 | 1.41 | 0.9965 | 3.892 | -381.5 | -176.9 | -87.42 | -54.21 |
| lat_tilt60_s41 | 1498 | -699.8 | 0.1964 | 0.0982 | 11.84 | 31.54 | 0.992 | 0.9967 | 3.757 | -381.6 | -177 | -87.38 | -53.91 |
| lat_tilt00_s42 | 1496 | -699 | 0.1812 | 0.2753 | 15.42 | 31.15 | 2.115 | 0.9964 | 3.859 | -381.4 | -176.9 | -87.4 | -53.32 |
| lat_tilt20_s42 | 1500 | -698.9 | 0.1015 | 0.2504 | 12.22 | 32.44 | 1.528 | 0.9967 | 3.812 | -381.6 | -177 | -87.2 | -53.12 |
| lat_tilt00_s41 | 1500 | -698.7 | 0.09311 | 0.1491 | 18.26 | 32.1 | 1.842 | 0.9967 | 3.802 | -381.6 | -176.9 | -87.15 | -52.97 |

Energy spread over restarts: 2.269 kBT

Local test of the best state (tiny downhill-only perturbations): {"attempts": 400, "step_deg": 0.2, "accepted_downhill": 51, "total_drop_kBT": 0.020574204279645233, "largest_single_drop_kBT": 0.0021700236264320526, "drift_kBT": 2.9428832608366896e-10}
