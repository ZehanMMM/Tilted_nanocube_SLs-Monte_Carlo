# V9.03 完整 2000-cycle 与正式多链实验

## 结论

已完成 1 条单链和 12 条正式多链。合并诊断中，12/12 个观测量未通过预设检查。
所有轨迹及合并诊断值均为有限数值，没有 smoke test 中因样本不足产生的 N/A。
结果仍受初态或链间差异影响，2000 cycles 不足以支持整体平衡结论。

## 运行设置

- 单链：初始化 seed=1，采样 seed=1，保持 notebook 默认方案。
- 正式多链：3 种初态 × seeds 1、2、3、4，各自使用独立 SeedSequence 流。
- 每条链 2000 cycles，舍弃前 500，保留 1500。合并数组为 12 × 1500。
- 每 cycle 5 个完整磁 sweep。没有运行 1/5/10 sweep 的优化试跑。
- 27 个 16 nm cubes，500 G，298.15 K，cubic_first_raw。
- vdW、磁作用、排斥和熵模型保持 a64ce8d 中的定义与参数。
- 本次执行 wall time：25.21 min，4 个工作进程。
- 单链复现结果不混入正式多链诊断。初始随机流不同，因此它与 compact seed 1 不应逐点相同。

## 合并诊断

| 观测量 | rank split R-hat | bulk ESS | tail ESS | mean ESS | MCSE(mean) | 状态 |
|---|---:|---:|---:|---:|---:|---|
| energy_kBT | 1.746 | 18.076 | 30.044 | 13.537 | 5.4886 | high_rhat,low_ess |
| local_beta_deg | 1.3226 | 28.743 | 31.402 | 22.549 | 0.66724 | high_rhat,low_ess |
| body_tilt_deg | 1.8298 | 17.446 | 26.688 | 16.904 | 1.5013 | high_rhat,low_ess,high_mcse |
| sl_pca_tilt_deg | 1.8071 | 17.654 | 26.059 | 16.425 | 1.5605 | high_rhat,low_ess,high_mcse |
| body_sl_difference_deg | 1.46 | 23.894 | 43.41 | 23.351 | 0.19659 | high_rhat,low_ess |
| body_sl_axis_angle_deg | 1.5963 | 20.35 | 43.777 | 18.877 | 0.20902 | high_rhat,low_ess |
| magnetization | 1.2773 | 31.887 | 36.87 | 26.643 | 0.0029302 | high_rhat,low_ess |
| abs_muB | 1.2773 | 31.887 | 36.87 | 26.643 | 0.0029302 | high_rhat,low_ess |
| body_order | 2.7334 | 13.976 | 20.748 | 12.618 | 0.0031455 | high_rhat,low_ess |
| sl_pca_order | 2.9638 | 13.671 | 15.226 | 13.336 | 0.00097438 | high_rhat,low_ess |
| min_gap_nm | 1.1933 | 41.825 | 134.83 | 43.319 | 0.0028927 | high_rhat,low_ess |
| rg_nm | 2.8017 | 13.895 | 25.436 | 12.075 | 0.2981 | high_rhat,low_ess |

R-hat 目标小于 1.01，bulk/tail/mean ESS 最低 400。
body 与 SL tilt 的 MCSE 目标各为 0.5°。MCSE 采用观测量自身单位。
若链尚未混合，MCSE 不应当作可靠平衡估计的误差条。ESS 数值也依赖稳定采样的前提。

## 各初态内部的诊断

| 初态 | 观测量 | R-hat | bulk ESS | MCSE(mean) | 状态 |
|---|---|---:|---:|---:|---|
| compact_aligned | energy_kBT | 1.0633 | 41.655 | 0.95889 | high_rhat,low_ess |
| compact_aligned | body_tilt_deg | 1.3572 | 9.5921 | 0.71433 | high_rhat,low_ess,high_mcse |
| compact_aligned | sl_pca_tilt_deg | 1.3178 | 10.458 | 0.69373 | high_rhat,low_ess,high_mcse |
| compact_tilted | energy_kBT | 1.0612 | 51.99 | 0.83112 | high_rhat,low_ess |
| compact_tilted | body_tilt_deg | 1.5704 | 6.8861 | 1.0578 | high_rhat,low_ess,high_mcse |
| compact_tilted | sl_pca_tilt_deg | 1.5244 | 7.1636 | 1.0511 | high_rhat,low_ess,high_mcse |
| expanded_tilted | energy_kBT | 1.2176 | 13.046 | 3.2276 | high_rhat,low_ess |
| expanded_tilted | body_tilt_deg | 1.4398 | 8.2088 | 2.8671 | high_rhat,low_ess,high_mcse |
| expanded_tilted | sl_pca_tilt_deg | 1.4696 | 7.8738 | 3.0039 | high_rhat,low_ess,high_mcse |

## 每条链的 warm-up 后均值

以下是轨迹描述统计，不自动解释为平衡期望值。

| 初态 | seed | E / kBT | body tilt / deg | SL tilt / deg | local beta / deg | signed magnetization |
|---|---:|---:|---:|---:|---:|---:|
| compact_aligned | 1 | -637.87 | 3.0968 | 3.1842 | 14.024 | 0.96221 |
| compact_aligned | 2 | -637.72 | 4.7137 | 5.0316 | 13.807 | 0.96005 |
| compact_aligned | 3 | -637.25 | 5.304 | 4.9727 | 13.987 | 0.96146 |
| compact_aligned | 4 | -638.43 | 6.3618 | 6.3232 | 13.875 | 0.96081 |
| compact_tilted | 1 | -636.87 | 6.7267 | 6.4935 | 14.239 | 0.95889 |
| compact_tilted | 2 | -638.03 | 3.0444 | 3.0612 | 13.999 | 0.96164 |
| compact_tilted | 3 | -637.29 | 7.2742 | 6.8904 | 14.171 | 0.95937 |
| compact_tilted | 4 | -636.97 | 4.0847 | 4.3368 | 13.908 | 0.96142 |
| expanded_tilted | 1 | -601.58 | 11.424 | 11.433 | 16.947 | 0.94745 |
| expanded_tilted | 2 | -598.65 | 11.317 | 11.889 | 17.589 | 0.94519 |
| expanded_tilted | 3 | -599.98 | 13.445 | 12.917 | 16.83 | 0.94397 |
| expanded_tilted | 4 | -592.06 | 16.288 | 17.386 | 19.944 | 0.93552 |

## 独立单链复现

- mean energy：-637.279 kBT。
- mean body tilt：6.27305°。
- mean SL PCA tilt：6.47478°。
- mean local beta：14.018°。

## 数值检查

- 13 条链全部通过轨迹长度、有限性与最终完整能量核对。
- 最终全量能量减去累计能量的最大绝对误差：6.821e-13 kBT。
- 所有 cycle 中最小的现有投影间隙：1.79999 nm。
- 最终磁矩长度和姿态矩阵正交性通过检查。
- 现有投影间隙检查不替代真实凸体碰撞验证。
- source SHA256 与 manifest 一致。原始模型和采样代码在运行中没有变化。

## 文件

- V903_MultiChain_Report.pdf：多链配置、均值、诊断与分开的 body/SL 轨迹。
- V903_MultiChain_diagnostics.csv：所有初态合并及每种初态内部的诊断。
- V903_MultiChain_chains.csv：链身份、随机流、时长与描述统计。
- V903_MultiChain_trajectories.npz：完整 2000-cycle 轨迹及最终状态。
- V903_MultiChain_metadata.json：实际运行配置。
- single_seed1/：独立单链报告、四图、原始数据与验证。
- 各初态 seed 目录：逐链 checkpoint、进度日志和验证结果。

## 后续判断

先根据 R-hat、跨初态的均值差异与轨迹漂移判断慢变量，再决定是否延长采样或改 proposal。
本次没有据此调整磁作用、vdW 或构型熵，也没有自动延长链或将 5 sweeps 宣称为最优。
