# V9.03 采样实验与诊断

本次更新只修改采样、初态设计、诊断与输出。vdW、磁相互作用和排斥势的公式及参数保持不变，未添加构型熵能量项。已有 PDF 和 SI 是历史结果，没有重新计算。

## 1. 磁更新数量

`N_MAG_PER_CYCLE = 5` 表示每个 MC cycle 做 5 个完整磁 sweep。每个 sweep 都按新的随机排列访问全部 27 个颗粒，每个颗粒尝试一次磁矩旋转。原先按有放回抽样选颗粒，单个颗粒可能在一个 cycle 内被访问多次或一次也没有被访问。

磁更新现在只计算 Zeeman、各向异性和偶极能量。磁矩改变时，位置与本体姿态不变，因此 vdW 与排斥项在能量差中抵消。这个优化不改变磁相互作用模型。

5 是暂定值，尚未通过正式多链试跑证明最优。设置 `RUN_MAGNETIC_PILOT = True`，可比较 `MAGNETIC_SWEEP_CANDIDATES = [1, 5, 10]`。试跑使用独立结果容器，不与正式样本合并，也不自动修改生产参数。

只有全部监测变量通过诊断的候选才会参与推荐。评分采用各变量中最低的 bulk ESS/秒，时间包括 warm-up。若均不合格，输出空推荐，不能仅凭速度决定更新数量。正式比较还应考虑运行时间波动，必要时重复试跑。

## 2. 初态与随机种子

默认采用 3 种初态与 4 个 seed 的笛卡尔积，共 12 条链：

- `compact_aligned`：紧凑 cluster，body 和 SL 初始 tilt 均为 0°。
- `compact_tilted`：整体共同旋转 20°，颗粒相对位置不变。
- `expanded_tilted`：初始中心坐标扩大 8%，整体共同旋转 40°。
- `MULTI_SEEDS = [1, 2, 3, 4]`。

第三种初态改变了颗粒间距，因此不只是整体姿态不同。这些初态仍属于当前 27 颗粒 cluster 类型，不能代表所有组装拓扑。所有链的 Hamiltonian 相同，膨胀仅作用于初始坐标，不改变 `A_NM` 或 vdW 截断。

初始化与采样分别使用 NumPy `SeedSequence([seed, structure_id, stream_id])`。其中 structure_id 是结构在 `INITIAL_STRUCTURES` 中的位置，stream_id 分别为 0 和 1。实际序列参数保存在 chain CSV 中，可以复现。

多链默认关闭。开启 `RUN_MULTI_SEED_SCAN` 才会运行完整实验。不同结构一起比较时，两类机械自由度必须都允许运动，否则不同链可能对应不同的条件分布。

## 3. Metropolis–Hastings

目标权重为 `exp(-U/kBT)`。接受概率是：

`min(1, exp(-delta_U/kBT + log_q_reverse - log_q_forward))`

当前随机平移、随机轴小角旋转和全局旋转 proposal 都是对称的，因此 log proposal ratio 为 0。这时 MH 就退化为原来的 Metropolis 接受率，使用 MH 名称不意味着目标分布发生变化。若以后加入不对称 proposal，必须传入正确的反向/正向密度比。

代码在 log 空间判断接受，避免直接计算巨大指数。当前没有在正式采样期间自适应修改步长。连续被拒绝的状态仍按 cycle 保存，不能删除拒绝样本。

## 4. 三类统计量的含义

### Rank-normalized split R-hat

正确名称是 rank-normalized，而非 rand-normalized。它先把每条链分成前后两半，再对合并样本做秩正态化，比较链内与链间差异。ArviZ 的 rank 方法还包含 folded 检查，可以发现分布宽度的差异。

- 接近 1：目前没有发现明显的链间不一致。
- 明显大于 1：不同链或前后半链仍有差异，可能尚未混合。
- 默认门槛：小于 1.01。

该门槛不是收敛证明。所有链困在同一个未被识别的局部构型，也可能得到接近 1 的结果。

### ESS

ESS 是有效样本数。MC 样本彼此相关，1500 个 cycle 不一定有 1500 个独立样本的信息量。

对稳定链，直观近似是：

`ESS = 样本总数 / (1 + 2 * 自相关系数之和)`

本项目分别输出：

- bulk ESS：分布主体的采样信息量，采用秩转换。
- tail ESS：5% 与 95% 分位附近的采样信息量，报告两者中较小值。
- mean ESS：针对原始变量均值的有效样本数，用于解释均值精度。

三者不能互相替代。默认要求至少 400，但真正需要多少取决于目标精度。ESS 可以因负自相关略大于原始样本数，不应机械截断。

### MCSE

MCSE 是 Monte Carlo standard error，即由有限采样导致的估计误差。本项目计算均值的 MCSE，其直观形式为：

`MCSE(mean) ≈ 样本标准差 / sqrt(mean ESS)`

例如平均 tilt 为 12°，MCSE 为 0.4°，表示这个均值估计的采样误差约为 0.4°。不表示所有颗粒都在 12±0.4° 之间，也不包括物理参数、几何近似或未收敛造成的误差。

默认要求 body tilt 和 SL tilt 的 MCSE 不超过 0.5°。这是可修改的精度目标，不是已达到的实验结果。

## 5. 数据如何组织

每个变量保留 `(chain, draw)` 数组。每条链分别去掉前 `MULTI_EQUIL` 个样本，不进行 thinning。不能把不同链首尾相接后当作一条链，也不能只拿每条链的均值计算 R-hat。

诊断分别在以下两层计算：

1. 所有初态和 seed 合并。
2. 每一种初态内部的各个 seed。

第一层用于检查初态依赖，第二层辅助区分初态差异和同一初态下的采样波动。结果保存在 diagnostics CSV 的 scope 字段。

少于两条链、每链少于 100 个保留样本、恒定链或非有限值会产生明确标记，诊断值记为不可用。至少四条链才允许通过最终检查。恒定链不会被报告为“MCSE 为零，因此已经精确收敛”。

NPZ 保存所有 cycle 的轨迹和最终状态。JSON 保存实际运行时的参数。最后一个 notebook cell 根据内存中的 `multi_bundle` 重新导出，采用原始运行设置，不会因之后修改前部参数而伪装成另一个实验。

## 6. 修正的观测量与控制行为

- local beta 现在是磁矩与自身 body [111] 的有向夹角，范围为 0–180°。旧实现实际测量 body–field 夹角。
- body–SL mismatch 现在是 coherent body tilt 减去 SL PCA tilt，可能为负。
- 两条轴的实际夹角作为独立字段保留，不能由 tilt 差值代替。
- 有符号磁化与绝对磁场对齐量同时保存。
- 位置或姿态开关关闭时，全局 move 也遵守限制。
- 中心位置仍固定，但中心颗粒的本体方向现在可以局部旋转。
- `include_vdw=False` 时，初始化与所有能量差使用同一开关。

这些修改使当前输出与定义一致，但不改变相互作用势。旧报告中的字段不应直接与新字段混用。

## 7. 本次验证与后续实验

回归测试覆盖磁能量局部差值与完整差值一致性、每个 sweep 的更新次数、开关、种子复现、不对称 proposal 的 MH 分布、短链状态、独立正态样本、偏移链、自相关链、报告导出和 notebook 格式。

短测试只证明实现的这些性质，不能证明 27 颗粒体系已经平衡，也不能决定最佳磁 sweep 数量。下一步正式实验应先检查初态依赖，再比较 ESS/秒和目标 MCSE。达到给定 cycle 数量本身不是停止依据。

当前模型仍没有显式有限容器。上述诊断不替代后续的 vdW、磁相互作用、构型熵与体系边界验证。MC cycle 始终不是物理时间。

方法来源：

- [Vehtari et al., 2021](https://arxiv.org/abs/1903.08008)
- [ArviZ 0.22 R-hat](https://python.arviz.org/en/v0.22.0/api/generated/arviz.rhat.html)
- [ArviZ 0.22 ESS](https://python.arviz.org/en/v0.22.0/api/generated/arviz.ess.html)
- [ArviZ 0.22 MCSE](https://python.arviz.org/en/v0.22.0/api/generated/arviz.mcse.html)
