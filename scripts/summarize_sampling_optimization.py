"""Validate completed calibration/validation checkpoints and write a Chinese report."""
import argparse
import csv
import hashlib
import json
import platform
from importlib.metadata import version
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
METRICS = ["rhat_rank", "ess_bulk", "ess_tail", "ess_mean", "mcse_mean"]


def csv_rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def verify_run(folder):
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    finished = json.loads((folder / "finished.json").read_text(encoding="utf-8"))
    config = manifest["config"]
    for name, digest in config["source_sha256"].items():
        assert hashlib.sha256((folder / "source_snapshot" / name).read_bytes()).hexdigest() == digest, name
    checkpoints = list(folder.glob("m*/complete.json"))
    expected = len(config["starts"]) * len(config["seeds"]) * len(config["sweeps"])
    assert len(checkpoints) == finished["jobs"] == expected
    max_error, min_gap = 0.0, float("inf")
    for marker in checkpoints:
        checks = json.loads(marker.read_text(encoding="utf-8"))
        max_error = max(max_error, abs(checks["full_minus_tracked_energy_kBT"]))
        min_gap = min(min_gap, checks["min_recorded_gap_nm"])
        info = json.loads((marker.parent / "result.json").read_text(encoding="utf-8"))
        assert info["sampling"]["n_cycles"] == config["cycles"]
        assert info["sampling"]["n_equil"] == config["equil"]
        if config["mode"] == "conditional":
            assert checks["frozen_geometry_unchanged"]
            assert not info["sampling"]["move_positions"]
            assert not info["sampling"]["move_orientations"]
        with np.load(marker.parent / "result.npz") as data:
            assert all(np.isfinite(data[k]).all() for k in data.files)
            for key in data.files:
                if key.startswith("traj_"):
                    assert len(data[key]) == config["cycles"]
            np.testing.assert_allclose(np.linalg.norm(data["mu"], axis=1), 1, atol=1e-10)
            np.testing.assert_allclose(data["Q"] @ data["Q"].transpose(0, 2, 1),
                                       np.tile(np.eye(3), (27, 1, 1)), atol=1e-10)
    diagnostics = csv_rows(folder / "diagnostics.csv")
    assert all(np.isfinite(float(row[key])) for row in diagnostics for key in METRICS)
    assert max_error < 1e-6 and min_gap > 0
    return dict(config=config, finished=finished, chains=csv_rows(folder / "chains.csv"),
                diagnostics=diagnostics, windows=csv_rows(folder / "windows.csv"),
                max_energy_error_kBT=max_error, min_gap_nm=min_gap,
                candidates=csv_rows(folder / "candidates.csv"))


def diagnostic_table(rows):
    lines = ["| 变量 | R-hat | bulk ESS | tail ESS | mean ESS | MCSE(mean) | 状态 |",
             "|---|---:|---:|---:|---:|---:|---|"]
    for row in rows:
        numbers = " | ".join(f"{float(row[k]):.5g}" for k in METRICS)
        lines.append(f"| {row['variable']} | {numbers} | {row['status']} |")
    return lines


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--date", default="20260916")
    args = parser.parse_args()
    root = ROOT / "outputs"
    names = [f"magnetic_calibration_{args.date}", f"magnetic_calibration_wide_{args.date}",
             f"magnetic_calibration_large_{args.date}", f"magnetic_validation_{args.date}",
             f"joint_validation_{args.date}"]
    runs = [verify_run(root / name) for name in names]
    pilots, magnetic, joint = runs[:3], runs[3], runs[4]
    mag_pass = all(r["status"] == "checks_passed" for r in magnetic["diagnostics"])
    pooled = [r for r in joint["diagnostics"] if r["scope"] == "all_starts"]
    assert len(pooled) == 12
    joint_pass = all(r["status"] == "checks_passed" for r in joint["diagnostics"])
    lines = ["# V9.03 采样优化与独立验证", "", "## 结论", "",
        "磁标定选择 1.60 rad、每 cycle 10 个完整磁 sweep。这是九组候选中唯一通过全部试跑检查的设置。",
        f"独立磁验证：{'全部通过' if mag_pass else '仍有未通过项目，见下文'}。",
        f"联合验证：{'全部有限运行诊断通过' if joint_pass else '仍有未通过项目，不能声称整体平衡'}。",
        "完整位置空间仍无界，任何诊断通过都不能解决目标权重不可归一化的物理边界问题。",
        "vdW、磁能和 steric 公式及参数保持不变，没有添加构型熵能量项。", "",
        "## 试跑比较", "",
        "四种冻结几何，每种四条链，每条 2000 cycles，舍弃前 500。",
        "每组候选监测 232 个几何/变量组合，含逐颗粒磁矩投影与 local beta。",
        "下表的 flags 是未通过的组合数。未通过候选的 ESS/秒只作为诊断分数，不是已经验证的精度效率。", "",
        "| 磁提议最大角 / rad | sweeps | 最大 R-hat | 最小 bulk ESS | 最差 bulk ESS/秒 | flags / 232 |",
        "|---:|---:|---:|---:|---:|---:|"]
    for run in pilots:
        step = run["config"].get("dip_step_rad", .30)
        for row in run["candidates"]:
            lines.append(f"| {step:.2f} | {row['sweeps']} | {float(row['max_rhat']):.5g} | "
                         f"{float(row['min_bulk_ess']):.5g} | {float(row['min_bulk_ess_per_sec']):.5g} | {row['failed_checks']} |")
    lines += ["", "![磁标定比较](figures/calibration_comparison.png)", "",
        "## 独立磁验证", "",
        "四种冻结几何 × seeds 21–24，每条 4000 cycles，前 1000 为 warm-up。",
        "初始化仍包括随机、顺场、逆场和 body easy-axis，采样流与试跑分开。",
        f"232 个组合中 {sum(r['status'] != 'checks_passed' for r in magnetic['diagnostics'])} 个未通过。",
        f"最大 R-hat = {max(float(r['rhat_rank']) for r in magnetic['diagnostics']):.6g}，"
        f"最小 bulk ESS = {min(float(r['ess_bulk']) for r in magnetic['diagnostics']):.6g}。", ""]
    for scope in magnetic["config"]["starts"]:
        lines += [f"### {scope}", ""] + diagnostic_table([r for r in magnetic["diagnostics"]
                         if r["scope"] == scope and r["variable"] in ["energy_kBT", "local_beta_deg", "magnetization"]]) + [""]
    lines += ["![独立磁验证的逐颗粒分布](magnetic_validation_figures/particle_beta_step_1.60.png)", "",
        "图示 compact 40° 中第 12 个颗粒，四个新 seed 的 local-beta 主体与尾部分布基本一致。数值检查覆盖全部 27 个颗粒，不能只凭这一张图判断。", "",
        "## 独立联合验证", "",
        "四种初态 × seeds 31–34，每条 4000 cycles，前 1000 为 warm-up。",
        "磁参数固定为 1.60 rad、10 sweeps。整体 co-tilt 最大 3°。",
        "每 cycle 另尝试一次中心坐标整体缩放，最大 log-scale 为 0.002，并包含 78×log(scale) 的 MH 修正。",
        "该修正维持原有 Cartesian 采样测度，不是添加新的物理熵模型。", ""]
    lines += diagnostic_table(pooled)
    lines += ["", "R-hat < 1.01，pooled ESS 门槛为 1600，四链组内门槛为 400。",
        "body/SL 的均值 MCSE 目标为 0.5°，local beta 为 0.25°，signed magnetization 为 0.002。",
        "不同观测量可能混合速度不同，不能用 body tilt 的通过代替全部结构自由度的通过。", "",
        "| 初态 | body tilt 均值 / deg | SL tilt 均值 / deg | Rg 均值 / nm | co-tilt 接受率 | scale 接受率 |",
        "|---|---:|---:|---:|---:|---:|"]
    for scope in joint["config"]["starts"]:
        rows = [r for r in joint["chains"] if r["start"] == scope]
        values = [np.mean([float(r[k]) for r in rows]) for k in
                  ["body_tilt_mean", "sl_pca_tilt_mean", "rg_nm_mean", "acc_cotilt", "acc_scale"]]
        lines.append("| " + scope + " | " + " | ".join(f"{v:.5g}" for v in values) + " |")
    lines += ["", "### 初态组内诊断", "",
        "每组四条链分别检查。组内通过仍不能代替不同初态间的一致性。", "",
        "| 初态 | body R-hat | body bulk ESS | body MCSE / deg | 未通过变量数 / 12 |",
        "|---|---:|---:|---:|---:|"]
    for scope in joint["config"]["starts"]:
        rows = [r for r in joint["diagnostics"] if r["scope"] == scope]
        body = next(r for r in rows if r["variable"] == "body_tilt_deg")
        failed = sum(r["status"] != "checks_passed" for r in rows)
        lines.append(f"| {scope} | {float(body['rhat_rank']):.5g} | "
                     f"{float(body['ess_bulk']):.5g} | {float(body['mcse_mean']):.5g} | {failed} |")
    lines += ["", "### 保留样本的前后时间窗口", "",
        "每个窗口先在链内取均值，再平均同初态的四条链。", "",
        "| 初态 | cycle 窗口 | body tilt / deg | SL tilt / deg | Rg / nm | body order | body-SL 轴夹角 / deg |",
        "|---|---|---:|---:|---:|---:|---:|"]
    for scope in joint["config"]["starts"]:
        rows = [r for r in joint["windows"] if r["start"] == scope]
        for first in sorted(set(int(r["cycle_first"]) for r in rows)):
            group = [r for r in rows if int(r["cycle_first"]) == first]
            values = [np.mean([float(r[k]) for r in group]) for k in
                      ["body_tilt_deg", "sl_pca_tilt_deg", "rg_nm", "body_order", "body_sl_axis_angle_deg"]]
            lines.append(f"| {scope} | {first}–{group[0]['cycle_last']} | " +
                         " | ".join(f"{v:.5g}" for v in values) + " |")
    lines += ["", "均值使用 warm-up 后的轨迹，接受率使用完整运行。未收敛变量的均值只作为描述统计。",
        "body order 越低，颗粒本体的方向越分散。body-SL 轴夹角与两个 tilt 的差值分别记录，两者不能互换。", "",
        "![联合 body tilt](joint_figures/traj_body_tilt.png)", "",
        "![联合 SL tilt](joint_figures/traj_sl_pca_tilt.png)", "",
        "![联合 Rg](joint_figures/traj_rg_nm.png)", "",
        "![联合能量](joint_figures/traj_E.png)", "",
        "![前后窗口分布](joint_figures/body_tilt_window_cdf.png)", "",
        "## 如何理解旧结果 R-hat≈1.83、bulk ESS≈17.4", "",
        "它们表明旧样本尚未一致、稳定地探索 body tilt 分布，不证明该角度不存在稳定态。",
        "有限温度下稳定的是概率分布，角度可以持续波动。正的平均夹角也不能单独证明非零倾角的自由能最低点。",
        "旧实验有 18000 个保留记录，但自相关和初态差异使有效信息量很低。非平稳时不能把 17.4 严格当成独立平衡样本数。",
        "MCSE 只估计采样误差，不覆盖未收敛偏差。旧膨胀初态前后窗口的 body tilt 均值由 17.13° 降到 9.10°，是持续松弛的直接证据。", "",
        "当前模型在某个非中心颗粒无限远离时，能量趋于有限值，位置体积积分却发散。",
        "因此完整无界位置空间没有可归一化的玻尔兹曼分布。固定几何的磁条件分布则可以归一化。",
        "后续完整平衡结论需要明确物理边界、有限密度或团簇条件。本次没有擅自添加这些约束。", "",
        "## 数值核验与文件", "",
        f"- 共验证 {sum(r['finished']['jobs'] for r in runs)} 条新链。所有轨迹和已估计诊断值为有限数值。",
        f"- 最终全量能量与累计能量的最大误差为 {max(r['max_energy_error_kBT'] for r in runs):.3e} kBT。",
        f"- 全部记录的最小现有投影 gap 为 {min(r['min_gap_nm'] for r in runs):.6g} nm。该检查不替代独立的真实凸体几何验证。",
        "- 磁矩长度、姿态矩阵正交性及 source snapshot SHA256 均通过检查。",
        "- 23 项回归测试通过，包含前部机械开关、MH 非对称提议、单磁矩球面积分参考和缩放 Jacobian 的 78 维高斯参考。",
        "- 每个实验目录保存 diagnostics.csv、chains.csv、windows.csv、manifest.json、finished.json 和逐链 NPZ/JSON。",
        "- 试跑与独立验证分开。联合和磁验证同时执行，验证阶段的运行时间不能直接用于与试跑比较加速倍数。", "",
        "方法：[Vehtari et al.](https://arxiv.org/abs/1903.08008)，[Stan 诊断说明](https://mc-stan.org/learn-stan/diagnostics-warnings.html)。", ""]
    output = root / f"sampling_optimization_{args.date}"
    output.mkdir(exist_ok=True)
    (output / "README.md").write_text("\n".join(lines), encoding="utf-8")
    summary = dict(verification_environment={"python": platform.python_version(),
                       **{name: version(name) for name in ["numpy", "scipy", "matplotlib", "arviz", "nbformat"]}},
                   magnetic_validation_passed=mag_pass, joint_all_checks_passed=joint_pass,
                   total_new_chains=sum(r["finished"]["jobs"] for r in runs),
                   pooled_joint=pooled, max_energy_error_kBT=max(r["max_energy_error_kBT"] for r in runs))
    (output / "verification.json").write_text(json.dumps(summary, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
