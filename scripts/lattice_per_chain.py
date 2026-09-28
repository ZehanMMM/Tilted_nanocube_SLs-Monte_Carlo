"""Every rigid-lattice Monte Carlo chain, one row each: the interval it stays in.

    python scripts/lattice_per_chain.py OUT_DIR

Angles are all relative to the external field B (lab x) unless stated:

  SL axis tilt     angle between the SL three-fold axis G x and B, 0-90 deg
  SL azimuth       direction in which that axis leans, around B (from G)
  body tilt        angle between the mean cube [111] and B
  cube twist       signed rotation of each cube about the lattice axis,
                   relative to its reference orientation in the lattice
  cube off-axis    angle between each cube's [111] and the lattice axis

Sampling chains: the first run (cycles 501-2000, after its 500-cycle warm-up)
and its continuation (cycles 2001-3500) are listed separately and together.
Nothing here depends on convergence: min / max / 2.5-97.5 % quantiles are
the bounds each chain actually visited.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from v904.model import LATTICE_BODY_REF  # noqa: E402

RUN1 = ROOT / "outputs" / "lattice_production"
RUN2 = ROOT / "outputs" / "lattice_production_cont"
ANNEAL = ROOT / "outputs" / "lattice_anneal"
WARM = 500
COLORS = {"lat_tilt00": "#2a6f97", "lat_tilt20": "#55a868", "lat_tilt40": "#c44e52",
          "lat_tilt60": "#dd8452"}


def interval(x):
    x = np.asarray(x, float)
    return dict(mean=float(x.mean()), sd=float(x.std()), min=float(x.min()),
                q025=float(np.quantile(x, 0.025)), q975=float(np.quantile(x, 0.975)),
                max=float(x.max()))


def snapshot_angles(path, first_cycle=0):
    """Per snapshot: SL azimuth, mean/min/max cube twist, mean cube off-axis angle."""
    out = dict(cycle=[], azimuth=[], tilt=[], twist_mean=[], twist_all=[], offaxis=[])
    with np.load(path) as d:
        for k in range(len(d["cycle"])):
            G, Q = d["G"][k], d["Q"][k]
            axis = G[:, 0]
            if axis[0] < 0:
                axis = -axis
            out["cycle"].append(int(d["cycle"][k]) + first_cycle)
            out["tilt"].append(float(np.degrees(np.arccos(np.clip(abs(axis[0]), 0, 1)))))
            out["azimuth"].append(float(np.degrees(np.arctan2(axis[2], axis[1]))))
            rel = np.einsum("ba,nbc,dc->nad", G, Q, LATTICE_BODY_REF)
            rv = Rotation.from_matrix(rel).as_rotvec()
            tw = np.degrees(rv[:, 0])
            out["twist_mean"].append(float(tw.mean()))
            out["twist_all"].append(tw)
            out["offaxis"].append(float(np.degrees(np.linalg.norm(rv[:, 1:], axis=1)).mean()))
    return out


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "outputs" / "lattice_per_chain"
    out.mkdir(parents=True, exist_ok=True)
    rows, traces, snaps = [], {}, {}
    for f1 in sorted(RUN1.iterdir()):
        if not (f1 / "complete.json").exists():
            continue
        info = json.loads((f1 / "chain.json").read_text(encoding="utf-8"))
        start, seed = info["start"], info["seed"]
        f2 = RUN2 / f"{start}_s{seed + 20}"
        with np.load(f1 / "trajectory.npz") as a, np.load(f2 / "trajectory.npz") as b:
            t = {k: np.r_[a[k], b[k]] for k in ("sl_axis_tilt_deg", "body_tilt_deg",
                                                "magnetization", "energy_kBT")}
        n1 = 2000
        s1 = snapshot_angles(f1 / "snapshots.npz")
        s2 = snapshot_angles(f2 / "snapshots.npz", first_cycle=n1)
        s = {k: s1[k] + s2[k] for k in s1}
        keep = np.array(s["cycle"]) > WARM
        tw_all = np.concatenate([x for x, k in zip(s["twist_all"], keep) if k])
        name = f"{start}_s{seed}"
        traces[name] = (start, t["sl_axis_tilt_deg"])
        snaps[name] = (start, s, keep)
        for label, sl in (("run 1 (501-2000)", slice(WARM, n1)),
                          ("continuation (2001-3500)", slice(n1, None)),
                          ("both (501-3500)", slice(WARM, None))):
            sk = keep & ((np.array(s["cycle"]) <= n1) if label.startswith("run 1") else
                         (np.array(s["cycle"]) > n1) if label.startswith("cont") else True)
            row = dict(chain=name, start_tilt_deg=float(start[-2:]), segment=label,
                       n_samples=int(len(t["sl_axis_tilt_deg"][sl])))
            for key, lab in (("sl_axis_tilt_deg", "sl_tilt"), ("body_tilt_deg", "body_tilt"),
                             ("magnetization", "m_x"), ("energy_kBT", "U_kBT")):
                for stat, v in interval(t[key][sl]).items():
                    row[f"{lab}_{stat}"] = v
            tm = np.array(s["twist_mean"])[sk]
            oa = np.array(s["offaxis"])[sk]
            tw = np.concatenate([x for x, k in zip(s["twist_all"], sk) if k])
            row.update(twist_mean_of_cubes_mean=float(tm.mean()), twist_mean_of_cubes_min=float(tm.min()),
                       twist_mean_of_cubes_max=float(tm.max()),
                       twist_single_cube_min=float(tw.min()), twist_single_cube_max=float(tw.max()),
                       offaxis_mean=float(oa.mean()), offaxis_max=float(oa.max()))
            rows.append(row)
        _ = tw_all
    with open(out / "per_chain_sampling.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # all chains together, after warm-up
    tilt_all = np.concatenate([tr[WARM:] for _, tr in traces.values()])
    tw_pool = np.concatenate([np.concatenate([x for x, k in zip(s["twist_all"], keep) if k])
                              for _, s, keep in snaps.values()])
    twm_pool = np.concatenate([np.array(s["twist_mean"])[keep] for _, s, keep in snaps.values()])
    overall = dict(sl_tilt=interval(tilt_all), cube_twist_single=interval(tw_pool),
                   cube_twist_mean=interval(twm_pool))

    # annealing restarts
    arows = []
    for fa in sorted(ANNEAL.iterdir()):
        if not (fa / "complete.json").exists():
            continue
        info = json.loads((fa / "chain.json").read_text(encoding="utf-8"))
        with np.load(fa / "best_state.npz") as d:
            G, Q = d["G"], d["Q"]
        axis = G[:, 0] * np.sign(G[0, 0])
        rel = np.einsum("ba,nbc,dc->nad", G, Q, LATTICE_BODY_REF)
        rv = Rotation.from_matrix(rel).as_rotvec()
        tw = np.degrees(rv[:, 0])
        with np.load(fa / "trajectory.npz") as d:
            late = d["sl_axis_tilt_deg"][-300:]           # last 300 cycles, kT <= 0.03
        arows.append(dict(restart=fa.name, start_tilt_deg=float(info["start"][-2:]),
                          best_energy_kBT=info["best_energy_kBT"],
                          sl_tilt_best_deg=float(np.degrees(np.arccos(min(1.0, abs(axis[0]))))),
                          sl_tilt_last300_min=float(late.min()), sl_tilt_last300_max=float(late.max()),
                          twist_mean_deg=float(tw.mean()), twist_min_deg=float(tw.min()),
                          twist_max_deg=float(tw.max()), n_cubes_twisted_gt20=int(np.sum(tw > 20))))
    with open(out / "per_restart_annealing.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(arows[0]))
        w.writeheader()
        w.writerows(arows)
    (out / "overall_intervals.json").write_text(json.dumps(overall, indent=2), encoding="utf-8")

    # --------------------------------------------------------------- report
    both = [r for r in rows if r["segment"].startswith("both")]
    L = ["# Every rigid-lattice Monte Carlo chain (298 K) and every annealing restart\n",
         "All angles in degrees.  SL tilt = angle between the SL three-fold axis and the field B.",
         "Sampling: cycles 501-3500 of each chain (run 1 after its warm-up + continuation).\n",
         "| chain | start tilt | SL tilt mean ± sd | SL tilt 2.5-97.5 % | SL tilt min-max | body tilt 2.5-97.5 % "
         "| mean cube twist (range over snapshots) | single-cube twist min-max | cube off-axis mean | m_x 2.5-97.5 % |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in both:
        L.append(f"| {r['chain']} | {r['start_tilt_deg']:.0f} | {r['sl_tilt_mean']:.2f} ± {r['sl_tilt_sd']:.2f} "
                 f"| {r['sl_tilt_q025']:.2f}-{r['sl_tilt_q975']:.2f} | {r['sl_tilt_min']:.2f}-{r['sl_tilt_max']:.2f} "
                 f"| {r['body_tilt_q025']:.2f}-{r['body_tilt_q975']:.2f} "
                 f"| {r['twist_mean_of_cubes_mean']:+.1f} ({r['twist_mean_of_cubes_min']:+.1f} to {r['twist_mean_of_cubes_max']:+.1f}) "
                 f"| {r['twist_single_cube_min']:+.1f} to {r['twist_single_cube_max']:+.1f} "
                 f"| {r['offaxis_mean']:.1f} | {r['m_x_q025']:.3f}-{r['m_x_q975']:.3f} |")
    o = overall
    L += ["", "## All 16 chains together (cycles 501-3500, 48 000 samples of the SL tilt)", "",
          f"- SL tilt: mean {o['sl_tilt']['mean']:.2f}, 2.5-97.5 % {o['sl_tilt']['q025']:.2f}-{o['sl_tilt']['q975']:.2f}, "
          f"min-max {o['sl_tilt']['min']:.2f}-{o['sl_tilt']['max']:.2f}",
          f"- single-cube twist: 2.5-97.5 % {o['cube_twist_single']['q025']:+.1f} to {o['cube_twist_single']['q975']:+.1f}, "
          f"min-max {o['cube_twist_single']['min']:+.1f} to {o['cube_twist_single']['max']:+.1f}",
          f"- mean twist of the 27 cubes: {o['cube_twist_mean']['min']:+.1f} to {o['cube_twist_mean']['max']:+.1f}",
          "", "## Annealing restarts (minimum of the total energy)", "",
          "| restart | start tilt | best U / kBT | SL tilt at best | SL tilt, last 300 cycles | mean cube twist | cube twist min-max | cubes twisted > 20 deg |",
          "|---|---|---|---|---|---|---|---|"]
    for a in arows:
        L.append(f"| {a['restart']} | {a['start_tilt_deg']:.0f} | {a['best_energy_kBT']:.2f} | {a['sl_tilt_best_deg']:.2f} "
                 f"| {a['sl_tilt_last300_min']:.2f}-{a['sl_tilt_last300_max']:.2f} | {a['twist_mean_deg']:+.1f} "
                 f"| {a['twist_min_deg']:+.1f} to {a['twist_max_deg']:+.1f} | {a['n_cubes_twisted_gt20']} |")
    (out / "PER_CHAIN_RESULTS.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L))

    # --------------------------------------------------------------- figures
    fig, ax = plt.subplots(figsize=(12, 4.2))
    for name, (start, tr) in traces.items():
        ax.plot(np.arange(1, len(tr) + 1), tr, lw=0.35, color=COLORS[start], alpha=0.8)
    ax.axhspan(o["sl_tilt"]["q025"], o["sl_tilt"]["q975"], color="k", alpha=0.08,
               label="2.5-97.5 % of all chains after cycle 500")
    ax.axvline(WARM, color="k", ls="--", lw=0.8)
    ax.axvline(2000, color="k", lw=0.8)
    ax.set(xlabel="MC cycle (solid line: start of the continuation)", ylabel="SL axis tilt from B / deg",
           title="every chain: SL three-fold axis relative to the field")
    handles = [plt.Line2D([], [], color=c, label=f"start {s[-2:]} deg") for s, c in COLORS.items()]
    ax.legend(handles=handles + [plt.Rectangle((0, 0), 1, 1, color="k", alpha=0.08,
                                               label="2.5-97.5 % after warm-up")], fontsize=8, ncol=5)
    fig.tight_layout()
    fig.savefig(out / "sl_tilt_every_chain.png", dpi=140)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    names = list(traces)
    axes[0].boxplot([traces[n][1][WARM:] for n in names], vert=False, whis=(2.5, 97.5),
                    showfliers=False)
    axes[0].set_yticks(range(1, len(names) + 1), names, fontsize=7)
    axes[0].set(xlabel="SL axis tilt from B / deg (box: quartiles, whiskers: 2.5-97.5 %)",
                title="per-chain interval of the SL tilt")
    for name, (start, s, keep) in snaps.items():
        th = np.radians(np.array(s["tilt"])[keep])
        ph = np.radians(np.array(s["azimuth"])[keep])
        axes[1].scatter(np.degrees(th) * np.cos(ph), np.degrees(th) * np.sin(ph), s=4,
                        color=COLORS[start], alpha=0.6)
    for r in (5, 10, 15):
        c = plt.Circle((0, 0), r, fill=False, color="0.6", lw=0.6, ls=":")
        axes[1].add_patch(c)
    axes[1].set_aspect("equal")
    lim = 1.1 * max(15, o["sl_tilt"]["max"])
    axes[1].set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="tilt x cos(azimuth) / deg",
                ylabel="tilt x sin(azimuth) / deg",
                title="where the SL axis points (B at the centre; circles 5, 10, 15 deg)")
    axes[2].boxplot([np.concatenate([x for x, k in zip(snaps[n][1]["twist_all"], snaps[n][2]) if k])
                     for n in names], vert=False, whis=(2.5, 97.5), showfliers=False)
    axes[2].set_yticks(range(1, len(names) + 1), names, fontsize=7)
    axes[2].axvline(0, color="k", lw=0.6)
    axes[2].set(xlabel="single-cube twist about the lattice axis / deg",
                title="per-chain interval of the cube twist")
    fig.tight_layout()
    fig.savefig(out / "per_chain_intervals.png", dpi=140)


if __name__ == "__main__":
    main()
