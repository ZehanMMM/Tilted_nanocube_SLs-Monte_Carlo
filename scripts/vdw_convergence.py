"""Voxel-refinement study: the V10 sharp-cube sum converges, and to what.

    python scripts/vdw_convergence.py [outputs/vdw_convergence]

Three levels of evidence, each with n^3 voxels per cube and no d^2 floor:
  A  co-oriented face pairs at 1.8 nm (V9 steric gap) and 3.0 nm (V10 gap),
     via V10's displacement-multiplicity sum up to n = 128;
  B  single rotated pairs taken from V9.03 and near-contact geometries,
     brute-force n^6 sums up to n = 32;
  C  the total vdW of the whole 27-cube starting cluster, all 351 pairs.
Each is compared with the fast evaluator used in the MC.
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from v904 import geometry  # noqa: E402
from v904.model import Hamiltonian, INITIAL_STRUCTURES, initial_state  # noqa: E402
from v904.vdw import (HamakerEvaluator, parallel_voxel_energy, richardson,  # noqa: E402
                      uniform_voxel_energy)

OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "outputs" / "vdw_convergence"
OUT.mkdir(parents=True, exist_ok=True)
H = Hamiltonian()
KT, L = H.kT, H.L
FAST = H.hamaker
I3 = np.eye(3)


def write_csv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    rows, summary = [], {}

    # ---------------------------------------------------------------- A
    ns_a = [4, 8, 12, 16, 24, 32, 48, 64, 96, 128]
    for gap_nm in (1.8, 3.0):
        r = np.array([L + gap_nm * 1e-9, 0, 0])
        fast = FAST.pair_energy(np.zeros(3), I3, r, I3)
        vals = [parallel_voxel_energy(r, n) for n in ns_a]
        for k, (n, v) in enumerate(zip(ns_a, vals)):
            rows.append(dict(case=f"A_face_{gap_nm}nm", n=n, energy_kBT=v / KT,
                             rel_err=v / fast - 1,
                             richardson_rel_err=(richardson(vals[:k + 1], ns_a[:k + 1]) / fast - 1)
                             if k else np.nan, fast_kBT=fast / KT))
        summary[f"face_{gap_nm}nm"] = dict(fast_J=fast, fast_kBT=fast / KT,
                                          n4_J=vals[0], ratio_fast_over_n4=fast / vals[0])

    # ---------------------------------------------------------------- B
    s0 = initial_state(**INITIAL_STRUCTURES["compact_aligned"], dipole_seed=1)
    i, j = np.triu_indices(27, 1)
    e_all = FAST.pair_energies(s0["pos"][i], s0["Q"][i], s0["pos"][j], s0["Q"][j])
    k = int(np.argmin(e_all))
    rng = np.random.default_rng(11)
    from scipy.spatial.transform import Rotation
    cases = {"B_V903_strongest_pair": (s0["pos"][i[k]], s0["Q"][i[k]], s0["pos"][j[k]], s0["Q"][j[k]])}
    for label, lo, hi in (("B_tilted_core_2nm", 1.5, 2.5), ("B_tilted_core_0.9nm", 0.7, 1.1)):
        while True:
            Qa = Rotation.from_rotvec(np.radians(rng.normal(0, 6, 3))).as_matrix()
            Qb = Rotation.from_rotvec(np.radians(rng.normal(0, 6, 3))).as_matrix()
            r = np.array([L + rng.uniform(1, 5) * 1e-9, *rng.uniform(-4e-9, 4e-9, 2)])
            d = geometry.core_distance(np.zeros(3), Qa, r, Qb, L) * 1e9
            gap = geometry.support_gap(r, Qa, Qb, L, H.r_round)[0] * 1e9
            if lo <= d <= hi and gap >= 1.8:
                cases[label] = (np.zeros(3), Qa, r, Qb)
                break
    ns_b = [4, 8, 12, 16, 24, 32]
    for label, (pa, Qa, pb, Qb) in cases.items():
        fast = FAST.pair_energy(pa, Qa, pb, Qb)
        vals = []
        for n in ns_b:
            t = time.time()
            vals.append(uniform_voxel_energy(pa, Qa, pb, Qb, n))
            rows.append(dict(case=label, n=n, energy_kBT=vals[-1] / KT, rel_err=vals[-1] / fast - 1,
                             richardson_rel_err=(richardson(vals, ns_b[:len(vals)]) / fast - 1)
                             if len(vals) > 1 else np.nan, fast_kBT=fast / KT))
            print(label, n, f"{time.time() - t:.1f}s", flush=True)
        summary[label] = dict(fast_kBT=fast / KT,
                              core_distance_nm=geometry.core_distance(pa, Qa, pb, Qb, L) * 1e9,
                              inherited_n4_floor_kBT=uniform_voxel_energy(
                                  pa, Qa, pb, Qb, 4, d2_floor_m2=1e-19) / KT)

    # ---------------------------------------------------------------- C
    ns_c = [4, 8, 12, 16]
    for name in ("compact_aligned", "expanded_tilted"):
        s = initial_state(**INITIAL_STRUCTURES[name], dipole_seed=1)
        fast = float(FAST.pair_energies(s["pos"][i], s["Q"][i], s["pos"][j], s["Q"][j]).sum())
        vals = []
        for n in ns_c:
            t = time.time()
            vals.append(sum(uniform_voxel_energy(s["pos"][a], s["Q"][a], s["pos"][b], s["Q"][b], n)
                            for a, b in zip(i, j)))
            rows.append(dict(case=f"C_cluster_{name}", n=n, energy_kBT=vals[-1] / KT,
                             rel_err=vals[-1] / fast - 1,
                             richardson_rel_err=(richardson(vals, ns_c[:len(vals)]) / fast - 1)
                             if len(vals) > 1 else np.nan, fast_kBT=fast / KT))
            print(name, n, f"{time.time() - t:.1f}s", flush=True)
        inherited = Hamiltonian(__import__("v904.model", fromlist=["PhysicalParameters"])
                                .PhysicalParameters(vdw_model="inherited"))
        summary[f"cluster_{name}"] = dict(
            fast_kBT=fast / KT,
            inherited_V903_kBT=inherited.energy_terms(s["pos"], s["Q"], s["mu"])["VdW"] / KT)

    write_csv(OUT / "vdw_voxel_convergence.csv", rows)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    # ---------------------------------------------------------------- figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    colors = ["#2a6f97", "#c44e52", "#55a868", "#8172b2", "#dd8452", "#64b5cd", "#937860"]
    for c, case in enumerate(dict.fromkeys(r["case"] for r in rows)):
        sub = [r for r in rows if r["case"] == case]
        n = np.array([r["n"] for r in sub])
        err = np.abs([r["rel_err"] for r in sub])
        axes[0].loglog(n, err, "o-", color=colors[c % len(colors)], label=case, lw=1.6, ms=4)
        rich = np.abs([r["richardson_rel_err"] for r in sub])
        axes[1].loglog(n[1:], rich[1:], "s--", color=colors[c % len(colors)], lw=1.2, ms=4)
    ref = np.array([8, 128])
    axes[0].loglog(ref, 0.3 * (ref / 8.0) ** -2, "k:", lw=1, label=r"$\propto n^{-2}$ (midpoint rule)")
    axes[0].set(xlabel="voxels per cube edge n", ylabel="|voxel sum / converged - 1|",
                title="n^3 voxel sum (V10 method, no floor)")
    axes[1].set(xlabel="finer n of the pair", ylabel="|Richardson / converged - 1|",
                title="h^2 Richardson extrapolation")
    axes[0].legend(fontsize=7)
    for ax in axes:
        ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "vdw_voxel_convergence.png", dpi=160)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
