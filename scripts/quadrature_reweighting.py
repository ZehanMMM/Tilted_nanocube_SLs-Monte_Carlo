"""Does the vdW quadrature error bias the sampled averages?  Reweighting check.

    python scripts/quadrature_reweighting.py RUN_DIR --equil 1000 [--config CONFIG.json]

With --config, the "production" quadrature is the one that run used (e.g.
the looser rigid-lattice setting); otherwise the package default.

MC samples pi_q ~ exp(-U_q / kT), where U_q is the energy with the production
quadrature.  The exact target has U_exact instead.  For any observable f

    E_exact[f] = E_q[f w] / E_q[w],     w = exp(-(U_exact - U_q) / kT),

so re-evaluating saved snapshots with a much tighter quadrature (a stand-in
for U_exact) measures the bias directly, with no new sampling.  If every
|U_tight - U_q| << kT, the weights are ~1 and the bias is negligible by
construction; the script reports the largest difference, the weight spread
and the reweighted shift of each observable.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from v904.model import Hamiltonian, PhysicalParameters  # noqa: E402
from v904.vdw import HamakerSettings  # noqa: E402

TIGHT = HamakerSettings(q_near=10, tol_kBT=1e-8, q_mid=16, q_far=10, far_pair_L=3.5,
                        volume_order=6)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run")
    p.add_argument("--equil", type=int, required=True)
    p.add_argument("--max-snapshots", type=int, default=12, help="per chain")
    p.add_argument("--config", default=None)
    a = p.parse_args()
    run = Path(a.run)
    if a.config:
        sys.path.insert(0, str(ROOT / "scripts"))
        from run_chains import build
        prod = build(json.loads(Path(a.config).read_text(encoding="utf-8")))[0]
    else:
        prod = Hamiltonian()
    tight = Hamiltonian(PhysicalParameters(vdw_settings=TIGHT))
    i, j = np.triu_indices(27, 1)
    rows = []
    for folder in sorted(run.iterdir()):
        snap = folder / "snapshots.npz"
        if not (folder / "complete.json").exists() or not snap.exists():
            continue
        with np.load(snap) as d:
            keep = np.nonzero(d["cycle"] > a.equil)[0]
            keep = keep[np.linspace(0, len(keep) - 1, min(a.max_snapshots, len(keep))).astype(int)]
            for k in keep:
                pos, Q = d["pos"][k], d["Q"][k]
                u_q = prod.vdw_pairs(pos[i], Q[i], pos[j], Q[j]).sum()
                u_t = tight.vdw_pairs(pos[i], Q[i], pos[j], Q[j]).sum()
                rows.append(dict(chain=folder.name, cycle=int(d["cycle"][k]),
                                 vdw_kBT=u_q / prod.kT, diff_kBT=(u_t - u_q) / prod.kT))
    diff = np.array([r["diff_kBT"] for r in rows])
    vdw = np.array([r["vdw_kBT"] for r in rows])
    w = np.exp(-(diff - diff.mean()))
    w /= w.mean()
    out = dict(snapshots=len(rows), production_settings=prod.phys.vdw_settings.__dict__,
               max_abs_diff_kBT=float(np.max(np.abs(diff))),
               mean_diff_kBT=float(diff.mean()), sd_diff_kBT=float(diff.std()),
               weight_min=float(w.min()), weight_max=float(w.max()),
               vdw_mean_kBT=float(vdw.mean()),
               vdw_mean_reweighted_kBT=float(np.sum(w * (vdw + diff)) / np.sum(w)),
               tight_settings=TIGHT.__dict__)
    (run / "quadrature_reweighting.json").write_text(json.dumps(dict(summary=out, rows=rows),
                                                               indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
