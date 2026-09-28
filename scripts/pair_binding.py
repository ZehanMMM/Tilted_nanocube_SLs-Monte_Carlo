"""Is a single pair bound at 298 K?  Two-cube MC in the same cluster ensemble.

    python scripts/pair_binding.py OUT_DIR [--cycles 20000]

Full Hamiltonian, all moves, cube 0 pinned, cube 1 inside the bond sphere.
A pair is "in contact" when the steric support gap is below 2.5 nm (the
wall is at 1.8 nm).  With P = P(contact):

    K = V_free * P / (1 - P),   V_free = (4/3) pi R^3 - V_ex

is the association volume (bound configurations weighted by exp(-U/kT)),
independent of the shell size if contact really is a separate state.  K
compared with the volume per particle tells whether contacts survive: the
bound fraction at number density c is ~ K c / (1 + K c).
V_ex is the orientation-averaged excluded volume of two convex bodies
(Isihara-Kihara: V1 + V2 + (S1 M2 + S2 M1) / 4 pi; a cube of edge l has
S = 6 l^2 and integrated mean curvature M = 3 pi l, so V_ex = 11 l^3), taken
for the steric body l = L + gap.  It is an estimate; the check that matters
is whether K comes out the same for two bond radii, i.e. whether contact is
a genuine bound state rather than something set by the box.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from v904 import geometry  # noqa: E402
from v904.diagnostics import batch_means_mcse  # noqa: E402
from v904.model import ClusterEnsemble, E111, Hamiltonian, PhysicalParameters  # noqa: E402
from v904.sampler import Chain, MoveSettings  # noqa: E402
from scipy.spatial.transform import Rotation  # noqa: E402


def run(bond_nm, cycles, seed):
    ham = Hamiltonian(PhysicalParameters(), ClusterEnsemble(bond_nm=bond_nm))
    moves = json.loads((ROOT / "configs" / "production.json").read_text())["moves"]
    moves = replace(MoveSettings(), **moves)
    Q = Rotation.align_vectors([[1.0, 0, 0]], [E111])[0].as_matrix()
    # start bound: face-to-face along the body x axis at the steric wall
    state = dict(pos=np.array([[0, 0, 0], (Q @ np.array([1.0, 0, 0])) * (16e-9 + 1.85e-9)]),
                 Q=np.array([Q, Q]), mu=np.array([[1.0, 0, 0], [1.0, 0, 0]]))
    chain = Chain(ham, state, moves, np.random.default_rng(seed))
    gap, vdw, dd = np.empty(cycles), np.empty(cycles), np.empty(cycles)
    for c in range(cycles):
        chain.translate(1)
        chain.rotate(0)
        chain.rotate(1)
        chain.gamma()
        chain.scale()
        for _ in range(moves.n_mag_sweeps):
            chain.dipole(0)
            chain.dipole(1)
        r = chain.pos[1] - chain.pos[0]
        gap[c] = geometry.support_gap(r, chain.Q[0], chain.Q[1], ham.L, ham.r_round)[0] * 1e9
        vdw[c] = chain.W[0, 1] / ham.kT
        dd[c] = chain.E_dd / ham.kT
    keep = slice(cycles // 10, None)
    contact = (gap[keep] < 2.5).astype(float)
    p = float(contact.mean())
    R = ham.bond_m * 1e9
    l_nm = (ham.L + ham.gap0) * 1e9
    v_free = 4.0 / 3.0 * np.pi * R ** 3 - 11.0 * l_nm ** 3
    return dict(bond_nm=bond_nm, cycles=cycles, p_contact=p,
                p_contact_mcse=batch_means_mcse(contact), v_free_nm3=v_free,
                K_nm3=v_free * p / (1 - p) if p < 1 else float("inf"),
                vdw_mean_kBT=float(vdw[keep].mean()), dipole_mean_kBT=float(dd[keep].mean()),
                vdw_in_contact_kBT=float(vdw[keep][contact > 0].mean()) if p > 0 else None,
                acceptance=chain.acceptance(), drift_kBT=chain.drift_kBT())


def main():
    a = argparse.ArgumentParser()
    a.add_argument("out")
    a.add_argument("--cycles", type=int, default=20000)
    args = a.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = [run(b, args.cycles, 900 + k) for k, b in enumerate((30.0, 36.0))]
    (out / "pair_binding.json").write_text(json.dumps(rows, indent=2, default=float),
                                           encoding="utf-8")
    print(json.dumps(rows, indent=2, default=float))


if __name__ == "__main__":
    main()
