"""Continuation states for production v2 from the archived v1 run.

    python scripts/make_v2_sources.py outputs/production_v1_energy_hole outputs/v2_sources

v1 sampled a target without the D0 core minimum, so it is not a sample of
the corrected target, but its chains have relaxed for 3000 cycles from the
four lattice starts.  Each chain's final state is kept if the corrected
target allows it; otherwise the latest saved snapshot that it allows is used
(only compact_tilted40_s1, which fell into the energy hole, needs this).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from v904.model import Hamiltonian  # noqa: E402


def main():
    src, out = Path(sys.argv[1]), Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    h = Hamiltonian()
    i, j = np.triu_indices(27, 1)
    log = {}
    for folder in sorted(src.iterdir()):
        if not (folder / "complete.json").exists():
            continue
        with np.load(folder / "final_state.npz") as d:
            state = dict(pos=d["pos"], Q=d["Q"], mu=d["mu"])
        used = "final_state_cycle_3000"
        if not h.allowed(state["pos"], state["Q"]):
            with np.load(folder / "snapshots.npz") as s:
                for k in range(len(s["cycle"]) - 1, -1, -1):
                    if h.allowed(s["pos"][k], s["Q"][k]):
                        state = dict(pos=s["pos"][k], Q=s["Q"][k], mu=s["mu"][k])
                        used = f"snapshot_cycle_{int(s['cycle'][k])}"
                        break
                else:
                    raise RuntimeError(f"no allowed state in {folder.name}")
        vdw = h.vdw_pairs(state["pos"][i], state["Q"][i], state["pos"][j], state["Q"][j]).sum()
        (out / folder.name).mkdir(exist_ok=True)
        np.savez(out / folder.name / "final_state.npz", **state)
        log[folder.name] = dict(source=used, vdw_kBT=float(vdw / h.kT))
        print(folder.name, used, round(vdw / h.kT, 2))
    (out / "SOURCES.json").write_text(json.dumps(log, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
