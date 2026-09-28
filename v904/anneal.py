"""Simulated annealing: the MC kernels used as an optimiser of the total energy.

Sampling and optimisation ask different questions and obey different rules:

* sampling (T fixed, steps frozen) estimates averages over exp(-U/kT); it
  needs detailed balance, a fixed kernel and convergence diagnostics.
* annealing looks for argmin U.  The acceptance test is the same, but kT is
  lowered along a schedule and step sizes MAY adapt, because nothing is
  averaged: the output is the best state found, not a sample.

Two caveats a sampler does not have.  (1) The guarantee of reaching the
global minimum needs a logarithmic schedule, T_k ~ c / log k, far too slow
to use; a geometric schedule is a heuristic, so the result is checked by
independent restarts from different starts: agreement of their best
energies and structures is the evidence for a global minimum.  (2) The
minimum of U is the T = 0 answer; at 298 K entropy also matters, so the
annealed structure is compared with the finite-temperature distribution,
not substituted for it.
"""
from __future__ import annotations

import time
from dataclasses import asdict, replace

import numpy as np

from .sampler import Chain, MoveSettings, observe

# move -> step field adapted during annealing
ADAPT = dict(rot="rot_step_deg", latrot="latrot_step_deg", cotilt="cotilt_step_deg",
             gamma="gamma_step_deg", dip="dip_step_rad")
LIMITS = dict(rot_step_deg=(1e-3, 60.0), latrot_step_deg=(1e-3, 30.0),
              cotilt_step_deg=(1e-3, 90.0), gamma_step_deg=(1e-3, 30.0),
              dip_step_rad=(1e-4, 3.0))


def geometric_schedule(n, t_start, t_end):
    return t_start * (t_end / t_start) ** (np.arange(n) / max(n - 1, 1))


def anneal(ham, state, moves: MoveSettings, rng, factors, adapt_every=20,
           target_acceptance=0.3, progress=None, lattice=True):
    """Anneal along `factors` (kT multipliers per cycle); returns the best state.

    Every `adapt_every` cycles each step is multiplied by
    exp(acceptance - target), clipped to LIMITS, so steps shrink as the
    landscape is resolved more finely at low temperature.
    """
    chain = Chain(ham, state, moves, rng)
    step = chain.cycle_lattice if lattice else chain.cycle
    n = len(factors)
    keys = None
    traj = {}
    best_E = np.inf
    best = None
    t0 = time.perf_counter()
    for c in range(n):
        chain.temperature_factor = float(factors[c])
        step()
        obs = observe(chain, lattice)
        if keys is None:
            keys = list(obs)
            traj = {k: np.empty(n) for k in keys}
            traj["temperature_factor"] = np.asarray(factors, float)
        for k in keys:
            traj[k][c] = obs[k]
        E = chain.energy_terms()["Total"]
        if E < best_E:
            best_E = E
            best = dict(cycle=c + 1, pos=chain.pos.copy(), Q=chain.Q.copy(),
                        mu=chain.mu.copy(), G=chain.G.copy())
        if adapt_every and (c + 1) % adapt_every == 0:
            acc = chain.acceptance()
            new = {}
            for move, field in ADAPT.items():
                a = acc.get(move)
                if a is None or not np.isfinite(a) or chain.tries[move] == 0:
                    continue
                lo, hi = LIMITS[field]
                new[field] = float(np.clip(getattr(chain.moves, field)
                                           * np.exp(a - target_acceptance), lo, hi))
            chain.moves = replace(chain.moves, **new)
            chain.reset_counters()
        if progress is not None:
            progress(c + 1, n, time.perf_counter() - t0)
    return dict(traj=traj, chain=chain, best=best, best_energy_kBT=best_E / ham.kT,
                final_moves=asdict(chain.moves), elapsed_sec=time.perf_counter() - t0,
                drift_kBT=chain.drift_kBT())
