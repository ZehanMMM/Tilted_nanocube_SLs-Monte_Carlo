"""Convergence and precision diagnostics for (chain, draw) arrays.

Each function answers one question; none of them proves convergence.

  rank-normalised split R-hat  do the chains (and their halves) agree?
  bulk / tail ESS              how many independent draws is this worth?
  MCSE of the mean             how uncertain is the estimate from finite sampling?
  integrated autocorr. time    how many cycles per independent draw? (n / ESS)
  batch-means MCSE             the same MCSE by a different estimator
  Geweke z                     is the start of the kept window unlike the end?
  MSER                         where does initialisation bias stop?

R-hat, ESS and MCSE come from ArviZ 0.22 (Vehtari et al. 2021).  tau_int,
batch means, Geweke and MSER are implemented here so that two independent
routes to the MCSE can be compared: if they disagree, the error bar is not
trustworthy.
"""
from __future__ import annotations

import numpy as np


def require_arviz():
    try:
        import arviz as az
    except ImportError as exc:                   # pragma: no cover
        raise ImportError("pip install arviz==0.22.0") from exc
    return az


# ---------------------------------------------------------------------------
# single-chain tools
# ---------------------------------------------------------------------------
def autocorrelation(x):
    """Normalised autocorrelation rho(t) via FFT, t = 0..n-1."""
    x = np.asarray(x, float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    return acf / acf[0] if acf[0] > 0 else np.zeros(n)


def integrated_time(x, c=5.0):
    """Sokal's self-consistent window: tau = 1 + 2 sum_{t<=M} rho(t), M >= c tau.

    Returns (tau, window, reliable); reliable requires n > 50 tau, the usual
    rule for the estimate itself to be trustworthy.
    """
    rho = autocorrelation(x)
    taus = 2.0 * np.cumsum(rho) - 1.0
    m = np.arange(len(taus))
    ok = m >= c * taus
    window = int(np.argmax(ok)) if np.any(ok) else len(taus) - 1
    tau = float(max(taus[window], 1e-12))
    return tau, window, bool(len(x) > 50 * tau)


def batch_means_mcse(x, n_batches=None, batch_size=None):
    """MCSE of the mean from non-overlapping batch means.

    Default batch size sqrt(n) (Flegal & Jones 2010).  Batch means is only
    unbiased when a batch is much longer than the autocorrelation time, so
    `diagnose` passes batch_size = max(sqrt(n), 5 tau) instead (bias of the
    batch variance ~ tau / batch length, so <~ 10 % in the MCSE).
    """
    x = np.asarray(x, float)
    n = len(x)
    if batch_size is None:
        n_batches = n_batches or max(int(np.sqrt(n)), 2)
        batch_size = n // n_batches
    else:
        n_batches = n // batch_size
    if batch_size < 1 or n_batches < 2:
        return np.nan
    means = x[: batch_size * n_batches].reshape(n_batches, batch_size).mean(axis=1)
    return float(np.std(means, ddof=1) / np.sqrt(n_batches))


def chain_tau(x):
    """tau = n / ESS_mean for one chain, ESS from ArviZ (Geyer's initial
    monotone sequence).  Less biased for short series than a Sokal window."""
    az = require_arviz()
    x = np.asarray(x, float)
    return float(len(x) / az.ess(x[None, :], method="mean"))


def geweke_z(x, first=0.1, last=0.5):
    """(mean of first 10 %) - (mean of last 50 %) in units of its standard error.

    Each segment's variance of the mean is var / ESS (Geweke 1992 uses the
    spectral density at zero; ESS is the same quantity).  Naive batch means
    or a Sokal window underestimate it for segments only a few tau long and
    inflate |z|: on stationary AR(1) data they gave sd(z) = 1.44, this gives
    1.06-1.10 (see docs).
    """
    x = np.asarray(x, float)
    n = len(x)
    a, b = x[: int(first * n)], x[int((1 - last) * n):]

    def var_of_mean(seg):
        return np.var(seg) * chain_tau(seg) / len(seg)

    se = np.sqrt(var_of_mean(a) + var_of_mean(b))
    return float((a.mean() - b.mean()) / se) if se > 0 else np.nan


def mser_burn_in(x, batch=5, max_fraction=0.5):
    """MSER-5 truncation point (White 1997): argmin_d var(x[d:]) / (n - d)^2.

    Computed on batch means of size 5 and searched over the first half.
    """
    x = np.asarray(x, float)
    nb = len(x) // batch
    y = x[: nb * batch].reshape(nb, batch).mean(axis=1)
    best, arg = np.inf, 0
    for d in range(int(max_fraction * nb)):
        tail = y[d:]
        score = np.sum((tail - tail.mean()) ** 2) / len(tail) ** 2
        if score < best:
            best, arg = score, d
    return arg * batch


# ---------------------------------------------------------------------------
# multi-chain table
# ---------------------------------------------------------------------------
def diagnose(draws: dict, rhat_max=1.01, ess_min=400, ess_per_chain=100,
             mcse_targets=None, min_chains=4, min_draws=100):
    """One row per variable; `status` lists every failed check.

    draws[name] has shape (chain, draw) after warm-up, never concatenated.
    """
    az = require_arviz()
    mcse_targets = mcse_targets or {}
    rows = []
    for name, values in draws.items():
        x = np.asarray(values, float)
        chains, n = x.shape
        row = dict(variable=name, chains=chains, draws=n, mean=float(np.mean(x)),
                   sd=float(np.std(x)), rhat=np.nan, ess_bulk=np.nan, ess_tail=np.nan,
                   ess_mean=np.nan, mcse_mean=np.nan, mcse_batch=np.nan, batch_size=np.nan,
                   tau_int_max=np.nan, tau_reliable=False, geweke_max_abs=np.nan,
                   ess_target=max(ess_min, ess_per_chain * chains),
                   mcse_target=mcse_targets.get(name, np.nan))
        reasons = []
        if not np.all(np.isfinite(x)):
            reasons.append("nonfinite")
        elif chains < 2 or n < min_draws:
            reasons.append("too_few_chains_or_draws")
        elif np.any(np.ptp(x, axis=1) == 0):
            reasons.append("constant_chain")
        else:
            row.update(rhat=float(az.rhat(x, method="rank")),
                       ess_bulk=float(az.ess(x, method="bulk")),
                       ess_tail=float(az.ess(x, method="tail")),
                       ess_mean=float(az.ess(x, method="mean")),
                       mcse_mean=float(az.mcse(x, method="mean")))
            # independent route: per-chain batch means (batches >= 5 tau long,
            # at least 10 of them), pooled as independent chains
            taus = [chain_tau(c) for c in x]
            size = max(int(np.sqrt(n)), int(np.ceil(5 * max(taus))))
            if n // size >= 10:
                per_chain = np.array([batch_means_mcse(c, batch_size=size) for c in x])
                row["mcse_batch"] = float(np.sqrt(np.sum(per_chain ** 2)) / chains)
            row["batch_size"] = size
            row["tau_int_max"] = float(max(taus))
            row["tau_reliable"] = bool(n > 50 * max(taus))
            row["geweke_max_abs"] = float(np.nanmax(np.abs([geweke_z(c) for c in x])))
            if chains < min_chains:
                reasons.append("fewer_than_4_chains")
            if row["rhat"] >= rhat_max:
                reasons.append("high_rhat")
            if min(row["ess_bulk"], row["ess_tail"]) < row["ess_target"]:
                reasons.append("low_ess")
            if np.isfinite(row["mcse_target"]) and row["mcse_mean"] > row["mcse_target"]:
                reasons.append("high_mcse")
        row["status"] = ",".join(reasons) if reasons else "passed"
        rows.append(row)
    return rows
