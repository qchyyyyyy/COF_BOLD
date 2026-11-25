"""Participant-equal paired effects with family-cluster uncertainty.

The primary P value uses an intercept-only CR1 covariance and a t(G-1)
reference. It is not a permutation P value. Whole-family bootstrap intervals
are pointwise percentile intervals, not simultaneous inference.
"""

from __future__ import annotations
import numpy as np
from scipy.stats import t as student_t


def bh(p):
    """Benjamini-Hochberg adjusted P values in input order; finite inputs only."""
    p = np.asarray(p, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or np.any((p < 0) | (p > 1)):
        raise ValueError("p must be a finite one-dimensional probability array")
    if not p.size:
        return p.copy()
    order = np.argsort(p, kind="stable")
    adjusted = np.minimum.accumulate(
        (p[order] * p.size / np.arange(1, p.size + 1))[::-1]
    )[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def _validate(d, family):
    d = np.asarray(d, dtype=float)
    family = np.asarray(family)
    if d.ndim != 1 or family.ndim != 1 or d.size != family.size:
        raise ValueError("d and family must be aligned one-dimensional arrays")
    if d.size < 3 or not np.isfinite(d).all():
        raise ValueError("at least three finite paired differences are required")
    for f in family:
        if f is None or str(f).strip().lower() in ("", "nan", "none", "nat"):
            raise ValueError(
                "missing family membership is not an independent singleton"
            )
        if isinstance(f, (float, np.floating)) and (not np.isfinite(f)):
            raise ValueError("nonfinite family identifier")
    try:
        _, inverse = np.unique(family, return_inverse=True)
    except TypeError as exc:
        raise ValueError(
            "family identifiers must have a consistent sortable type"
        ) from exc
    g = int(inverse.max()) + 1
    if g < 2:
        raise ValueError("at least two families are required")
    return (d, inverse, g)


def _percentile(values):
    finite = np.asarray(values)[np.isfinite(values)]
    if not finite.size:
        return (float("nan"), float("nan"))
    return tuple((float(x) for x in np.percentile(finite, [2.5, 97.5])))


def analyze(d, family, seed, bootstrap=10000):
    """Estimate the participant-equal mean of within-participant differences.

    Bootstrap draws sample G families uniformly with replacement, preserving
    all their members. Each draw's mean is its sum divided by its participant
    count, rather than the unweighted mean of family means. Its sample SD uses
    the resampled participant count minus one. Undefined nonzero/zero-SD dz
    draws are excluded from dz intervals and explicitly counted; the interval
    is conditional on finite draws. An all-zero draw has dz=0 by convention.

    Zero cluster variance with a nonzero effect has no calibrated finite t
    inference here: t is signed infinity and p is NaN, with explicit status.
    """
    d, inverse, g = _validate(d, family)
    if not isinstance(bootstrap, (int, np.integer)) or bootstrap < 1:
        raise ValueError("bootstrap must be a positive integer")
    n = d.size
    mean = float(d[0]) if np.all(d == d[0]) else float(d.mean())
    centered = d - mean
    counts = np.bincount(inverse).astype(float)
    sums = np.bincount(inverse, weights=centered)
    squares = np.bincount(inverse, weights=centered * centered)
    variance = float(np.dot(centered, centered) / (n - 1))
    sd = float(np.sqrt(variance))
    dz = mean / sd if sd > 0 else 0.0 if mean == 0 else float("nan")
    se = float(np.sqrt(g / (g - 1.0) * np.dot(sums, sums) / (n * n)))
    if se > 0:
        statistic = mean / se
        p = float(2 * student_t.sf(abs(statistic), g - 1))
        status = "ok"
    elif mean == 0:
        statistic, p = (0.0, 1.0)
        status = "all_zero" if np.all(d == 0) else "zero_cluster_variance_zero_mean"
    else:
        statistic, p = (float(np.copysign(np.inf, mean)), float("nan"))
        status = "zero_cluster_variance_nonzero_mean"
    rng = np.random.default_rng(seed)
    boot_mean = np.empty(bootstrap)
    boot_dz = np.empty(bootstrap)
    batch = max(1, min(256, 1000000 // g))
    for start in range(0, bootstrap, batch):
        end = min(bootstrap, start + batch)
        indices = rng.integers(0, g, size=(end - start, g))
        bn = counts[indices].sum(axis=1)
        bs = sums[indices].sum(axis=1)
        bq = squares[indices].sum(axis=1)
        bm = mean + bs / bn
        bv = np.maximum((bq - bs * bs / bn) / (bn - 1), 0.0)
        bsd = np.sqrt(bv)
        bd = np.full_like(bm, np.nan)
        np.divide(bm, bsd, out=bd, where=bsd > 0)
        bd[(bsd == 0) & (bm == 0)] = 0.0
        boot_mean[start:end], boot_dz[start:end] = (bm, bd)
    ci_low, ci_high = _percentile(boot_mean)
    dz_low, dz_high = _percentile(boot_dz)
    return {
        "n": int(n),
        "G": g,
        "mean_diff": mean,
        "dz": float(dz),
        "se": se,
        "t": float(statistic),
        "df": g - 1,
        "p": p,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "dz_ci_low": dz_low,
        "dz_ci_high": dz_high,
        "status": status,
        "bootstrap_draws": int(bootstrap),
        "bootstrap_dz_undefined": int((~np.isfinite(boot_dz)).sum()),
        "bootstrap_mean": float(boot_mean.mean()),
    }
