"""Holm adjustment for the full finite family of tests."""

import numpy as np


def holm(p):
    p = np.asarray(p, float)
    if p.ndim != 1 or np.any((p[np.isfinite(p)] < 0) | (p[np.isfinite(p)] > 1)):
        raise ValueError("Expected a vector of probabilities")
    result = np.full(p.shape, np.nan)
    finite = np.isfinite(p)
    values = p[finite]
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    adjusted[order] = np.minimum(
        1, np.maximum.accumulate(values[order] * (len(values) - np.arange(len(values))))
    )
    result[finite] = adjusted
    return result
