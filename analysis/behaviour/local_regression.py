"""Vectorised local family-CR1 regressions; nonfinite and constant predictors are unavailable."""

import numpy as np
from scipy import stats


def fit_maps(X, y, Z, groups):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    Z = np.asarray(Z, float)
    assert np.isfinite(y).all() and np.isfinite(Z).all()
    if np.linalg.matrix_rank(Z) != Z.shape[1]:
        raise ValueError("Rank-deficient covariates")
    if len(y) != X.shape[0] or len(groups) != len(y) or Z.shape[0] != len(y):
        raise ValueError("Misaligned local-regression inputs")
    Q = np.linalg.qr(Z, mode="reduced")[0]
    yr = y - Q @ (Q.T @ y)
    finite = np.isfinite(X).all(0)
    XR = np.full_like(X, np.nan)
    XR[:, finite] = X[:, finite] - Q @ (Q.T @ X[:, finite])
    ss = np.sum(XR * XR, axis=0)
    valid = finite & (ss > 1e-20)
    A = XR[:, valid]
    ss = ss[valid]
    beta = yr @ A / ss
    gs = np.unique(groups)
    n = len(y)
    k = Z.shape[1] + 1
    G = len(gs)
    assert G > k
    U = np.stack(
        [
            A[groups == g].T @ yr[groups == g]
            - beta * np.sum(A[groups == g] ** 2, axis=0)
            for g in gs
        ]
    )
    se = np.sqrt(G / (G - 1) * ((n - 1) / (n - k)) * np.sum(U * U, axis=0)) / ss
    t = beta / se
    p = 2 * stats.t.sf(abs(t), G - 1)
    c = stats.t.ppf(0.975, G - 1)
    result = {"valid": valid, "df": G - 1}
    for key, v in [
        ("beta", beta),
        ("se", se),
        ("t", t),
        ("p", p),
        ("ci_low", beta - c * se),
        ("ci_high", beta + c * se),
        ("r", yr @ A / np.sqrt(ss * (yr @ yr))),
    ]:
        arr = np.full(X.shape[1], np.nan)
        arr[valid] = v
        result[key] = arr
    return result
