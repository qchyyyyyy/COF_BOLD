"""OLS coefficients with the manuscript's family CR1 covariance."""

import numpy as np
from scipy import stats


def ols(X, y, clusters):
    X, y, clusters = (np.asarray(X, float), np.asarray(y, float), np.asarray(clusters))
    n, k = X.shape
    if (
        y.shape != (n,)
        or clusters.shape != (n,)
        or (not np.isfinite(X).all())
        or (not np.isfinite(y).all())
    ):
        raise ValueError("Invalid or misaligned regression inputs")
    if n <= k or np.linalg.matrix_rank(X) < k:
        raise ValueError("Rank-deficient model or insufficient observations")
    groups = np.unique(clusters)
    G = len(groups)
    if G <= k:
        raise ValueError("Insufficient independent clusters")
    inverse = np.linalg.inv(X.T @ X)
    beta = inverse @ X.T @ y
    residual = y - X @ beta
    scores = np.array([X[clusters == g].T @ residual[clusters == g] for g in groups])
    covariance = (
        inverse @ scores.T @ scores @ inverse * (G / (G - 1)) * ((n - 1) / (n - k))
    )
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    p = 2 * stats.t.sf(abs(beta / np.maximum(se, 1e-30)), G - 1)
    return dict(
        beta=beta,
        cov=covariance,
        se=se,
        p=p,
        df=G - 1,
        n=n,
        method="family_CR1",
        residual=residual,
    )


def design(x, age, sex, mu, sd, quadratic=True):
    z = (np.asarray(x) - mu) / sd
    return np.column_stack(
        [np.ones(len(z)), z] + ([z * z] if quadratic else []) + [age, sex]
    )
