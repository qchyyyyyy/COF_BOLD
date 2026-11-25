"""Model-equal, method-adjusted quadratic association and model bootstrap."""

import numpy as np
from scipy.stats import t as student_t


def intervals(v):
    return np.quantile(v, [0.025, 0.975], axis=0)


def normalized_fraction(value):
    x = float(value)
    return "NA" if not np.isfinite(x) else format(x, ".12g")


def model_setting_mean(rows, field):
    trials = {}
    for r in rows:
        trials.setdefault(r["trial_index"], []).append(float(r[field]))
    return float(np.mean([np.mean(v) for v in trials.values()]))


def curve_model(x, y, method, model, seed, B):
    n = len(x)
    codes, gi = np.unique(model, return_inverse=True)
    G = len(codes)
    counts = np.bincount(gi)
    w = 1 / counts[gi]
    mu = np.average(x, weights=w)
    sd = np.sqrt(np.average((x - mu) ** 2, weights=w))
    assert sd > 0
    z = (x - mu) / sd
    labels = sorted(set(method))
    X = np.column_stack(
        [
            np.ones(n),
            z,
            z * z,
            *[(method == label).astype(float) for label in labels[1:]],
        ]
    )
    k = X.shape[1]
    assert G > k and np.linalg.matrix_rank(X) == k
    A = np.array([X[gi == g].T @ (w[gi == g, None] * X[gi == g]) for g in range(G)])
    b = np.array([X[gi == g].T @ (w[gi == g] * y[gi == g]) for g in range(G)])
    inv = np.linalg.inv(A.sum(axis=0))
    beta = inv @ b.sum(axis=0)
    res = y - X @ beta
    scores = np.array([X[gi == g].T @ (w[gi == g] * res[gi == g]) for g in range(G)])
    cov = G / (G - 1) * (n - 1) / (n - k) * inv @ scores.T @ scores @ inv
    se = np.sqrt(np.diag(cov))
    t = beta[2] / se[2]
    p = 2 * student_t.sf(abs(t), G - 1)
    rng = np.random.default_rng(seed)
    draws = []
    failed = 0
    for _ in range(B):
        m = rng.multinomial(G, np.ones(G) / G)
        aa = np.einsum("g,gij->ij", m, A)
        bb = m @ b
        if np.linalg.matrix_rank(aa) != k:
            failed += 1
            continue
        draws.append(np.linalg.solve(aa, bb)[2])
    assert len(draws) >= 0.99 * B
    return {
        "n_model_settings": n,
        "n_models": G,
        "quadratic_beta_standardized_x": float(beta[2]),
        "se_CR1": float(se[2]),
        "t": float(t),
        "df": G - 1,
        "p": float(p),
        "bootstrap_ci": intervals(np.array(draws)).tolist(),
        "bootstrap_rank_failures": failed,
        "condition_number": float(np.linalg.cond(X)),
        "method_levels": labels,
        "predictor_mean": float(mu),
        "predictor_sd": float(sd),
    }
