import numpy as np
from scipy import sparse
from scipy.optimize import linear_sum_assignment


def operator(A, mask):
    lookup = np.full(mask.size, -1, int)
    lookup[mask.ravel(order="F")] = np.arange(mask.sum())
    lookup = lookup.reshape(mask.shape, order="F")
    rows = []
    cols = []
    weights = []
    for v in np.flatnonzero(A["support"]):
        x, y, _ = A["flat"][v]
        hi = int(v >= 32492)
        gx = (x - (-499 if hi == 0 else 1)) / 2
        gy = (y + 199) / 2
        j = int(np.floor(gx)) + hi * 250
        i = int(np.floor(gy))
        dx = gx - np.floor(gx)
        dy = gy - i
        pts = [(i, j), (i, j + 1), (i + 1, j), (i + 1, j + 1)]
        assert all((mask[z] for z in pts))
        for z, w in zip(
            pts, [(1 - dx) * (1 - dy), dx * (1 - dy), (1 - dx) * dy, dx * dy]
        ):
            rows.append(v)
            cols.append(lookup[z])
            weights.append(w)
    W = sparse.csr_matrix(
        (weights, (rows, cols)), shape=(len(A["support"]), int(mask.sum()))
    )
    rows = []
    cols = []
    weights = []
    for k, l in enumerate(A["parcel_ids"]):
        ix = np.flatnonzero((A["labels"] == l) & A["support"])
        w = A["area"][ix] / A["area"][ix].sum()
        rows.extend([k] * len(ix))
        cols.extend(ix)
        weights.extend(w)
    G = sparse.csr_matrix(
        (weights, (rows, cols)), shape=(len(A["parcel_ids"]), len(A["support"]))
    )
    V = G @ W
    assert np.allclose(np.asarray(V.sum(1)).ravel(), 1)
    return (W, V)


def spin_index(centers, hemi, R):
    out = np.empty(len(centers), int)
    reflect = np.diag([-1.0, 1.0, 1.0])
    cost = []
    for h, rot in enumerate([R, reflect @ R @ reflect]):
        ix = np.flatnonzero(hemi == h)
        cc = centers[ix]
        dist = np.linalg.norm((cc @ rot)[:, None] - cc[None], axis=2)
        i, j = linear_sum_assignment(dist)
        out[ix[i]] = ix[j]
        cost.extend(dist[i, j])
    return (out, float(np.mean(cost)))


def corrs(X, y):
    X = X - X.mean(0)
    y = y - y.mean()
    return X.T @ y / np.sqrt((X * X).sum(0) * (y @ y))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    v = np.minimum(1, np.maximum.accumulate(p[o] * (len(p) - np.arange(len(p)))))
    a = np.empty_like(p)
    a[o] = v
    return a
