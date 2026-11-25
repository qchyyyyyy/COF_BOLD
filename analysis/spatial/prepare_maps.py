"""Prepare manuscript regional maps and the twelve fixed-ROI family-CR1 contrasts."""

import numpy as np
import pandas as pd
import h5py
from scipy.io import loadmat
from scipy import stats
from analysis.common.io import path, output_directory, Inputs, dump
from analysis.common.regression import ols
from .operators import operator, holm, corrs


def run(C):
    inputs = Inputs(C)
    out = output_directory(C)
    A = np.load(inputs.add(path(C, "geometry_npz")))
    F = loadmat(inputs.add(path(C, "task_effects_mat")))
    mask = F["mask"].astype(bool)
    assert mask.shape == (200, 500) and mask.sum() == 45990
    W, V = operator(A, mask)
    term = np.flatnonzero(A["terms"] == C["roi"])
    assert len(term) == 1
    support = A["support"]
    roi = (A["B"][:, term[0]] > 0) & support
    other = ~roi & support
    assert roi.sum() > 100 and other.sum() > 100
    w = A["area"] * roi
    w = w / w.sum()
    v = A["area"] * other
    v = v / v.sum()
    roiweight = np.asarray((w - v) @ W).ravel()
    T = path(C, "wm_task_root")
    sp = T / "group/tfMRI_WM_LR_SubList.mat"
    S = loadmat(inputs.add(sp), simplify_cells=True)
    ids = S.get("subjectList", S.get("SubjectList"))
    assert len(ids) == C["expected_n"] and len(set(ids)) == len(ids)
    FT = (
        pd.read_csv(inputs.add(path(C, "family_csv")))
        .set_index("Subject")
        .reindex([int(i) for i in ids])
    )
    assert FT.Family_ID.notna().all()
    families = pd.factorize(FT.Family_ID)[0]
    cache = np.stack(
        [
            F["maps"][:, :, k, c].ravel(order="F")[mask.ravel(order="F")]
            for c in range(3)
            for k in range(4)
        ],
        axis=1,
    ).reshape(45990, 3, 4)
    mean = np.zeros((45990, 3, 4))
    m2 = np.zeros_like(mean)
    values = []
    for i, sid in enumerate(ids):
        sid = str(sid)
        density_path = (
            T
            / "participant"
            / sid
            / (sid + "_tfMRI_WM_LR_numPattern_level1_LR_COF.mat")
        )
        op = (
            T
            / "participant"
            / sid
            / (sid + "_tfMRI_WM_LR_orderParam_level1_LR_COF.mat")
        )
        x = loadmat(inputs.add(density_path), variable_names=["cope"])["cope"].reshape(
            45990, 3, 4, order="F"
        )[:, :, :2]
        o = loadmat(inputs.add(op), variable_names=["cope"])["cope"].reshape(
            45990, 3, order="F"
        )
        xx = np.stack([x[:, :, 0], x[:, :, 1], x.sum(2), o], axis=2)
        assert np.isfinite(xx).all()
        delta = xx - mean
        mean += delta / (i + 1)
        m2 += delta * (xx - mean)
        values.append(np.einsum("v,vcm->cm", roiweight, xx))
        if (i + 1) % 100 == 0:
            print("COPE", i + 1, "/" + str(len(ids)), flush=True)
    dz = mean / (np.sqrt(m2 / (len(ids) - 1)) + np.finfo(float).eps)
    assert np.allclose(dz[:, :, [0, 1, 3]], cache[:, :, :3], atol=1e-10, rtol=1e-10)
    reg = V @ dz.reshape(45990, 12)
    bold = V @ cache[:, :, 3]
    assert np.isfinite(reg).all() and np.isfinite(bold).all()
    with h5py.File(inputs.add(path(C, "mode_render_mat"))) as h:
        refs = h["args"][:, 0]
        M = h[refs[0]]
        geo = h[refs[1]]
        assert np.array_equal(np.asarray(geo["V_mask"]).T.astype(bool), mask)
        names = [
            "".join((chr(int(t)) for t in np.asarray(h[r]).ravel()))
            for r in M["method"][:, 0]
        ]
        ix = names.index("COF")
        div = np.asarray(h[M["Div"][ix, 0]]).T
        grad = np.asarray(h[refs[3]]).T
        assert np.isfinite(div[mask]).all() and np.isfinite(grad[mask]).all()
        x = div.ravel(order="F")[mask.ravel(order="F")]
        y = grad.ravel(order="F")[mask.ravel(order="F")]
        tail = np.where(abs(x) >= 2, np.sign(x), 0)
        fig1 = V @ np.column_stack([x, tail])
        gradient = V @ y
        saved_div_r = float(np.asarray(h[M["div_r"][ix, 0]]).squeeze())
        assert abs(corrs(tail[:, None], y)[0] - saved_div_r) < 1e-10
    vv = np.array(values)
    roirows = []
    for c in range(3):
        for m, name in enumerate(["D", "C", "D+C", "OP"]):
            f = ols(np.ones((len(ids), 1)), vv[:, c, m], families)
            t = stats.t.ppf(0.975, f["df"])
            roirows.append(
                {
                    "contrast": c,
                    "metric": name,
                    "mean_ROI_minus_complement": float(f["beta"][0]),
                    "ci": [
                        float(f["beta"][0] - t * f["se"][0]),
                        float(f["beta"][0] + t * f["se"][0]),
                    ],
                    "p": float(f["p"][0]),
                    "families": len(np.unique(families)),
                    "n": len(ids),
                    "df": f["df"],
                }
            )
    for r, q in zip(roirows, holm([r["p"] for r in roirows])):
        r["p_holm12"] = float(q)
    np.savez_compressed(
        out / "prepared_maps.npz",
        regional=reg.reshape(360, 3, 4),
        bold=bold,
        fig1=fig1,
        gradient=gradient,
        centers=A["centers"],
        hemi=A["parcel_hemi"],
        parcel_ids=A["parcel_ids"],
    )
    dump(out / "roi_tests.json", roirows)
    inputs.finish(out)
