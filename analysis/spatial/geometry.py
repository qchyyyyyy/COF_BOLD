"""Prepare atlas geometry, vertex area and the fixed working-memory annotation."""

import numpy as np
import h5py
from scipy.io import loadmat
from scipy import sparse
from analysis.common.io import path, output_directory, Inputs, dump


def run(C):
    import nibabel as nib

    T = path(C, "templates_root")
    inputs = Inputs(C)
    for filename in [
        T / "BAT/cifti_database_data.mat",
        T / "coordinate_cii_v2.mat",
        T
        / "HCP_S1200_Atlas_Z4_pkXDZ/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii",
    ]:
        inputs.add(filename)
    for hemisphere in ["L", "R"]:
        for stem in ["flat", "sphere", "midthickness_MSMAll"]:
            inputs.add(
                T
                / "HCP_S1200_Atlas_Z4_pkXDZ"
                / f"S1200.{hemisphere}.{stem}.32k_fs_LR.surf.gii"
            )
    L = T / "HCP_S1200_Atlas_Z4_pkXDZ"
    parts = []
    with h5py.File(T / "coordinate_cii_v2.mat") as h:
        for side in ["left", "right"]:
            a = np.asarray(h["cii_mask_" + side]).T
            parts.append(
                (a[:-1, :-1] == 1)
                & (a[1:, :-1] == 1)
                & (a[:-1, 1:] == 1)
                & (a[1:, 1:] == 1)
            )
    mask = np.hstack(parts)
    assert mask.shape == (200, 500) and mask.sum() == 45990
    N = 32492
    D = loadmat(T / "BAT/cifti_database_data.mat", simplify_cells=True)
    terms = [d["name"].split("_pAgF")[0] for d in D["dir_path"]]
    B = np.zeros((2 * N, len(terms)))
    agree = []
    for j in range(len(terms)):
        v = np.atleast_2d(D["term_surface_map_res_all"][j])
        ix = v[:, 0].astype(int) - 1
        assert (
            np.all(v[:, 0] == ix + 1)
            and len(np.unique(ix)) == len(ix)
            and (ix.min() >= 0)
            and (ix.max() < 2 * N)
        )
        B[ix, j] = v[:, 1]
        for hi, side in enumerate(["left", "right"]):
            w = np.atleast_2d(D["term_surface_map_res_" + side][j])
            ii = w[:, 0].astype(int) - 1
            if hi == 1 and ii.min() >= N:
                ii = ii - N
            assert ii.min() >= 0 and ii.max() < N
            tmp = np.zeros(N)
            tmp[ii] = w[:, 1]
            agree.append(bool(np.array_equal(tmp, B[hi * N : (hi + 1) * N, j])))
    assert all(agree), "BAT all/left/right index disagreement"
    flat = []
    sphere = []
    faces = []
    areas = []
    geometry = []
    with h5py.File(T / "coordinate_cii_v2.mat") as f:
        for hi, (hemi, side) in enumerate([("L", "left"), ("R", "right")]):
            base = T / "HCP_S1200_Atlas_Z4_pkXDZ"
            x, ff = nib.load(base / f"S1200.{hemi}.flat.32k_fs_LR.surf.gii").agg_data()
            s, sf = nib.load(
                base / f"S1200.{hemi}.sphere.32k_fs_LR.surf.gii"
            ).agg_data()
            mid, mf = nib.load(
                base / f"S1200.{hemi}.midthickness_MSMAll.32k_fs_LR.surf.gii"
            ).agg_data()
            assert np.array_equal(sf, mf)
            st = {tuple(sorted(t)) for t in sf}
            assert all((tuple(sorted(t)) in st for t in ff))
            bad = np.asarray(f["i_invalid_" + side]).ravel().astype(int) - 1
            valid = np.setdiff1d(np.arange(N), bad)
            xx = np.asarray(f["x_" + side]).ravel()
            yy = np.asarray(f["y_" + side]).ravel()
            err = min(
                np.max(np.abs(xx - x[valid, 0])),
                np.max(np.abs(xx - (x[valid, 0] + (-250 if hi == 0 else 250)))),
            )
            assert err < 0.001 and np.max(np.abs(yy - x[valid, 1])) < 0.001, (
                "coordinate order mismatch",
                hi,
                err,
            )
            area = np.zeros(N)
            ta = (
                np.linalg.norm(
                    np.cross(
                        mid[mf[:, 1]] - mid[mf[:, 0]], mid[mf[:, 2]] - mid[mf[:, 0]]
                    ),
                    axis=1,
                )
                / 6
            )
            for k in range(3):
                np.add.at(area, mf[:, k], ta)
            x = x.astype(float)
            x[:, 0] += -250 if hi == 0 else 250
            flat.append(x)
            sphere.append(s / np.linalg.norm(s, axis=1, keepdims=True))
            faces.append(sf)
            areas.append(area)
            geometry.append(
                {
                    "hemisphere": hemi,
                    "coordinate_order_max_error": float(err),
                    "flat_faces_subset_of_sphere": True,
                    "invalid_vertices": len(bad),
                }
            )
    flat = np.vstack(flat)
    sphere = np.vstack(sphere)
    area = np.concatenate(areas)
    labfile = (
        L
        / "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii"
    )
    ni = nib.load(labfile)
    lab = np.zeros(2 * N, int)
    val = ni.get_fdata()[0]
    for nm, sl, bm in ni.header.get_axis(1).iter_structures():
        lab[bm.vertex + (0 if "LEFT" in nm else N)] = val[sl].astype(int)
    labels = sorted(np.unique(lab[lab > 0]))
    assert len(labels) == 360
    rows = []
    cols = []
    weights = []
    support = np.zeros(2 * N, bool)
    lookup = np.full(mask.shape, -1, int)
    lookup[mask] = np.arange(mask.sum())
    for v, (x, y, _) in enumerate(flat):
        hi = int(v >= N)
        gx = (x - (-499 if hi == 0 else 1)) / 2
        gy = (y + 199) / 2
        j = int(np.floor(gx))
        i = int(np.floor(gy))
        j += hi * 250
        if i < 0 or i >= 199 or j < hi * 250 or (j > hi * 250 + 248) or (lab[v] == 0):
            continue
        pp = [(i, j), (i, j + 1), (i + 1, j), (i + 1, j + 1)]
        if not all((mask[z] for z in pp)):
            continue
        dx = gx - np.floor(gx)
        dy = gy - i
        ww = [(1 - dx) * (1 - dy), dx * (1 - dy), (1 - dx) * dy, dx * dy]
        support[v] = True
        for z, w in zip(pp, ww):
            rows.append(v)
            cols.append(lookup[z])
            weights.append(w)
    W = sparse.csr_matrix((weights, (rows, cols)), shape=(2 * N, mask.sum()))
    assert np.allclose(np.asarray(W.sum(1)).ravel()[support], 1)
    centers = []
    coverage = []
    rr = []
    cc = []
    ww = []
    for k, l in enumerate(labels):
        vv = np.flatnonzero((lab == l) & support)
        full = np.flatnonzero(lab == l)
        cov = area[vv].sum() / area[full].sum()
        assert cov >= 0.5, (l, cov)
        coverage.append(cov)
        w = area[vv] / area[vv].sum()
        rr.extend([k] * len(vv))
        cc.extend(vv)
        ww.extend(w)
        cen = np.average(sphere[vv], axis=0, weights=w)
        centers.append(cen / np.linalg.norm(cen))
    G = sparse.csr_matrix((ww, (rr, cc)), shape=(360, 2 * N))
    BW = G @ B
    V = G @ W
    out = output_directory(C)
    np.savez_compressed(
        out / "geometry.npz",
        B=B,
        sphere=sphere,
        area=area,
        support=support,
        labels=lab,
        centers=np.array(centers),
        parcel_ids=labels,
        parcel_hemi=np.array(
            [0 if np.all(np.flatnonzero(lab == label) < N) else 1 for label in labels]
        ),
        terms=np.array(terms),
        flat=flat,
    )
    dump(
        out / "summary.json",
        {
            "hemispheres": geometry,
            "support": int(support.sum()),
            "parcels": len(labels),
        },
    )
    inputs.finish(out)
