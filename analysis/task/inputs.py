"""Task labels, inclusive MATLAB block boundaries and condition means."""

from pathlib import Path
import hashlib
import numpy as np
from analysis.common.matlab import read_mat, scalar_string

METRIC_FIELDS = {
    "PQ_diverging": ("numPattern", "NumSource"),
    "PQ_converging": ("numPattern", "NumSink"),
    "PQ_cc_spiral": ("numPattern", "NumCCVortex"),
    "PQ_c_spiral": ("numPattern", "NumCVortex"),
    "PS_diverging": ("sizePattern", "SizeSource"),
    "PS_converging": ("sizePattern", "SizeSink"),
    "PS_cc_spiral": ("sizePattern", "SizeCCVortex"),
    "PS_c_spiral": ("sizePattern", "SizeCVortex"),
    "OP": ("orderParam", "order_param_global"),
}


def task_paths(cfg, task):
    run = "tfMRI_" + task["key"] + "_LR"
    group = Path(cfg["input_root"]) / run / "group"
    paths = {
        "subjects": group / (run + "_SubList.mat"),
        "events": group / (run + "_EVs.mat"),
    }
    for kind in ("numPattern", "sizePattern", "orderParam"):
        paths[kind] = group / (run + "_" + kind + "_global_LR_COF.mat")
    return paths


def load_task(cfg, task):
    paths = task_paths(cfg, task)
    sub = read_mat(paths["subjects"])
    ids = [
        scalar_string(x)
        for x in next((v for k, v in sub.items() if k.lower() == "subjectlist")).ravel()
    ]
    assert len(ids) == task["n"] and len(set(ids)) == len(
        ids
    ), "Subject count/uniqueness failed"
    ev = read_mat(paths["events"])
    events = np.asarray(ev["Events"], dtype=np.int16)
    stim = np.asarray(ev["Stim_Type_subs"], dtype=np.int64)
    names = [scalar_string(x) for x in ev["condition_names"].ravel()]
    assert names == task["ev_names"], ("Condition names changed", task["key"], names)
    assert events.shape == (task["n"], task["t"])
    assert set(np.unique(events)) == {0, 1, 2}
    assert stim.shape[1:] == (3, task["n"])
    reconstructed = np.zeros_like(events)
    occupancy = np.zeros_like(events)
    for i in range(task["n"]):
        for start, end, label in stim[:, :, i]:
            if label == 0:
                assert start == 0 and end == 0, "Nonzero unlabeled Stim row"
                continue
            assert label in (1, 2)
            assert end >= start or start > task["t"] or end < 1
            lo, hi = (max(1, int(start)), min(task["t"], int(end)))
            if lo <= hi:
                reconstructed[i, lo - 1 : hi] = label
                occupancy[i, lo - 1 : hi] += 1
    assert np.array_equal(events, reconstructed), (
        "Stored Events/Stim mismatch",
        task["key"],
    )
    metrics = {}
    for kind in ("numPattern", "sizePattern", "orderParam"):
        d = read_mat(paths[kind])
        for name in task["metrics"]:
            file_kind, field = METRIC_FIELDS[name]
            if kind == file_kind:
                x = np.asarray(d[field], dtype=float)
                assert x.shape == events.shape
                metrics[name] = x
    valid = np.ones_like(events, dtype=bool)
    valid[:, : cfg["window"]["trim_first"]] = False
    if cfg["window"]["trim_last"]:
        valid[:, -cfg["window"]["trim_last"] :] = False
    overlaps = occupancy > 1
    summary = {
        "task": task["key"],
        "n": len(ids),
        "frames": task["t"],
        "condition_names": names,
        "labels": {"0": "fixation", "1": task["control"], "2": task["task"]},
        "subject_order_sha256": hashlib.sha256(
            np.array(ids, dtype="<f8").tobytes()
        ).hexdigest(),
        "event_reconstruction_mismatch_frames": int(
            np.count_nonzero(events != reconstructed)
        ),
        "overlap_runs": int(np.any(overlaps, axis=1).sum()),
        "overlap_frames": int(overlaps.sum()),
        "condition_frame_minmax": {
            str(c): [
                int(np.sum((events == c) & valid, axis=1).min()),
                int(np.sum((events == c) & valid, axis=1).max()),
            ]
            for c in (0, 1, 2)
        },
        "participants_missing_condition_after_trim": {
            str(c): int(np.count_nonzero(np.sum((events == c) & valid, axis=1) == 0))
            for c in (0, 1, 2)
        },
        "metrics": {
            name: {
                "nonfinite": int((~np.isfinite(x)).sum()),
                "all_zero_rows": int(np.all(x == 0, axis=1).sum()),
            }
            for name, x in metrics.items()
        },
    }
    if task["key"] == "WM" and "wm_expected_overlap_runs" in cfg["window"]:
        assert summary["overlap_runs"] == cfg["window"]["wm_expected_overlap_runs"]
        assert summary["overlap_frames"] == cfg["window"]["wm_expected_overlap_frames"]
    return (metrics, events, valid, overlaps, summary)


def condition_means(x, events, valid):
    means = np.full((x.shape[0], 3), np.nan)
    counts = np.zeros_like(means, dtype=int)
    for c in (0, 1, 2):
        mask = (events == c) & valid
        counts[:, c] = mask.sum(axis=1)
        sums = np.where(mask, x, 0.0).sum(axis=1)
        np.divide(sums, counts[:, c], out=means[:, c], where=counts[:, c] > 0)
    return (means, counts)
