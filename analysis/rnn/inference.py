"""Fourteen paired scaling contrasts and two model-clustered curvature tests."""

from pathlib import Path
import csv
import numpy as np
from scipy.io import loadmat
from analysis.common.io import resolve_path, output_directory, Inputs, dump
from analysis.common.family import analyze
from analysis.common.multiple_testing import holm
from .models import intervals, normalized_fraction, model_setting_mean, curve_model


def run(cfg):
    inputs = Inputs(cfg)
    out = output_directory(cfg)
    B = cfg["bootstrap"]
    scaling = []
    for ii, item in enumerate(cfg["scaling"]["inputs"]):
        path = resolve_path(cfg, item["path"])
        inputs.add(path, item.get("sha256"))
        keys = (
            ["NumSource", "NumSink"]
            if item["metric"] == "PQ"
            else ["SizeSource", "SizeSink"]
        )
        mat = loadmat(path, variable_names=keys)
        a, b = (mat[keys[0]], mat[keys[1]])
        assert list(a.shape) == item["shape"] and a.shape == b.shape
        assert np.isfinite(a).all() and np.isfinite(b).all()
        values = (a.mean(axis=(0, 1)) + b.mean(axis=(0, 1))) / 2
        d = values[2] - values[1]
        r = analyze(d, np.arange(len(d)), cfg["seed"] + ii, B)
        pooled = values[1:].reshape(-1)
        s = pooled.std(ddof=1)
        r.update(
            size=item["size"],
            metric=item["metric"],
            mean_0back=float(values[1].mean()),
            mean_2back=float(values[2].mean()),
            standardized_difference=float(d.mean() / s),
            mean_diff_raw=float(d.mean()),
        )
        scaling.append(r)
    for r, q in zip(
        scaling, holm([r["p"] if np.isfinite(r["p"]) else 1 for r in scaling])
    ):
        r["p_holm"] = float(q)
    assert len(scaling) == cfg["expected_scaling_tests"]
    print("Scaling14 completed", flush=True)
    entry = cfg["input"]
    path = resolve_path(cfg, entry["path"])
    inputs.add(path, entry.get("sha256"))
    with path.open(encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == entry["rows"]
    models = sorted({r["subject"] for r in rows})
    assert len(models) == entry["unique_subjects"]
    grouped = {}
    seen = set()
    index_map = {}
    clean = {}
    fields = ["accuracy", "PQ_DC_mean", "PS_DC_mean"]
    for r in rows:
        r["trial_index"] = int(float(r["trial_index"]))
        setting = (
            r["method"],
            float(r["severity"]),
            normalized_fraction(r["fraction"]),
        )
        key = (r["subject"], r["trial_index"], setting, int(float(r["repeat"])))
        assert key not in seen
        seen.add(key)
        index_map.setdefault(r["subject"], set()).add(r["subject_num"])
        for f in fields:
            v, c, d = map(float, [r[f], r["clean_" + f], r["delta_" + f]])
            assert np.isfinite([v, c, d]).all()
            assert np.isclose(v - c, d, atol=1e-08, rtol=1e-08)
            ck = (r["subject"], r["trial_index"], f)
            if ck in clean:
                assert np.isclose(c, clean[ck], atol=1e-08, rtol=1e-08)
            clean[ck] = c
        grouped.setdefault((r["subject"], setting), []).append(r)
    assert all((len(v) == 1 for v in index_map.values()))
    settings = sorted({key[1] for key in grouped})
    S = len(settings)
    values = np.full((len(models), S, 3), np.nan)
    trial_counts = np.zeros((len(models), S), int)
    for i, m in enumerate(models):
        for j, s in enumerate(settings):
            rr = grouped.get((m, s), [])
            if not rr:
                continue
            trial_counts[i, j] = len({r["trial_index"] for r in rr})
            for k, f in enumerate(fields):
                values[i, j, k] = model_setting_mean(rr, f)
    valid = np.isfinite(values).all(axis=-1)
    mi, si = np.where(valid)
    method = np.array([settings[j][0] for j in si])
    model = np.array(models)[mi]
    curves = []
    for k, label in [(1, "PQ"), (2, "PS")]:
        rr = curve_model(
            values[mi, si, k],
            values[mi, si, 0],
            method,
            model,
            cfg["seed"] + 500 + k,
            cfg["regression_bootstrap"],
        )
        rr["metric"] = label
        curves.append(rr)
    for r, q in zip(curves, holm([r["p"] for r in curves])):
        r["p_holm"] = float(q)
    dump(
        out / "summary.json",
        {
            "scaling": scaling,
            "curvature": curves,
            "coverage": {
                "models": len(models),
                "settings": S,
                "rows": len(rows),
                "missing_model_settings": int((~valid).sum()),
            },
        },
    )
    inputs.finish(out)
