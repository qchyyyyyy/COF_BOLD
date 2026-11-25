"""Synthetic numerical checks independent of the released cohort results."""

import importlib
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import stats
from scipy.io import savemat
from scipy.spatial.transform import Rotation

from analysis.common.family import analyze, bh
from analysis.common.regression import ols
from analysis.common.multiple_testing import holm
from analysis.common.matlab import read_mat, scalar_string, subject_list
from analysis.common.io import load_config, output_directory, Inputs
from analysis.behaviour.motion import power_fd
from analysis.behaviour.local_regression import fit_maps
from analysis.rnn.models import curve_model, model_setting_mean
from analysis.spatial.operators import spin_index, corrs
from analysis.spatial.inference import spin_tests
from analysis.meg.paired import infer, contrast_table
from analysis.meg.inference import spatial_effects, descriptive_curves


def test_family_cr1_and_whole_family_bootstrap():
    d = np.array([1.0, 3.0, -2.0, 2.0, 4.0, 6.0, -1.0, 5.0])
    family = np.array([0, 0, 1, 2, 2, 2, 3, 3])
    result = analyze(d, family, 51, 91)
    mean = d.mean()
    sums = np.array([(d[family == g] - mean).sum() for g in range(4)])
    se = np.sqrt(4 / 3 * (sums @ sums) / len(d) ** 2)
    np.testing.assert_allclose(result["se"], se, rtol=1e-14)
    np.testing.assert_allclose(
        result["p"], 2 * stats.t.sf(abs(mean / se), 3), rtol=1e-14
    )
    sampled = np.random.default_rng(51).integers(0, 4, (91, 4))
    draws = [np.concatenate([d[family == g] for g in row]) for row in sampled]
    np.testing.assert_allclose(
        [result["ci_low"], result["ci_high"]],
        np.percentile([v.mean() for v in draws], [2.5, 97.5]),
        atol=1e-14,
    )
    finite_dz = [v.mean() / v.std(ddof=1) for v in draws if v.std(ddof=1) > 0]
    np.testing.assert_allclose(
        [result["dz_ci_low"], result["dz_ci_high"]],
        np.percentile(finite_dz, [2.5, 97.5]),
        atol=1e-14,
    )
    assert analyze(np.zeros(8), family, 1, 5)["p"] == 1
    with pytest.raises(ValueError):
        analyze(d, np.array([""] * 8), 1, 5)


def test_cr1_against_explicit_sandwich_and_local_partial_regression():
    rng = np.random.default_rng(12)
    n = 80
    groups = np.repeat(np.arange(20), 4)
    covariates = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
    predictors = rng.normal(size=(n, 7))
    predictors[:, 6] = 1
    outcome = predictors[:, 0] * 0.4 + rng.normal(size=n) + rng.normal(size=20)[groups]
    maps = fit_maps(predictors, outcome, covariates, groups)
    assert not maps["valid"][6]
    for column in range(6):
        X = np.column_stack([covariates, predictors[:, column]])
        fitted = ols(X, outcome, groups)
        beta = np.linalg.lstsq(X, outcome, rcond=None)[0]
        inverse = np.linalg.inv(X.T @ X)
        residual = outcome - X @ beta
        meat = sum(
            np.outer(
                X[groups == g].T @ residual[groups == g],
                X[groups == g].T @ residual[groups == g],
            )
            for g in range(20)
        )
        reference = inverse @ meat @ inverse * 20 / 19 * 79 / 75
        np.testing.assert_allclose(fitted["cov"], reference, rtol=1e-12, atol=1e-14)
        for key in ["beta", "se", "p"]:
            np.testing.assert_allclose(
                maps[key][column], fitted[key][-1], rtol=1e-10, atol=1e-12
            )
    with pytest.raises(ValueError, match="Rank-deficient"):
        fit_maps(predictors, outcome, np.ones((n, 2)), groups)


def test_multiplicity_and_motion_units():
    values = [0.01, 0.04, 0.03, 0.20]
    np.testing.assert_allclose(
        bh(values), [0.04, 0.05333333333333333, 0.05333333333333333, 0.2]
    )
    np.testing.assert_allclose(holm(values), [0.04, 0.09, 0.09, 0.2])
    motion = np.zeros((3, 12))
    motion[1:, 0] = 1
    motion[2, 3] = 180 / np.pi
    np.testing.assert_allclose(power_fd(motion), [1, 50])


def test_matlab_classic_and_reference_cells(tmp_path):
    classic = tmp_path / "classic.mat"
    savemat(
        classic,
        {
            "condition_names": np.array(["one", "two"], dtype=object),
            "x": np.arange(6).reshape(2, 3),
        },
    )
    hdf = tmp_path / "hdf.mat"
    with h5py.File(hdf, "w") as h:
        h["x"] = np.arange(6).reshape(2, 3).T
        a = h.create_dataset(
            "#refs#/a", data=np.array([[ord(c)] for c in "one"], dtype=np.uint16)
        )
        a.attrs["MATLAB_class"] = np.bytes_("char")
        b = h.create_dataset(
            "#refs#/b", data=np.array([[ord(c)] for c in "two"], dtype=np.uint16)
        )
        b.attrs["MATLAB_class"] = np.bytes_("char")
        cells = h.create_dataset("condition_names", (2, 1), dtype=h5py.ref_dtype)
        cells[0, 0], cells[1, 0] = a.ref, b.ref
    left, right = read_mat(classic), read_mat(hdf)
    np.testing.assert_array_equal(left["x"], right["x"])
    assert [scalar_string(v) for v in right["condition_names"].ravel()] == [
        "one",
        "two",
    ]


def test_subject_list_aliases(tmp_path):
    savemat(tmp_path / "lower.mat", {"subjectList": np.array([1, 2, 3])})
    np.testing.assert_array_equal(
        subject_list(tmp_path / "lower.mat").ravel(), [1, 2, 3]
    )
    savemat(
        tmp_path / "ambiguous.mat",
        {"subjectList": np.array([1, 2, 3]), "SubjectList": np.array([3, 2, 1])},
    )
    with pytest.raises(ValueError, match="Conflicting"):
        subject_list(tmp_path / "ambiguous.mat")


def test_model_equal_weighted_curvature():
    rng = np.random.default_rng(24)
    model = np.repeat(np.arange(30), 8)
    method = np.tile(np.repeat(["a", "b"], 4), 30)
    x = rng.normal(size=len(model))
    y = 0.2 * x * x + 0.1 * x + 0.3 * (method == "b") + rng.normal(size=len(model))
    fitted = curve_model(x, y, method, model, 53, 50)
    weights = 1 / np.bincount(model)[model]
    z = (x - np.average(x, weights=weights)) / np.sqrt(
        np.average((x - np.average(x, weights=weights)) ** 2, weights=weights)
    )
    X = np.column_stack([np.ones(len(x)), z, z * z, method == "b"])
    beta = np.linalg.solve(X.T @ (weights[:, None] * X), X.T @ (weights * y))
    np.testing.assert_allclose(
        fitted["quadratic_beta_standardized_x"], beta[2], rtol=1e-13
    )
    assert (
        model_setting_mean(
            [
                {"trial_index": 0, "v": 1},
                {"trial_index": 0, "v": 3},
                {"trial_index": 1, "v": 6},
            ],
            "v",
        )
        == 4
    )


def test_spin_null_and_hemisphere_assignment():
    rng = np.random.default_rng(75)
    centers = rng.normal(size=(16, 3))
    centers /= np.linalg.norm(centers, axis=1, keepdims=True)
    hemisphere = np.repeat([0, 1], 8)
    regional = rng.normal(size=(16, 3, 4))
    bold = rng.normal(size=(16, 3))
    modes = rng.normal(size=(16, 2))
    gradient = rng.normal(size=16)
    result = spin_tests(regional, bold, modes, gradient, centers, hemisphere, 11, 19)
    null = []
    for rotation in Rotation.random(19, random_state=11).as_matrix():
        ix, _ = spin_index(centers, hemisphere, rotation)
        np.testing.assert_array_equal(hemisphere[ix], hemisphere)
        assert len(set(ix)) == 16
        null.append(
            [
                stats.pearsonr(regional[:, c, m], bold[ix, c]).statistic
                for c in range(3)
                for m in range(4)
            ]
        )
    expected = (1 + (np.abs(null) >= abs(result["observed"].ravel())).sum(0)) / 20
    np.testing.assert_array_equal(result["p"].ravel(), expected)


def test_meg_random_stream_and_spatial_effects(tmp_path):
    rows = []
    rng = np.random.default_rng(81)
    for subject in range(12):
        for condition in ["fixation", "0-back", "2-back"]:
            rows.append(
                dict(
                    subject_id=str(subject),
                    branch="exact_project1",
                    condition=condition,
                    PQ=rng.normal(),
                    PS=rng.normal(),
                )
            )
    contrasts = contrast_table(pd.DataFrame(rows))
    cfg = {
        "inference": {"seed": 63, "signflip_permutations": 41, "bootstrap_samples": 43}
    }
    result = infer(contrasts, cfg)
    random = np.random.default_rng(63)
    for row, (_, group) in zip(
        result.itertuples(), contrasts.groupby(["metric", "contrast"])
    ):
        values = group.value.to_numpy()
        null = (random.choice([-1.0, 1.0], size=(41, len(values))) * values).mean(1)
        boot = random.choice(values, size=(43, len(values)), replace=True).mean(1)
        assert row.p_signflip == (1 + (abs(null) >= abs(values.mean())).sum()) / 42
        np.testing.assert_allclose(
            [row.ci95_low, row.ci95_high], np.quantile(boot, [0.025, 0.975])
        )
    values = rng.normal(size=(8, 9, 3, 4))
    files = []
    for i, array in enumerate(values):
        filename = tmp_path / f"{i}.mat"
        savemat(
            filename,
            {
                "participant": {
                    "cope": array,
                    "subject_id": str(i),
                    "V_mask": np.ones((3, 3)),
                    "mask_index": np.arange(9) + 1,
                    "pattern_names": np.array(["D", "C", "CC", "CW"], dtype=object),
                    "contrasts": np.array(
                        [{"name": name} for name in ["a", "b", "c"]], dtype=object
                    ),
                    "contrast_vector": np.eye(3),
                }
            },
        )
        files.append(filename)
    np.testing.assert_allclose(
        spatial_effects(files, [str(i) for i in range(8)]),
        values.mean(0) / (values.std(0, ddof=1) + np.finfo(float).eps),
        rtol=1e-13,
        atol=1e-13,
    )


def test_output_preservation_and_input_immutability(tmp_path):
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"output_dir": "result"}))
    loaded = load_config(config)
    out = output_directory(loaded)
    with pytest.raises(FileExistsError):
        output_directory(loaded)
    checkout = Path(__file__).resolve().parents[2]
    with pytest.raises(ValueError, match="outside"):
        output_directory({**loaded, "output_dir": str(checkout / "private_results")})
    inputs = Inputs(loaded)
    filename = tmp_path / "input.txt"
    filename.write_text("one")
    inputs.add(filename)
    filename.write_text("two")
    with pytest.raises(RuntimeError, match="changed"):
        inputs.finish(out)


def test_all_command_modules_import_without_running(tmp_path, monkeypatch):
    from analysis.__main__ import COMMANDS

    monkeypatch.chdir(tmp_path)
    for module in COMMANDS.values():
        importlib.import_module(module)
    assert list(tmp_path.iterdir()) == []


def test_spatial_identity_and_grid_contracts(tmp_path):
    values = np.random.default_rng(99).normal(size=(9, 3, 4))
    base = dict(
        cope=values,
        subject_id="a",
        V_mask=np.ones((3, 3)),
        mask_index=np.arange(9) + 1,
        pattern_names=np.array(["D", "C", "CC", "CW"], dtype=object),
        contrasts=np.array([{"name": v} for v in ["a", "b", "c"]], dtype=object),
        contrast_vector=np.eye(3),
    )
    first, second = tmp_path / "a.mat", tmp_path / "b.mat"
    savemat(first, {"participant": base})
    savemat(second, {"participant": {**base, "subject_id": "b"}})
    with pytest.raises(ValueError, match="Duplicate resolved"):
        spatial_effects([first, first], ["a", "b"])
    with pytest.raises(ValueError, match="identity mismatch"):
        spatial_effects([first, second], ["b", "a"])
    savemat(
        second,
        {"participant": {**base, "subject_id": "b", "mask_index": np.arange(9, 0, -1)}},
    )
    with pytest.raises(ValueError, match="mask_index"):
        spatial_effects([first, second], ["a", "b"])
    savemat(
        second,
        {
            "participant": {
                **base,
                "subject_id": "b",
                "pattern_names": np.array(["C", "D", "CC", "CW"], dtype=object),
            }
        },
    )
    with pytest.raises(ValueError, match="pattern order"):
        spatial_effects([first, second], ["a", "b"])
    hdf = tmp_path / "b_hdf.mat"
    with h5py.File(hdf, "w") as h:
        participant = h.create_group("participant")
        references = h.create_group("#refs#")

        def string(parent, key, value):
            item = parent.create_dataset(
                key, data=np.array([[ord(c)] for c in value], dtype=np.uint16)
            )
            item.attrs["MATLAB_class"] = np.bytes_("char")
            return item

        string(participant, "subject_id", "b")
        participant["cope"] = (values + 0.5).transpose(2, 1, 0)
        participant["V_mask"] = np.ones((3, 3)).T
        participant["mask_index"] = (np.arange(9) + 1)[None, :]
        participant["contrast_vector"] = np.eye(3)
        names = participant.create_dataset(
            "pattern_names", (4, 1), dtype=h5py.ref_dtype
        )
        for i, name in enumerate(["D", "C", "CC", "CW"]):
            names[i, 0] = string(references, f"p{i}", name).ref
        contrasts = participant.create_dataset(
            "contrasts", (3, 1), dtype=h5py.ref_dtype
        )
        for i, name in enumerate(["a", "b", "c"]):
            entry = references.create_group(f"c{i}")
            string(entry, "name", name)
            contrasts[i, 0] = entry.ref
    expected = np.stack([values, values + 0.5])
    np.testing.assert_allclose(
        spatial_effects([first, hdf], ["a", "b"]),
        expected.mean(0) / (expected.std(0, ddof=1) + np.finfo(float).eps),
        rtol=1e-12,
    )


def test_pinned_motion_hash_is_fatal(tmp_path):
    from analysis.behaviour.global_associations import run

    if sys.platform == "win32":
        tmp_path = Path("\\\\?\\" + str(tmp_path))
    pd.DataFrame({"Subject": [1, 2]}).to_csv(tmp_path / "participants.csv", index=False)
    pd.DataFrame({"Subject": [1, 2], "Family_ID": [1, 2]}).to_csv(
        tmp_path / "family.csv", index=False
    )
    filename = (
        tmp_path / "motion/1/MNINonLinear/Results/tfMRI_WM_LR/Movement_Regressors.txt"
    )
    filename.parent.mkdir(parents=True)
    np.savetxt(filename, np.zeros((3, 6)))
    config = tmp_path / "config.json"
    config.write_text(
        json.dumps(
            dict(
                participants_csv="participants.csv",
                family_csv="family.csv",
                motion_root="motion",
                output_dir="result",
                expected_n=2,
                minimum_matched_n=1,
                radius_mm=50,
                feature_frames=2,
                input_sha256={str(filename.resolve()): "0" * 64},
            )
        )
    )
    with pytest.raises(ValueError, match="checksum mismatch"):
        run(load_config(config))
