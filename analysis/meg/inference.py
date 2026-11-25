"""Manuscript MEG block tests and unthresholded participant spatial effects."""

import numpy as np
import pandas as pd
import h5py
from scipy import stats
from scipy.io import loadmat
from analysis.common.io import path, resolve_path, output_directory, Inputs, dump
from analysis.common.multiple_testing import holm
from .paired import contrast_table, infer, CONTRASTS, METRICS
from analysis.common.matlab import read_mat, scalar_string


def native_tests(participants, historical):
    rows = []
    for metric in METRICS:
        table = participants.pivot(
            index="subject_id", columns="condition", values=metric + "_mean"
        )
        for first, second in [
            ("0-back", "2-back"),
            ("0-back", "fixation"),
            ("2-back", "fixation"),
        ]:
            difference = (table[first] - table[second]).to_numpy()
            if not np.isfinite(difference).all():
                raise ValueError("Incomplete native MEG participant contrasts")
            contrast = first + "_vs_" + second
            reference = historical[
                (historical.metric == metric) & (historical.contrast == contrast)
            ]
            if len(reference) != 1:
                raise ValueError("Historical bootstrap interval missing or duplicated")
            reference = reference.iloc[0]
            test = stats.wilcoxon(difference, correction=False, method="approx")
            np.testing.assert_allclose(
                test.pvalue, reference.p_raw, rtol=1e-08, atol=1e-10
            )
            np.testing.assert_allclose(
                difference.mean(),
                reference.mean_difference_a_minus_b,
                rtol=0,
                atol=1e-10,
            )
            rows.append(
                dict(
                    metric=metric,
                    contrast=contrast,
                    n=len(difference),
                    mean_difference=difference.mean(),
                    cohens_dz=difference.mean() / difference.std(ddof=1),
                    p=test.pvalue,
                    ci=[reference.bootstrap_ci_low, reference.bootstrap_ci_high],
                    ci_source="original MATLAB participant bootstrap",
                )
            )
    for row, adjusted in zip(rows, holm([r["p"] for r in rows])):
        row["p_holm6"] = adjusted
        reference = historical[
            (historical.metric == row["metric"])
            & (historical.contrast == row["contrast"])
        ].iloc[0]
        np.testing.assert_allclose(adjusted, reference.p_holm, rtol=1e-08, atol=1e-10)
    return rows


def descriptive_curves(table):
    if table.duplicated(["subject_id", "condition", "time_index"]).any():
        raise ValueError("Duplicate participant time bins")
    rows = []
    for metric in METRICS:
        for (condition, index), group in table.groupby(["condition", "time_index"]):
            values = group[metric + "_mean"].to_numpy(float)
            if not np.isfinite(values).all():
                raise ValueError("Nonfinite time-course values")
            rows.append(
                dict(
                    metric=metric,
                    condition=condition,
                    time_index=int(index),
                    mean=values.mean(),
                    sem=values.std(ddof=1) / np.sqrt(len(values)),
                    n=len(values),
                )
            )
    return pd.DataFrame(rows)


def contrast_names(value):
    if isinstance(value, dict):
        return [
            scalar_string(v) for v in np.asarray(value["name"], dtype=object).ravel()
        ]
    return [scalar_string(v["name"]) for v in np.asarray(value, dtype=object).ravel()]


def spatial_effects(files, subject_ids, return_metadata=False):
    resolved = [filename.resolve() for filename in files]
    if len(set(resolved)) != len(files):
        raise ValueError("Duplicate resolved spatial participant paths")
    if len(subject_ids) != len(files) or len(set(subject_ids)) != len(files):
        raise ValueError("Invalid spatial participant identities")
    mean = m2 = reference = None
    for index, (filename, subject) in enumerate(zip(files, subject_ids)):
        participant = (
            read_mat(filename)["participant"]
            if h5py.is_hdf5(filename)
            else loadmat(filename, simplify_cells=True)["participant"]
        )
        if scalar_string(participant["subject_id"]) != str(subject):
            raise ValueError("Spatial participant identity mismatch")
        values = np.asarray(participant["cope"], float)
        if (
            values.ndim != 3
            or values.shape[1:] != (3, 4)
            or (not np.isfinite(values).all())
        ):
            raise ValueError("Invalid participant spatial COPE")
        mask = np.asarray(participant["V_mask"], bool)
        indices = np.asarray(participant["mask_index"]).ravel()
        patterns = [
            scalar_string(v)
            for v in np.asarray(participant["pattern_names"], dtype=object).ravel()
        ]
        contrasts = contrast_names(participant["contrasts"])
        if (
            mask.ndim != 2
            or int(mask.sum()) != values.shape[0]
            or (not np.array_equal(indices, np.flatnonzero(mask.ravel(order="F")) + 1))
        ):
            raise ValueError("Spatial mask_index does not match the MATLAB mask order")
        if len(patterns) != 4 or len(contrasts) != 3:
            raise ValueError("Unexpected spatial contrast/pattern metadata")
        metadata = dict(
            mask=mask,
            mask_index=indices,
            pattern_names=np.array(patterns),
            contrast_names=np.array(contrasts),
            contrast_vector=np.asarray(participant["contrast_vector"]),
        )
        if reference is None:
            reference = metadata
        elif any(
            (
                not np.array_equal(reference[key], value)
                for key, value in metadata.items()
            )
        ):
            raise ValueError("Inconsistent spatial masks, contrasts or pattern order")
        if mean is None:
            mean = np.zeros_like(values, dtype=float)
            m2 = np.zeros_like(values, dtype=float)
        delta = values - mean
        mean += delta / (index + 1)
        m2 += delta * (values - mean)
    if len(files) < 2:
        raise ValueError("At least two spatial participants required")
    effect = mean / (np.sqrt(m2 / (len(files) - 1)) + np.finfo(float).eps)
    return (effect, reference) if return_metadata else effect


def run(config):
    inputs = Inputs(config)
    native = pd.read_csv(
        inputs.add(path(config, "native_participants_csv")), dtype={"subject_id": str}
    )
    historical = pd.read_csv(inputs.add(path(config, "native_bootstrap_csv")))
    bridge = pd.read_csv(
        inputs.add(path(config, "boldlike_participants_csv")), dtype={"subject_id": str}
    )
    if native.subject_id.nunique() != config["expected_n"] or set(
        native.subject_id
    ) != set(bridge.subject_id):
        raise ValueError("Native/BOLD-like participant cohort mismatch")
    if (
        len(native) != 3 * config["expected_n"]
        or len(bridge) != 3 * config["expected_n"]
    ):
        raise ValueError("Expected exactly three conditions per participant")
    if (
        not (native[["n_scans_PQ", "n_scans_PS"]] == 2).all().all()
        or not (bridge.n_scans == 2).all()
    ):
        raise ValueError("Expected two scans per participant")
    if set(bridge.branch) != {"exact_project1"}:
        raise ValueError("Unexpected BOLD-like calculation branch")
    if set(native.condition) != {"fixation", "0-back", "2-back"} or set(
        bridge.condition
    ) != {"fixation", "0-back", "2-back"}:
        raise ValueError("Unexpected MEG condition labels")
    if (
        bridge.duplicated(["subject_id", "condition"]).any()
        or not np.isfinite(bridge[["PQ", "PS"]].to_numpy(float)).all()
    ):
        raise ValueError("Incomplete or duplicate BOLD-like participant conditions")
    native_results = native_tests(native, historical)
    bridge_results = infer(contrast_table(bridge), config)
    if len(bridge_results) != 6:
        raise ValueError("Expected exactly six BOLD-like comparisons")
    out = output_directory(config)
    dump(
        out / "summary.json",
        dict(native=native_results, boldlike=bridge_results.to_dict("records")),
    )
    if config.get("native_curves_csv"):
        table = pd.read_csv(
            inputs.add(path(config, "native_curves_csv")), dtype={"subject_id": str}
        )
        if set(table.subject_id) != set(native.subject_id):
            raise ValueError("Time curves use a different cohort")
        descriptive_curves(table).to_csv(out / "native_time_curves.csv", index=False)
        fixation = (
            table[table.condition == "fixation"]
            .groupby("subject_id")[[m + "_mean" for m in METRICS]]
            .mean()
        )
        dump(
            out / "native_fixation_reference.json",
            {
                m: dict(
                    mean=fixation[m + "_mean"].mean(),
                    sem=fixation[m + "_mean"].std(ddof=1) / np.sqrt(len(fixation)),
                    n=len(fixation),
                )
                for m in METRICS
            },
        )
    for branch in config.get("spatial", []):
        name = branch["name"]
        if name not in ["native", "boldlike"]:
            raise ValueError("Unrecognised MEG spatial branch")
        manifest = pd.read_csv(
            inputs.add(resolve_path(config, branch["participant_files_csv"])),
            dtype={"subject_id": str},
        )
        if (
            manifest.subject_id.duplicated().any()
            or len(manifest) != branch["expected_n"]
            or (not set(manifest.subject_id).issubset(set(native.subject_id)))
        ):
            raise ValueError("Invalid spatial participant cohort")
        files = [inputs.add(resolve_path(config, p)) for p in manifest.path]
        effect, metadata = spatial_effects(
            files, manifest.subject_id.tolist(), return_metadata=True
        )
        np.savez_compressed(
            out / (name + "_spatial_effects.npz"), cohens_d=effect, **metadata
        )
    inputs.finish(out)
