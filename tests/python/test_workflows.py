"""Execute the public CLI using generated inputs only."""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import pytest
from scipy import stats
from scipy.io import savemat

from analysis.common.multiple_testing import holm

REPO = Path(__file__).resolve().parents[2]


def execute(tmp_path, command, config):
    filename = tmp_path / (command + ".json")
    filename.write_text(json.dumps(config), encoding="utf-8")
    completed = subprocess.run(
        [sys.executable, "-m", "analysis", command, "--config", str(filename)],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    output = tmp_path / config["output_dir"]
    assert (output / "provenance.json").is_file()
    return output


def test_task_cli(tmp_path):
    rng = np.random.default_rng(71)
    n = 30
    t = 60
    folder = tmp_path / "task/tfMRI_WM_LR/group"
    folder.mkdir(parents=True)

    def save(suffix, values):
        savemat(folder / ("tfMRI_WM_LR_" + suffix + ".mat"), values)

    save("SubList", {"subjectList": np.arange(100, 100 + n)})
    events = np.tile(np.repeat([0, 1, 2], 20), (n, 1))
    stimulus = np.tile(np.array([[21, 40, 1], [41, 60, 2]])[:, :, None], (1, 1, n))
    save(
        "EVs",
        dict(
            Events=events,
            Stim_Type_subs=stimulus,
            condition_names=np.array(["a", "b"], dtype=object),
        ),
    )
    fields = {
        "numPattern": ["NumSource", "NumSink", "NumCCVortex", "NumCVortex"],
        "sizePattern": ["SizeSource", "SizeSink", "SizeCCVortex", "SizeCVortex"],
        "orderParam": ["order_param_global"],
    }
    for kind, names in fields.items():
        save(
            kind + "_global_LR_COF",
            {name: rng.normal(size=(n, t)) + 0.1 * events for name in names},
        )
    pd.DataFrame(
        {"Subject": np.arange(100, 100 + n), "Family_ID": np.repeat(np.arange(15), 2)}
    ).to_csv(tmp_path / "families.csv", index=False)
    metrics = [
        "PQ_diverging",
        "PQ_converging",
        "PQ_cc_spiral",
        "PQ_c_spiral",
        "PS_diverging",
        "PS_converging",
        "PS_cc_spiral",
        "PS_c_spiral",
        "OP",
    ]
    config = dict(
        input_root="task",
        family_csv="families.csv",
        output_dir="task_out",
        seed=20261119,
        bootstrap_samples=25,
        expected_tests=27,
        window=dict(
            trim_first=10,
            trim_last=10,
            wm_expected_overlap_runs=0,
            wm_expected_overlap_frames=0,
        ),
        tasks=[
            dict(
                key="WM",
                n=n,
                t=t,
                ev_names=["a", "b"],
                control="0-back",
                task="2-back",
                metrics=metrics,
            )
        ],
    )
    out = execute(tmp_path, "task", config)
    result = pd.read_csv(out / "contrasts.csv")
    assert len(result) == 27 and set(result.family) == {"WM27"}
    assert np.isfinite(result.p).all() and result.n.eq(n).all()


def test_behaviour_three_step_cli(tmp_path):
    if sys.platform == "win32":
        tmp_path = Path("\\\\?\\" + str(tmp_path))
    rng = np.random.default_rng(72)
    n = 40
    ids = np.arange(100, 100 + n)
    data = dict(Subject=ids, age=rng.uniform(20, 40, n), sex=np.tile([0, 1], n // 2))
    for condition in [0, 2]:
        for field in ["source", "sink", "boxcox", "accuracy", "rt"]:
            data[f"{field}_{condition}"] = rng.normal(size=n)
    pd.DataFrame(data).to_csv(tmp_path / "participants.csv", index=False)
    pd.DataFrame(dict(Subject=ids, Family_ID=np.repeat(np.arange(n // 2), 2))).to_csv(
        tmp_path / "family.csv", index=False
    )
    for subject in ids:
        folder = tmp_path / "motion" / str(subject) / "MNINonLinear/Results/tfMRI_WM_LR"
        folder.mkdir(parents=True)
        np.savetxt(folder / "Movement_Regressors.txt", rng.normal(size=(21, 12)))
    out = execute(
        tmp_path,
        "behaviour-global",
        dict(
            participants_csv="participants.csv",
            family_csv="family.csv",
            motion_root="motion",
            output_dir="global_out",
            expected_n=n,
            minimum_matched_n=n,
            radius_mm=50,
            feature_frames=20,
        ),
    )
    assert len(pd.read_csv(out / "associations.csv")) == 8
    assert not (out / "diagnostics.json").exists()
    savemat(
        tmp_path / "coordinate.mat",
        {"cii_mask_left": np.ones((5, 4)), "cii_mask_right": np.ones((5, 4))},
    )
    savemat(tmp_path / "subjects.mat", {"subjectList": ids})
    density = {
        f"{component}Map_{condition}B": rng.normal(size=(4, 6, n))
        for component in ["Source", "Sink"]
        for condition in [0, 2]
    }
    savemat(tmp_path / "density.mat", density)
    execute(
        tmp_path,
        "behaviour-features",
        dict(
            participants_csv="global_out/participants_with_fd.csv",
            coordinate_file="coordinate.mat",
            subject_list_file="subjects.mat",
            density_mat="density.mat",
            expected_locations=24,
            output_dir="features_out",
        ),
    )
    out = execute(
        tmp_path,
        "behaviour-local",
        dict(
            participants_csv="global_out/participants_with_fd.csv",
            family_csv="family.csv",
            density_npz="features_out/density_features.npz",
            output_dir="local_out",
            expected_n=n,
            expected_families=n // 2,
            expected_locations=24,
        ),
    )
    assert len(pd.read_csv(out / "summary.csv")) == 8
    with np.load(out / "maps.npz") as result:
        assert result["raw_accuracy_0_D_accuracy_beta"].shape == (4, 6)


def test_spatial_cli(tmp_path):
    rng = np.random.default_rng(73)
    n = 20
    centers = rng.normal(size=(n, 3))
    centers /= np.linalg.norm(centers, axis=1, keepdims=True)
    np.savez(
        tmp_path / "maps.npz",
        regional=rng.normal(size=(n, 3, 4)),
        bold=rng.normal(size=(n, 3)),
        fig1=rng.normal(size=(n, 2)),
        gradient=rng.normal(size=n),
        centers=centers,
        hemi=np.repeat([0, 1], n // 2),
    )
    out = execute(
        tmp_path,
        "spatial",
        dict(
            prepared_maps_npz="maps.npz",
            seed=20261150,
            n_spins=13,
            output_dir="spatial_out",
        ),
    )
    result = json.loads((out / "summary.json").read_text())
    assert len(result["fig2"]) == 12 and len(result["fig1"]) == 2


def test_meg_cli(tmp_path):
    rng = np.random.default_rng(74)
    n = 20
    native = []
    bridge = []
    times = []
    for subject in range(n):
        for condition in ["fixation", "0-back", "2-back"]:
            pq, ps = rng.normal(size=2)
            native.append(
                dict(
                    subject_id=str(subject),
                    condition=condition,
                    PQ_mean=pq,
                    PS_mean=ps,
                    n_scans_PQ=2,
                    n_scans_PS=2,
                )
            )
            bridge.append(
                dict(
                    subject_id=str(subject),
                    branch="exact_project1",
                    condition=condition,
                    PQ=pq + 0.1,
                    PS=ps + 0.2,
                    n_scans=2,
                )
            )
            for index in range(5):
                times.append(
                    dict(
                        subject_id=str(subject),
                        condition=condition,
                        time_index=index,
                        PQ_mean=rng.normal(),
                        PS_mean=rng.normal(),
                    )
                )
    native = pd.DataFrame(native)
    native.to_csv(tmp_path / "native.csv", index=False)
    pd.DataFrame(bridge).to_csv(tmp_path / "bridge.csv", index=False)
    pd.DataFrame(times).to_csv(tmp_path / "curves.csv", index=False)
    historical = []
    for metric in ["PQ", "PS"]:
        wide = native.pivot(
            index="subject_id", columns="condition", values=metric + "_mean"
        )
        for a, b in [
            ("0-back", "2-back"),
            ("0-back", "fixation"),
            ("2-back", "fixation"),
        ]:
            difference = (wide[a] - wide[b]).to_numpy()
            historical.append(
                dict(
                    metric=metric,
                    contrast=a + "_vs_" + b,
                    p_raw=stats.wilcoxon(
                        difference, correction=False, method="approx"
                    ).pvalue,
                    mean_difference_a_minus_b=difference.mean(),
                    bootstrap_ci_low=difference.mean() - 0.5,
                    bootstrap_ci_high=difference.mean() + 0.5,
                )
            )
    for row, p in zip(historical, holm([r["p_raw"] for r in historical])):
        row["p_holm"] = p
    pd.DataFrame(historical).to_csv(tmp_path / "intervals.csv", index=False)
    files = []
    for subject in range(n):
        filename = f"spatial_{subject}.mat"
        savemat(
            tmp_path / filename,
            {
                "participant": {
                    "cope": rng.normal(size=(12, 3, 4)),
                    "subject_id": str(subject),
                    "V_mask": np.ones((3, 4)),
                    "mask_index": np.arange(12) + 1,
                    "pattern_names": np.array(["D", "C", "CC", "CW"], dtype=object),
                    "contrasts": np.array(
                        [{"name": name} for name in ["a", "b", "c"]], dtype=object
                    ),
                    "contrast_vector": np.eye(3),
                }
            },
        )
        files.append(dict(subject_id=str(subject), path=filename))
    pd.DataFrame(files).to_csv(tmp_path / "spatial.csv", index=False)
    out = execute(
        tmp_path,
        "meg",
        dict(
            native_participants_csv="native.csv",
            native_bootstrap_csv="intervals.csv",
            boldlike_participants_csv="bridge.csv",
            native_curves_csv="curves.csv",
            expected_n=n,
            output_dir="meg_out",
            inference=dict(
                seed=20260910, signflip_permutations=31, bootstrap_samples=33
            ),
            spatial=[
                dict(name="native", expected_n=n, participant_files_csv="spatial.csv")
            ],
        ),
    )
    result = json.loads((out / "summary.json").read_text())
    assert len(result["native"]) == 6 and len(result["boldlike"]) == 6
    assert set(result) == {"native", "boldlike"}
    assert len(pd.read_csv(out / "native_time_curves.csv")) == 30
    assert (out / "native_spatial_effects.npz").exists()


def test_meg_scan_aggregation_cli(tmp_path):
    rows = []
    for subject in ["one", "two"]:
        for scan in ["6", "7"]:
            filename = f"{subject}_{scan}.csv"
            records = []
            for condition in ["fix15sec", "0back", "2back"]:
                count = 2 if scan == "6" else 3
                records.extend(
                    dict(
                        subject_id=subject,
                        scan_id=scan,
                        condition=condition,
                        PQ=int(scan) + repeat,
                        PS=10 + repeat,
                    )
                    for repeat in range(count)
                )
            pd.DataFrame(records).to_csv(tmp_path / filename, index=False)
            rows.append(dict(subject_id=subject, scan_id=scan, path=filename))
    pd.DataFrame(rows).to_csv(tmp_path / "blocks.csv", index=False)
    out = execute(
        tmp_path,
        "meg-aggregate",
        dict(
            block_files_csv="blocks.csv",
            scan_ids=["6", "7"],
            expected_n=2,
            output_dir="aggregate_out",
        ),
    )
    table = pd.read_csv(out / "condition_participant.csv")
    assert len(table) == 6 and table.n_scans.eq(2).all()
    # Equal-weight means of scan means, despite unequal block counts.
    np.testing.assert_allclose(table.PQ, (6.5 + 8) / 2)
    np.testing.assert_allclose(table.PS, (10.5 + 11) / 2)


def test_rnn_cli(tmp_path):
    rng = np.random.default_rng(76)
    scaling = []
    for size in [144, 256, 400, 576, 784, 1024, 1296]:
        for metric in ["PQ", "PS"]:
            prefix = "Num" if metric == "PQ" else "Size"
            filename = f"{size}_{metric}.mat"
            shape = [8, 3, 3, 30]
            savemat(
                tmp_path / filename,
                {
                    prefix + "Source": rng.normal(size=shape),
                    prefix + "Sink": rng.normal(size=shape),
                },
            )
            scaling.append(dict(path=filename, size=size, metric=metric, shape=shape))
    rows = []
    for subject in range(30):
        for method in ["a", "b"]:
            for severity in [0.0, 0.5, 1.0]:
                for trial in range(3):
                    row = dict(
                        subject=f"model_{subject}",
                        subject_num=subject,
                        trial_index=trial,
                        method=method,
                        severity=severity,
                        fraction=0.5,
                        repeat=1,
                    )
                    for field in ["accuracy", "PQ_DC_mean", "PS_DC_mean"]:
                        value = rng.normal() + severity
                        clean = 0.3 + trial * 0.1
                        row.update(
                            {
                                field: value,
                                "clean_" + field: clean,
                                "delta_" + field: value - clean,
                            }
                        )
                    rows.append(row)
    pd.DataFrame(rows).to_csv(tmp_path / "trials.csv", index=False)
    out = execute(
        tmp_path,
        "rnn",
        dict(
            output_dir="rnn_out",
            seed=20261146,
            bootstrap=31,
            regression_bootstrap=31,
            expected_scaling_tests=14,
            scaling={"inputs": scaling},
            input=dict(path="trials.csv", rows=len(rows), unique_subjects=30),
        ),
    )
    result = json.loads((out / "summary.json").read_text())
    assert len(result["scaling"]) == 14 and len(result["curvature"]) == 2


def test_spatial_preparation_cli(tmp_path):
    if sys.platform == "win32":
        tmp_path = Path("\\\\?\\" + str(tmp_path))
    rng = np.random.default_rng(79)
    n = 8
    locations = 45990
    mask = np.zeros((200, 500), bool)
    mask.ravel()[:locations] = True
    flat = np.array(
        [[-498 + 2 * (i % 180), -198 + 2 * (i // 180), 0] for i in range(360)]
    )
    centers = rng.normal(size=(360, 3))
    centers /= np.linalg.norm(centers, axis=1, keepdims=True)
    np.savez(
        tmp_path / "geometry.npz",
        support=np.ones(360, bool),
        flat=flat,
        parcel_ids=np.arange(360),
        labels=np.arange(360),
        area=np.ones(360),
        B=np.r_[np.ones(180), np.zeros(180)][:, None],
        terms=np.array(["working_memory"]),
        centers=centers,
        parcel_hemi=np.repeat([0, 1], 180),
    )
    group = tmp_path / "task/group"
    group.mkdir(parents=True)
    ids = np.arange(100, 100 + n)
    savemat(group / "tfMRI_WM_LR_SubList.mat", {"SubjectList": ids})
    pd.DataFrame(dict(Subject=ids, Family_ID=np.arange(n))).to_csv(
        tmp_path / "family.csv", index=False
    )
    mean = np.zeros((locations, 3, 4))
    m2 = mean.copy()
    for i, subject in enumerate(ids):
        folder = tmp_path / "task/participant" / str(subject)
        folder.mkdir(parents=True)
        density = rng.normal(size=(locations, 3, 4))
        op = rng.normal(size=(locations, 3))
        savemat(
            folder / (str(subject) + "_tfMRI_WM_LR_numPattern_level1_LR_COF.mat"),
            {"cope": density},
        )
        savemat(
            folder / (str(subject) + "_tfMRI_WM_LR_orderParam_level1_LR_COF.mat"),
            {"cope": op},
        )
        xx = np.stack(
            [density[:, :, 0], density[:, :, 1], density[:, :, :2].sum(2), op], axis=2
        )
        delta = xx - mean
        mean += delta / (i + 1)
        m2 += delta * (xx - mean)
    dz = mean / (np.sqrt(m2 / (n - 1)) + np.finfo(float).eps)
    maps = np.zeros((200, 500, 4, 3))
    for c in range(3):
        for k, source in enumerate([0, 1, 3, None]):
            grid = np.zeros(mask.size)
            grid[mask.ravel(order="F")] = (
                rng.normal(size=locations) if source is None else dz[:, c, source]
            )
            maps[:, :, k, c] = grid.reshape(mask.shape, order="F")
    savemat(tmp_path / "effects.mat", dict(mask=mask, maps=maps))
    div = rng.normal(size=mask.shape) * 3
    gradient = rng.normal(size=mask.shape)
    x = div.ravel(order="F")[mask.ravel(order="F")]
    y = gradient.ravel(order="F")[mask.ravel(order="F")]
    tail = np.where(abs(x) >= 2, np.sign(x), 0)
    with h5py.File(tmp_path / "mode.mat", "w") as h:
        mode = h.create_group("mode")
        geo = h.create_group("geometry")
        geo["V_mask"] = mask.T
        name = h.create_dataset("name", data=np.array([ord(c) for c in "COF"]))
        divergence = h.create_dataset("divergence", data=div.T)
        correlation = h.create_dataset(
            "correlation", data=stats.pearsonr(tail, y).statistic
        )
        for key, value in [
            ("method", name),
            ("Div", divergence),
            ("div_r", correlation),
        ]:
            item = mode.create_dataset(key, (1, 1), dtype=h5py.ref_dtype)
            item[0, 0] = value.ref
        grad = h.create_dataset("gradient", data=gradient.T)
        refs = h.create_dataset("args", (4, 1), dtype=h5py.ref_dtype)
        refs[:, 0] = [mode.ref, geo.ref, geo.ref, grad.ref]
    out = execute(
        tmp_path,
        "spatial-prepare",
        dict(
            geometry_npz="geometry.npz",
            task_effects_mat="effects.mat",
            mode_render_mat="mode.mat",
            wm_task_root="task",
            family_csv="family.csv",
            expected_n=n,
            roi="working_memory",
            output_dir="prepared_out",
        ),
    )
    assert len(json.loads((out / "roi_tests.json").read_text())) == 12
    with np.load(out / "prepared_maps.npz") as result:
        assert result["regional"].shape == (360, 3, 4)


@pytest.fixture
def geometry_directory(tmp_path):
    # Nibabel normalises Windows extended-length paths to POSIX slashes.
    # Keep its synthetic file tree short, and remove it when the test ends.
    if sys.platform == "win32" and len(str(tmp_path)) > 140:
        with tempfile.TemporaryDirectory(
            prefix=".geometry-test-", dir=REPO.parent
        ) as folder:
            yield Path(folder)
    else:
        yield tmp_path


def test_gifti_cifti_geometry_cli(geometry_directory):
    tmp_path = geometry_directory
    nib = pytest.importorskip("nibabel")
    from nibabel.cifti2 import cifti2_axes

    n = 32492
    width = 181
    index = np.arange(n)
    row = index // width
    column = index % width
    vertices = np.column_stack([-248 + column, -198 + row, np.zeros(n)]).astype(
        np.float32
    )
    faces = []
    for i in range(n - width - 1):
        if i % width < width - 1:
            faces.extend([[i, i + 1, i + width], [i + 1, i + width + 1, i + width]])
    faces = np.array(faces, np.int32)
    rng = np.random.default_rng(80)
    sphere = rng.normal(size=(n, 3)).astype(np.float32)
    sphere /= np.linalg.norm(sphere, axis=1, keepdims=True)
    templates = tmp_path / "templates"
    base = templates / "HCP_S1200_Atlas_Z4_pkXDZ"
    base.mkdir(parents=True)
    for hemisphere in ["L", "R"]:
        for stem, coords in [
            ("flat", vertices),
            ("sphere", sphere),
            ("midthickness_MSMAll", vertices),
        ]:
            image = nib.gifti.GiftiImage(
                darrays=[
                    nib.gifti.GiftiDataArray(coords, intent="NIFTI_INTENT_POINTSET"),
                    nib.gifti.GiftiDataArray(faces, intent="NIFTI_INTENT_TRIANGLE"),
                ]
            )
            nib.save(image, base / f"S1200.{hemisphere}.{stem}.32k_fs_LR.surf.gii")
    with h5py.File(templates / "coordinate_cii_v2.mat", "w") as h:
        for side in ["left", "right"]:
            mask = np.zeros((201, 251))
            mask[:116, :201] = 1
            mask[115, :5] = 0
            h["cii_mask_" + side] = mask.T
            h["i_invalid_" + side] = np.array([], int)
            h["x_" + side] = vertices[:, 0]
            h["y_" + side] = vertices[:, 1]
    axis = cifti2_axes.BrainModelAxis.from_surface(
        np.arange(n), n, "CortexLeft"
    ) + cifti2_axes.BrainModelAxis.from_surface(np.arange(n), n, "CortexRight")
    labels = np.r_[np.arange(n) % 180 + 1, np.arange(n) % 180 + 181]
    label_axis = cifti2_axes.LabelAxis(
        ["atlas"], [{i: (str(i), (1.0, 0.0, 0.0, 1.0)) for i in range(1, 361)}]
    )
    image = nib.Cifti2Image(
        labels[None, :], header=nib.Cifti2Header.from_axes([label_axis, axis])
    )
    nib.save(
        image,
        base
        / "Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii",
    )
    bat = templates / "BAT"
    bat.mkdir()
    arrays = {
        name: np.empty(2, dtype=object)
        for name in [
            "term_surface_map_res_all",
            "term_surface_map_res_left",
            "term_surface_map_res_right",
        ]
    }
    for j in range(2):
        values = (np.arange(2 * n) % 2 == j).astype(float)
        arrays["term_surface_map_res_all"][j] = np.column_stack(
            [np.arange(2 * n) + 1, values]
        )
        arrays["term_surface_map_res_left"][j] = np.column_stack(
            [np.arange(n) + 1, values[:n]]
        )
        arrays["term_surface_map_res_right"][j] = np.column_stack(
            [np.arange(n) + 1, values[n:]]
        )
    savemat(
        bat / "cifti_database_data.mat",
        dict(
            dir_path=np.array(
                [{"name": "working_memory_pAgF"}, {"name": "control_pAgF"}],
                dtype=object,
            ),
            **arrays,
        ),
    )
    out = execute(
        tmp_path,
        "spatial-geometry",
        dict(templates_root="templates", output_dir="geometry_out"),
    )
    with np.load(out / "geometry.npz") as result:
        assert result["centers"].shape == (360, 3)
        assert result["support"].sum() > n
