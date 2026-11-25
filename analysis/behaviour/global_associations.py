"""Figure 4 global D+C behaviour associations, adjusted for scan FD."""

import numpy as np
import pandas as pd
from scipy import stats
from analysis.common.io import path, output_directory, Inputs, dump, family_labels
from analysis.common.regression import ols, design
from analysis.common.family import bh
from .motion import power_fd


def run(config):
    inputs = Inputs(config)
    data = pd.read_csv(inputs.add(path(config, "participants_csv")))
    if not data.Subject.is_unique or len(data) != config["expected_n"]:
        raise ValueError("Participant identity/count contract failed")
    family, _ = family_labels(inputs.add(path(config, "family_csv")), data.Subject)
    out = output_directory(config)
    fd, motion_status = ([], [])
    for subject in data.Subject:
        filename = (
            path(config, "motion_root")
            / str(int(subject))
            / "MNINonLinear/Results/tfMRI_WM_LR/Movement_Regressors.txt"
        )
        value, status = (np.nan, "missing")
        if filename.is_file():
            inputs.add(filename)
            try:
                motion = np.loadtxt(filename)
                displacement = power_fd(motion, config["radius_mm"])
                if len(displacement) != config["feature_frames"]:
                    raise ValueError("Unexpected number of adjacent motion differences")
                value, status = (displacement.mean(), "ok")
            except (ValueError, OSError) as error:
                status = "invalid: " + str(error)
        fd.append(value)
        motion_status.append(dict(Subject=int(subject), status=status))
    data["fd_scan"] = fd
    data.to_csv(out / "participants_with_fd.csv", index=False)
    good = np.isfinite(data.fd_scan)
    data, family = (data.loc[good].reset_index(drop=True), family[good])
    if len(data) < config["minimum_matched_n"]:
        raise ValueError("Insufficient same-scan motion coverage")
    associations, curves = ([], [])
    fd = data.fd_scan.to_numpy()
    if fd.std(ddof=1) <= 0:
        raise ValueError("Scan FD has zero variance")
    fz = (fd - fd.mean()) / fd.std(ddof=1)
    for condition in [0, 2]:
        pq = (data[f"source_{condition}"] + data[f"sink_{condition}"]).to_numpy()
        mu, sd = (pq.mean(), pq.std(ddof=1))
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError("D+C predictor has zero variance")
        age, sex = (data.age.to_numpy(), data.sex.to_numpy())
        grid = np.linspace(pq.min(), pq.max(), 160)
        for outcome in ["boxcox", "rt"]:
            y = data[f"{outcome}_{condition}"].to_numpy()
            for quadratic in [False, True]:
                X = np.column_stack([design(pq, age, sex, mu, sd, quadratic), fz])
                P = np.column_stack(
                    [
                        design(
                            grid,
                            np.full(160, age.mean()),
                            np.full(160, sex.mean()),
                            mu,
                            sd,
                            quadratic,
                        ),
                        np.zeros(160),
                    ]
                )
                fit = ols(X, y, family)
                index = 2 if quadratic else 1
                critical = stats.t.ppf(0.975, fit["df"])
                row = dict(
                    condition=condition,
                    outcome=outcome,
                    model="quadratic" if quadratic else "linear",
                    n=len(y),
                    families=len(np.unique(family)),
                    beta=fit["beta"][index],
                    ci_low=fit["beta"][index] - critical * fit["se"][index],
                    ci_high=fit["beta"][index] + critical * fit["se"][index],
                    p=fit["p"][index],
                    df=fit["df"],
                )
                associations.append(row)
                prediction = P @ fit["beta"]
                band = critical * np.sqrt(
                    np.maximum(np.einsum("ij,jk,ik->i", P, fit["cov"], P), 0)
                )
                curves.append(
                    dict(
                        condition=condition,
                        outcome=outcome,
                        model=row["model"],
                        x=grid,
                        y=prediction,
                        low=prediction - band,
                        high=prediction + band,
                    )
                )
    for row, q in zip(associations, bh([r["p"] for r in associations])):
        row["q_BH8"] = q
    pd.DataFrame(associations).to_csv(out / "associations.csv", index=False)
    dump(out / "curves.json", curves)
    dump(
        out / "summary.json",
        dict(
            matched_n=len(data),
            excluded_n=int((~good).sum()),
            motion_status=motion_status,
        ),
    )
    inputs.finish(out)
