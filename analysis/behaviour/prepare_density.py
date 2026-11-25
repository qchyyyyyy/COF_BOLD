"""Extract local D/C density matrices in the original participant and grid order."""

import hashlib
import numpy as np
import pandas as pd
from analysis.common.io import path, output_directory, Inputs
from analysis.common.matlab import readvar, scalar_string, subject_list


def run(config):
    inputs = Inputs(config)
    data = pd.read_csv(inputs.add(path(config, "participants_csv")))
    coordinate = inputs.add(path(config, "coordinate_file"))
    parts = []
    for side in ["left", "right"]:
        a = readvar(coordinate, "cii_mask_" + side)
        parts.append(
            (a[:-1, :-1] == 1)
            & (a[1:, :-1] == 1)
            & (a[:-1, 1:] == 1)
            & (a[1:, 1:] == 1)
        )
    mask = np.hstack(parts)
    if int(mask.sum()) != config["expected_locations"]:
        raise ValueError("Unexpected density-grid support")
    if (
        config.get("mask_sha256")
        and hashlib.sha256(mask.tobytes()).hexdigest() != config["mask_sha256"]
    ):
        raise ValueError("Density mask checksum mismatch")
    subjects = subject_list(inputs.add(path(config, "subject_list_file")))
    ids = [int(float(scalar_string(s))) for s in subjects.ravel()]
    lookup = {s: i for i, s in enumerate(ids)}
    if len(lookup) != len(ids) or not data.Subject.is_unique:
        raise ValueError("Duplicate participants")
    index = np.array([lookup[s] for s in data.Subject])
    density_file = inputs.add(path(config, "density_mat"))
    payload = {}
    for key in ["SourceMap_0B", "SinkMap_0B", "SourceMap_2B", "SinkMap_2B"]:
        values = readvar(density_file, key)
        if values.shape[:2] != mask.shape or values.shape[2] != len(ids):
            raise ValueError("Unexpected density-array shape: " + key)
        payload[key] = values[mask][:, index].T.astype(float)
        expected = config.get("feature_sha256", {}).get(key)
        if expected and hashlib.sha256(payload[key].tobytes()).hexdigest() != expected:
            raise ValueError("Density-array checksum mismatch: " + key)
    out = output_directory(config)
    np.savez_compressed(
        out / "density_features.npz",
        mask=mask,
        Subject=data.Subject.to_numpy(),
        **payload
    )
    inputs.finish(out)
