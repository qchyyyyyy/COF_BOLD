"""Separate D/C linear maps for raw accuracy and RT; Holm over all eight maps."""

import hashlib
import numpy as np
import pandas as pd
from analysis.common.io import path, output_directory, Inputs, dump, family_labels
from analysis.common.multiple_testing import holm
from .local_regression import fit_maps


def run(config):
    inputs = Inputs(config)
    data = pd.read_csv(inputs.add(path(config, "participants_csv")))
    family, _ = family_labels(inputs.add(path(config, "family_csv")), data.Subject)
    with np.load(
        inputs.add(path(config, "density_npz")), allow_pickle=False
    ) as archive:
        if (
            len(data) != config["expected_n"]
            or not data.Subject.is_unique
            or (not np.array_equal(archive["Subject"], data.Subject))
        ):
            raise ValueError("Density and participant identity/order mismatch")
        mask = archive["mask"].astype(bool)
        if mask.sum() != config["expected_locations"]:
            raise ValueError("Unexpected local-map support")
        if (
            config.get("mask_sha256")
            and hashlib.sha256(mask.tobytes()).hexdigest() != config["mask_sha256"]
        ):
            raise ValueError("Local-map mask checksum mismatch")
        if len(np.unique(family)) != config["expected_families"]:
            raise ValueError("Unexpected family count")
        Z = np.column_stack(
            [np.ones(len(data))]
            + [data[k].to_numpy(float) for k in ["age", "sex", "fd_scan"]]
        )
        fitted = {}
        for condition in [0, 2]:
            for component, source in [("D", "Source"), ("C", "Sink")]:
                X = archive[f"{source}Map_{condition}B"]
                expected = config.get("feature_sha256", {}).get(
                    f"{source}Map_{condition}B"
                )
                if expected and hashlib.sha256(X.tobytes()).hexdigest() != expected:
                    raise ValueError("Local-map feature checksum mismatch")
                for outcome in ["accuracy", "rt"]:
                    fitted[f"{condition}_{component}_{outcome}"] = fit_maps(
                        X, data[f"{outcome}_{condition}"].to_numpy(), Z, family
                    )
    all_p = np.concatenate([r["p"] for r in fitted.values()])
    corrected = holm(all_p)
    payload, summary, offset = ({"mask": mask}, [], 0)
    for name, result in fitted.items():
        count = len(result["p"])
        result["p_holm"] = corrected[offset : offset + count]
        offset += count
        for key, value in result.items():
            if isinstance(value, np.ndarray):
                grid = np.full(mask.shape, np.nan)
                grid[mask] = value
                payload[f"raw_accuracy_{name}_{key}"] = grid
        summary.append(
            dict(
                map=name,
                n=len(data),
                families=len(np.unique(family)),
                valid_locations=int(result["valid"].sum()),
                family_tests=int(np.isfinite(all_p).sum()),
                nominal=int((result["p"] < 0.05).sum()),
                Holm=int((result["p_holm"] < 0.05).sum()),
            )
        )
    out = output_directory(config)
    np.savez_compressed(out / "maps.npz", **payload)
    pd.DataFrame(summary).to_csv(out / "summary.csv", index=False)
    dump(out / "summary.json", summary)
    inputs.finish(out)
