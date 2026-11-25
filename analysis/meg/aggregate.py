"""Equal-weight scan aggregation of the adopted BOLD-like block metrics."""

import numpy as np
import pandas as pd
from analysis.common.io import path, resolve_path, output_directory, Inputs, dump

METRICS = ("PQ", "PS")


def canonical(value: object) -> str:
    text = str(value).strip().lower()
    if text in {"fix15sec", "fixation", "rest", "resting fixation"}:
        return "fixation"
    if text in {"0-back", "0back", "0b"}:
        return "0-back"
    if text in {"2-back", "2back", "2b"}:
        return "2-back"
    return text


def participant_table(blocks: pd.DataFrame, branch: str) -> pd.DataFrame:
    scan = blocks.groupby(["subject_id", "scan_id", "condition"], as_index=False)[
        list(METRICS)
    ].mean()
    counts = (
        scan.groupby(["subject_id", "condition"])["scan_id"]
        .nunique()
        .rename("n_scans")
        .reset_index()
    )
    participant = scan.groupby(["subject_id", "condition"], as_index=False)[
        list(METRICS)
    ].mean()
    participant = participant.merge(counts, on=["subject_id", "condition"])
    participant = participant[participant["n_scans"] == 2].copy()
    participant["branch"] = branch
    return participant


def run(config):
    inputs = Inputs(config)
    manifest = pd.read_csv(
        inputs.add(path(config, "block_files_csv")),
        dtype={"subject_id": str, "scan_id": str},
    )
    required_scans = set((str(s) for s in config["scan_ids"]))
    if required_scans != {"6", "7"}:
        raise ValueError("The manuscript BOLD-like analysis requires scans 6 and 7")
    if (
        manifest.duplicated(["subject_id", "scan_id"]).any()
        or manifest.subject_id.nunique() != config["expected_n"]
    ):
        raise ValueError("Duplicate scan files or unexpected participant count")
    if (
        set(manifest.scan_id) != required_scans
        or len(manifest) != 2 * config["expected_n"]
    ):
        raise ValueError("Exactly two scan files are required per participant")
    blocks = []
    for row in manifest.itertuples():
        table = pd.read_csv(
            inputs.add(resolve_path(config, row.path)),
            dtype={"subject_id": str, "scan_id": str},
        )
        if not {"subject_id", "scan_id", "condition", "PQ", "PS"}.issubset(table):
            raise ValueError("Missing required BOLD-like block columns")
        if set(table.subject_id) != {row.subject_id} or set(table.scan_id) != {
            row.scan_id
        }:
            raise ValueError("Block file identity differs from its manifest row")
        table["condition"] = table.condition.map(canonical)
        table = table[table.condition.isin(["fixation", "0-back", "2-back"])]
        if (
            set(table.condition) != {"fixation", "0-back", "2-back"}
            or not np.isfinite(table[["PQ", "PS"]]).all().all()
        ):
            raise ValueError("Missing conditions or nonfinite BOLD-like block measures")
        blocks.append(table[["subject_id", "scan_id", "condition", "PQ", "PS"]])
    participant = participant_table(
        pd.concat(blocks, ignore_index=True), "exact_project1"
    )
    if (
        len(participant) != 3 * config["expected_n"]
        or participant.subject_id.nunique() != config["expected_n"]
    ):
        raise ValueError("Participant aggregation lost a required scan/condition")
    out = output_directory(config)
    participant.to_csv(out / "condition_participant.csv", index=False)
    dump(
        out / "summary.json",
        {
            "participants": config["expected_n"],
            "scan_aggregation": "equal-weight scan 6/7 after within-scan block mean",
        },
    )
    inputs.finish(out)
