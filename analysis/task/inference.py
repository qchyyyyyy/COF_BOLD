"""Participant-equal task contrasts with family CR1 and family bootstrap."""

import numpy as np
import pandas as pd
from analysis.common.io import path, output_directory, Inputs, dump, family_labels
from analysis.common.family import analyze, bh
from analysis.common.matlab import read_mat, scalar_string
from .inputs import task_paths, load_task, condition_means


def run(config):
    config = dict(config)
    config["input_root"] = str(path(config, "input_root"))
    inputs = Inputs(config)
    family_file = inputs.add(path(config, "family_csv"))
    out = output_directory(config)
    rows, coverage = ([], [])
    for task_index, task in enumerate(config["tasks"]):
        for filename in task_paths(config, task).values():
            inputs.add(filename)
        metrics, events, valid, _, summary = load_task(config, task)
        subjects = read_mat(task_paths(config, task)["subjects"])
        subjects = next((v for k, v in subjects.items() if k.lower() == "subjectlist"))
        ids = [scalar_string(v) for v in subjects.ravel()]
        family, covered = family_labels(family_file, ids, require_all=False)
        summary["family_missing_n"] = int((~covered).sum())
        coverage.append(summary)
        names = (
            ["0back_minus_fixation", "2back_minus_fixation", "2back_minus_0back"]
            if task["key"] == "WM"
            else ["control_minus_fixation", "task_minus_fixation", "task_minus_control"]
        )
        for metric_index, (metric, values) in enumerate(metrics.items()):
            means, counts = condition_means(values, events, valid)
            for contrast_index, (first, second) in enumerate(((1, 0), (2, 0), (2, 1))):
                good = (
                    covered
                    & np.isfinite(means[:, first])
                    & np.isfinite(means[:, second])
                )
                row = dict(
                    task=task["key"],
                    metric=metric,
                    contrast=names[contrast_index],
                    family="WM27" if task["key"] == "WM" else "OTHER75",
                    available_n=len(ids),
                    family_missing_n=int((~covered).sum()),
                    excluded_n=int((~good).sum()),
                )
                if good.sum() < 3 or len(np.unique(family[good])) < 2:
                    row.update(
                        status="insufficient_pairs_or_families",
                        n=int(good.sum()),
                        p=np.nan,
                    )
                else:
                    seed = int(
                        np.random.SeedSequence(
                            [
                                config["seed"],
                                task_index,
                                0,
                                metric_index,
                                contrast_index,
                            ]
                        ).generate_state(1)[0]
                    )
                    row.update(
                        analyze(
                            means[good, first] - means[good, second],
                            family[good],
                            seed,
                            config["bootstrap_samples"],
                        )
                    )
                    row.update(
                        seed=seed,
                        min_frames_a=int(counts[good, first].min()),
                        min_frames_b=int(counts[good, second].min()),
                    )
                rows.append(row)
    if len(rows) != config["expected_tests"]:
        raise ValueError(
            "Task/metric configuration does not match the specified test family"
        )
    for family_name in sorted({r["family"] for r in rows}):
        selected = [r for r in rows if r["family"] == family_name]
        for row, q in zip(
            selected, bh([r["p"] if np.isfinite(r["p"]) else 1 for r in selected])
        ):
            row["q"] = float(q) if np.isfinite(row["p"]) else np.nan
    pd.DataFrame(rows).to_csv(out / "contrasts.csv", index=False)
    dump(out / "summary.json", {"contrasts": rows, "coverage": coverage})
    inputs.finish(out)
