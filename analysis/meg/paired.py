"""Selected manuscript MEG paired contrasts and sign-flip inference."""

import numpy as np
import pandas as pd

CONTRASTS = (("2-back", "0-back"), ("2-back", "fixation"), ("0-back", "fixation"))
METRICS = ("PQ", "PS")


def contrast_table(participant: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for branch in sorted(participant.branch.unique()):
        subset = participant[participant.branch == branch]
        for metric in METRICS:
            wide = subset.pivot(index="subject_id", columns="condition", values=metric)
            for first, second in CONTRASTS:
                if first not in wide or second not in wide:
                    continue
                diff = wide[first] - wide[second]
                for subject, value in diff.dropna().items():
                    rows.append(
                        {
                            "subject_id": subject,
                            "branch": branch,
                            "metric": metric,
                            "contrast": f"{first}_minus_{second}",
                            "value": float(value),
                        }
                    )
    return pd.DataFrame(rows)


def bh(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    order = np.argsort(values)
    ranked = np.empty_like(values)
    ranked[order] = np.minimum.accumulate(
        (values[order] * len(values) / np.arange(1, len(values) + 1))[::-1]
    )[::-1]
    return np.minimum(ranked, 1.0)


def infer(contrasts: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    rng = np.random.default_rng(int(cfg["inference"]["seed"]))
    rows = []
    primary = contrasts[contrasts.branch == "exact_project1"]
    for (metric, contrast), group in primary.groupby(["metric", "contrast"]):
        values = group.value.to_numpy(float)
        n = len(values)
        mean = float(values.mean())
        sd = float(values.std(ddof=1))
        signs = rng.choice(
            (-1.0, 1.0), size=(int(cfg["inference"]["signflip_permutations"]), n)
        )
        null = (signs * values).mean(axis=1)
        p = (1 + np.count_nonzero(np.abs(null) >= abs(mean))) / (len(null) + 1)
        boot = rng.choice(
            values, size=(int(cfg["inference"]["bootstrap_samples"]), n), replace=True
        ).mean(axis=1)
        rows.append(
            {
                "metric": metric,
                "contrast": contrast,
                "n_participants": n,
                "mean_difference": mean,
                "sd_difference": sd,
                "cohens_dz": mean / sd if sd > 0 else np.nan,
                "ci95_low": float(np.quantile(boot, 0.025)),
                "ci95_high": float(np.quantile(boot, 0.975)),
                "p_signflip": p,
            }
        )
    out = pd.DataFrame(rows)
    out["q_bh_6_tests"] = bh(out.p_signflip.to_numpy(float))
    return out
