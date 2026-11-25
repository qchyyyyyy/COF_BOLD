"""Hemisphere-constrained bijective atlas spins: Holm 12 and Holm 2."""

import numpy as np
from scipy.spatial.transform import Rotation
from analysis.common.io import path, output_directory, Inputs, dump
from .operators import spin_index, corrs, holm


def spin_tests(regional, bold, modes, gradient, centers, hemisphere, seed, n_spins):
    if (
        n_spins < 1
        or not np.isfinite(regional).all()
        or (not np.isfinite(bold).all())
        or (not np.isfinite(modes).all())
        or (not np.isfinite(gradient).all())
    ):
        raise ValueError("Invalid spatial input")
    observed = np.array([corrs(regional[:, c], bold[:, c]) for c in range(3)])
    observed_modes = corrs(modes, gradient)
    null = np.empty((n_spins, 3, 4))
    null_modes = np.empty((n_spins, 2))
    rotations = Rotation.random(n_spins, random_state=seed).as_matrix()
    costs = []
    for index, rotation in enumerate(rotations):
        permutation, cost = spin_index(centers, hemisphere, rotation)
        assert len(np.unique(permutation)) == len(centers)
        assert np.array_equal(hemisphere[permutation], hemisphere)
        costs.append(cost)
        null[index] = np.array(
            [corrs(regional[:, c], bold[permutation, c]) for c in range(3)]
        )
        null_modes[index] = corrs(modes, gradient[permutation])
    p = (1 + (abs(null) >= abs(observed)).sum(0)) / (n_spins + 1)
    p_modes = (1 + (abs(null_modes) >= abs(observed_modes)).sum(0)) / (n_spins + 1)
    return dict(
        observed=observed,
        p=p,
        holm=holm(p.ravel()).reshape(3, 4),
        observed_modes=observed_modes,
        p_modes=p_modes,
        holm_modes=holm(p_modes),
        assignment_distance_quantiles=np.quantile(costs, [0, 0.5, 1]),
    )


def run(config):
    inputs = Inputs(config)
    with np.load(
        inputs.add(path(config, "prepared_maps_npz")), allow_pickle=False
    ) as maps:
        result = spin_tests(
            maps["regional"],
            maps["bold"],
            maps["fig1"],
            maps["gradient"],
            maps["centers"],
            maps["hemi"],
            config["seed"],
            config["n_spins"],
        )
    out = output_directory(config)
    rows = [
        dict(
            contrast=c,
            metric=name,
            r=result["observed"][c, m],
            p_spin=result["p"][c, m],
            p_holm12=result["holm"][c, m],
        )
        for c in range(3)
        for m, name in enumerate(["D", "C", "D+C", "OP"])
    ]
    modes = [
        dict(
            definition=name,
            r_signed=result["observed_modes"][i],
            r_absolute=abs(result["observed_modes"][i]),
            p_spin=result["p_modes"][i],
            p_holm2=result["holm_modes"][i],
        )
        for i, name in enumerate(["continuous_divergence", "signed_two_sd_tails"])
    ]
    dump(
        out / "summary.json",
        dict(
            fig2=rows,
            fig1=modes,
            n_spins=config["n_spins"],
            seed=config["seed"],
            assignment_distance_quantiles=result["assignment_distance_quantiles"],
        ),
    )
    inputs.finish(out)
