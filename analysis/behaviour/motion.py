"""Power framewise displacement: translations in mm, rotations in degrees."""

import numpy as np


def power_fd(motion, radius_mm=50):
    motion = np.asarray(motion, float)
    if (
        motion.ndim != 2
        or motion.shape[1] < 6
        or motion.shape[0] < 2
        or (not np.isfinite(motion[:, :6]).all())
    ):
        raise ValueError("Invalid motion matrix")
    delta = np.diff(motion[:, :6], axis=0)
    return np.abs(delta[:, :3]).sum(1) + radius_mm * np.pi / 180 * np.abs(
        delta[:, 3:6]
    ).sum(1)
