"""Classic and MATLAB v7.3 arrays, cells and scalar structs."""

import h5py
import numpy as np
from scipy.io import loadmat


def read_mat(filename):
    if not h5py.is_hdf5(filename):
        return {k: v for k, v in loadmat(filename).items() if not k.startswith("__")}
    with h5py.File(filename, "r") as handle:

        def decode(item):
            if isinstance(item, h5py.Group):
                return {k: decode(v) for k, v in item.items()}
            value = np.asarray(item)
            if h5py.check_dtype(ref=item.dtype) is not None:
                result = np.empty(value.shape, dtype=object)
                for index in np.ndindex(value.shape):
                    result[index] = (
                        decode(handle[value[index]]) if value[index] else None
                    )
                return result.transpose(tuple(range(result.ndim))[::-1])
            value = value.transpose(tuple(range(value.ndim))[::-1])
            kind = item.attrs.get("MATLAB_class", b"")
            if kind == b"char":
                return np.array(
                    [
                        "".join((chr(int(c)) for c in row))
                        for row in np.atleast_2d(value)
                    ]
                )
            return value

        return {k: decode(v) for k, v in handle.items() if not k.startswith("#")}


def readvar(filename, key):
    if h5py.is_hdf5(filename):
        with h5py.File(filename, "r") as handle:
            item = handle[key]
            return np.asarray(item).transpose(tuple(range(item.ndim))[::-1])
    return loadmat(filename, variable_names=[key])[key]


def scalar_string(value):
    array = np.asarray(value)
    if array.dtype.kind in "US":
        return "".join(array.ravel().tolist())
    if array.size == 1 and array.dtype.kind == "O":
        return scalar_string(array.ravel()[0])
    return str(array.ravel()[0])


def subject_list(filename):
    values = read_mat(filename)
    found = [values[key] for key in ("subjectList", "SubjectList") if key in values]
    if not found:
        raise ValueError("Missing subjectList/SubjectList field")
    if len(found) == 2:
        first = [scalar_string(v) for v in found[0].ravel()]
        second = [scalar_string(v) for v in found[1].ravel()]
        if first != second:
            raise ValueError("Conflicting subjectList and SubjectList fields")
    return found[0]
