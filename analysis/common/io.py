"""Explicit paths, immutable inputs and fresh output directories."""

from pathlib import Path
import hashlib
import json
import platform
import numpy as np
import pandas as pd


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def plain(value):
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, float) and (not np.isfinite(value)):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def dump(path, value):
    Path(path).write_text(
        json.dumps(plain(value), indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def resolve_path(config, value):
    path = Path(value).expanduser()
    return (
        path.resolve()
        if path.is_absolute()
        else (Path(config["_base"]) / path).resolve()
    )


def path(config, key):
    return resolve_path(config, config[key])


def load_config(filename):
    filename = Path(filename).resolve()
    config = json.loads(filename.read_text(encoding="utf-8-sig"))
    config["_base"] = str(filename.parent)
    config["_config"] = str(filename)
    config["_config_sha256"] = sha(filename)
    return config


def output_directory(config):
    out = path(config, "output_dir")
    checkout = Path(__file__).resolve().parents[2]
    if out.is_relative_to(checkout):
        raise ValueError(
            "Output directories must be outside the public code checkout: " + str(out)
        )
    if out.exists():
        raise FileExistsError(
            "Choose a new output_dir; existing runs are never overwritten: " + str(out)
        )
    out.mkdir(parents=True)
    return out


class Inputs:

    def __init__(self, config):
        self.config = config
        self.records = {config["_config"]: config["_config_sha256"]}
        self.code_root = Path(__file__).resolve().parents[1]
        self.code = {
            str(p.relative_to(self.code_root)): sha(p)
            for p in self.code_root.rglob("*.py")
        }

    def add(self, filename, expected=None):
        filename = Path(filename).resolve()
        digest = sha(filename)
        expected = expected or self.config.get("input_sha256", {}).get(str(filename))
        if expected and digest != expected:
            raise ValueError("Input checksum mismatch: " + str(filename))
        self.records[str(filename)] = digest
        return filename

    def finish(self, out):
        for filename, digest in self.records.items():
            if sha(filename) != digest:
                raise RuntimeError("Input changed during analysis: " + filename)
        for filename, digest in self.code.items():
            if sha(self.code_root / filename) != digest:
                raise RuntimeError("Analysis code changed during the run: " + filename)
        dump(
            out / "provenance.json",
            {
                "inputs": self.records,
                "configuration": self.config,
                "python": platform.python_version(),
                "numpy": np.__version__,
                "config_sha256": self.config["_config_sha256"],
                "code_sha256": self.code,
            },
        )


def family_labels(filename, subjects, require_all=True):
    table = pd.read_csv(filename, dtype={"Subject": str, "Family_ID": str})
    table["Subject"] = table.Subject.map(lambda value: str(int(float(value))))
    if table.Subject.duplicated().any():
        raise ValueError("Duplicate subjects in family metadata")
    family = table.set_index("Subject").Family_ID.reindex(
        [str(int(float(s))) for s in subjects]
    )
    missing = family.isna() | family.fillna("").str.strip().eq("")
    if require_all and missing.any():
        raise ValueError("Missing family metadata for analysis participants")
    return (family.fillna("").to_numpy(), (~missing).to_numpy())
