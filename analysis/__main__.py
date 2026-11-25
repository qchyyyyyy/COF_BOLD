"""Run one manuscript analysis using an explicit JSON configuration."""

import argparse
import importlib
from analysis.common.io import load_config

COMMANDS = {
    "task": "analysis.task.inference",
    "behaviour-global": "analysis.behaviour.global_associations",
    "behaviour-features": "analysis.behaviour.prepare_density",
    "behaviour-local": "analysis.behaviour.local_associations",
    "spatial-geometry": "analysis.spatial.geometry",
    "spatial-prepare": "analysis.spatial.prepare_maps",
    "spatial": "analysis.spatial.inference",
    "meg": "analysis.meg.inference",
    "meg-aggregate": "analysis.meg.aggregate",
    "rnn": "analysis.rnn.inference",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis", choices=COMMANDS)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    if not __debug__:
        parser.error(
            "Do not use python -O: scientific input assertions must remain enabled"
        )
    module = importlib.import_module(COMMANDS[args.analysis])
    module.run(load_config(args.config))


if __name__ == "__main__":
    main()
