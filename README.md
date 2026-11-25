# Cortical propagation motifs

Computation and statistical analysis for **Brain-wide reconfiguration of
cortical propagation motifs supports human cognition**.

This repository follows the manuscript saved on 25 November 2026. Retained
historical results use their original estimators; revised statistics are
included where adopted in the paper. Code is organized by scientific function,
independently of the authors' internal experiment folders.

| Directory | Contents |
| --- | --- |
| `matlab/core` | Shared velocity estimation, preprocessing, motif extraction and numerical functions |
| `matlab/workflows` | Simulations, resting fMRI, task fMRI, MEG, behaviour, clinical analyses, RNN and extended-data calculations |
| `analysis` | Python statistical modules, with one command-line interface |
| `python/rnn` | RNN training and hidden-state computation |
| `configs` | Configuration templates grouped by analysis purpose |
| `third_party` | Licensed CIFTI/GIfTI/XML readers, nanconv and optional RainCloudPlots helpers |
| `tests` | Synthetic tests; no empirical result tables or participant data |
| `docs` | Manuscript map, input contracts, dependencies and verification |

## Start here

Python 3.10 or newer:

```sh
python -m venv .venv
# Activate the environment using your shell's activation command.
python -m pip install -r requirements-test.txt
python tools/check_source.py --manifest
python -m pytest tests/python
python -m analysis --help
```

Copy an example configuration outside this checkout, supply external input
paths, and choose a fresh output directory. For example:

```sh
python -m analysis task --config /path/to/private/task.json
```

MATLAB, from this directory:

```matlab
setup('clinical');                     % select one scientific workflow
addpath('tools'); check_dependencies;  % bundled IO and toolbox checks
addpath('tests/matlab'); test_core_clinical; test_clinical_inputs;
% clinical_run('/path/to/private/zib.json');
```

Other workflow names are `simulations`, `rest`, `fmri`, `meg`, `behaviour`,
`rnn` and `extended_data`. `setup` adds shared functions and only the selected
workflow. Required MATLAB toolboxes, SPM and MEG software are listed in
[DEPENDENCIES.md](docs/DEPENDENCIES.md).

See [the manuscript map](docs/MANUSCRIPT_MAP.md), [computation guide](docs/COMPUTATION.md),
[statistics guide](docs/STATISTICS.md) and [clinical guide](docs/CLINICAL.md)
for actual entry points and schemas. Optional RNN training uses
`requirements-rnn.txt`.

## Reproducibility and data

Human imaging data, participant identifiers, phenotype/family tables, weights,
and generated empirical results are not distributed. Obtain datasets and
templates under their original access conditions. Keep private configurations
and derived participant data outside this checkout; see [DATA.md](docs/DATA.md).

Tests cover numerical contracts, original-versus-refactored estimators and
synthetic workflows. They do not constitute a new full-cohort reproduction.
[VERIFICATION.md](docs/VERIFICATION.md) distinguishes tested components from
historical input/version gaps and external dependencies.

`SOURCE_MANIFEST.json` records every distributed file; run the source checker
after extracting the archive. `docs/SOURCE_PROVENANCE.json` records source
snapshots and refactoring changes. Figure assembly and most rendering code
are omitted; optional licensed plotting helpers are available separately.

## Licence and citation

Author code uses [MIT](LICENSE). Bundled third-party files retain their
original licences, including GPL terms for the CIFTI/GIfTI dependencies;
see [third-party notices](third_party/README.md). Use [CITATION.cff](CITATION.cff)
to cite the associated manuscript.

## 中文说明

代码按计算用途重新组织，文件名与入口函数对应。统计模块不依赖内部实验目录。
临床部分仅保留论文采用的组间比较、HC 参考四象限和 disorganization 趋势；
未采用的统计过程、额外表型分析及其结果不在发布版中。
