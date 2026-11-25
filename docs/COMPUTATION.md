# Computation workflows

Run MATLAB from the repository root and select one domain with `setup('rest')`,
`setup('fmri')`, `setup('meg')`, `setup('behaviour')`, `setup('rnn')`, or
`setup('extended_data')`. This adds the shared core and the selected workflow.
Workflow functions do not reset the MATLAB path. Numerical variants used by MEG
and the historical surrogate calculation have distinct `meg_` and `surrogate_`
names, so they cannot shadow the shared fMRI solver.

Copy the relevant example JSON to a local configuration and supply external
input/output paths. Template images, subject manifests, imaging arrays and
phenotypes remain outside Git. Use fresh output locations. The source-file
names inside historical MAT inputs are part of their data contract and are
preserved even when the public function has a new descriptive name.

## Rest and simulations

| Stage | Entry point | Input/output contract |
|---|---|---|
| Analytic wave examples | `simulate_wavefields(configFile)` | `configs/simulations/simulation.example.json`; packet, radial, standing or Gaussian wavefields; numerical arrays only |
| HCP preprocessing and velocity | `preprocess_hcp(configFile)` | `configs/rest/preprocess.example.json`; CIFTI scans, flat surfaces and coordinate masks; preprocessed movie and bilateral velocity MAT files |
| Participant and group motifs | `rest_motif_metrics(stage,batch,batchSize,cfg)` | `stage='participant'` then `'group'`; fields documented in `rest_config`; archived subject/run filenames |
| Motif density | `rest_motif_densities(cfg)` | Complete participant/group density maps; average maps and numerical plot data |
| Principal modes | `rest_principal_modes(cfg)` | Velocity fields, masks and the external gradient map; spatial modes, singular values and temporal coefficients |
| Global states | `rest_global_states(cfg)` | Complete group motif/order-parameter arrays; descriptive state summaries |
| Surrogate velocity and metrics | `rest_surrogate_velocity(batch,batchSize,cfg)`, `rest_surrogate_metrics(batch,batchSize,cfg)` | Existing preprocessed movies; numerical surrogate derivatives and summaries |

Select simulations with `setup('simulations')`. Rest configuration accepts a
MATLAB struct; explicit `dataRoot`, `vvfDataDir`, `participantMetricDir`,
`groupMetricDir`, `figureDir` and `figureDataDir` overrides keep all derived
files in the chosen external locations. `COF_DATA_ROOT`, `COF_PROCESSED_ROOT`
and `COF_OUTPUT_ROOT` are optional path defaults. Numerical plot bundles can
be retained without installing an Illustrator or figure-rendering workflow.

The Gaussian simulation's original random seed was not recorded. Its example
configuration fixes a new seed for a reproducible reconstruction. It does not
claim recovery of the exact originally plotted noise realization.

## Task fMRI

1. Run `preprocess_hcp('configs/fmri/preprocess.local.json')` with the task
   preprocessing configuration.
2. Set `HCP_TASK_RAW_ROOT` to the original task-directory tree and
   `COF_DATA_ROOT` to the derived scan tree, then run `fmri_task_events`.
   `fmri_gambling_events` and `fmri_development_events` retain the original
   task-specific event exports. Their input event labels and timing remain
   unchanged.
3. Set `COF_RESULTS_ROOT` to the metadata/group tree containing the participant
   lists and events, and `COF_TASK_OUTPUT_ROOT` to a different fresh output
   tree. Run `fmri_task_metrics('participant',1,1,flags)`, then
   `fmri_task_metrics('group',[],[],flags)`. Load `flags` from
   `configs/fmri/feature_flags.json` to select the manuscript features.
4. `fmri_bold_copes('participant',1,1)` and then
   `fmri_bold_copes('group',[],[])` produce the original BOLD comparison maps.

Metadata files are named `<experiment>_SubList.mat` (`subjectList`) and
`<experiment>_EVs.mat` (`Events`, participant by frame). Scan inputs use
`VelField_lh_<method>.mat` and `VelField_rh_<method>.mat`, each containing
`VelField.Ux/Uy`, together with the preprocessed signal where required.
An enabled group output requires every selected participant file. Shorter
participant rows cannot be silently extended with zeros. Group statistics
for the main manuscript are exposed separately through `python -m analysis`.

Historical extended-data task statistics are intentionally separate from
those revised main-text models:

| Function | Historical calculation |
|---|---|
| `task_condition_statistics(configFile)` | The original trimmed condition comparisons and permutation procedure; explicit branch identifies the retained extended-data panel |
| `task_spatial_statistics(patterns,orderParameter,outputFile)` | Original paired spatial tests, including the historically misleading input filename suffix |
| `activation_correspondence(configFile)` | Original activation/motif-map correspondence |
| `task_activation_comparison(time_r,seed,outputFile)` | Original correlation transformation and permutation comparison; requires the archived `time_r` input |
| `surrogate_historical_summary(configFile)` | Original matched real/surrogate histogram, pooled variance test and density summary; exact solver argument order retained |

The exact upstream construction of `time_r` and the original matched
real/surrogate movie pair is not recovered by these functions. Supply the
archived inputs used for the paper. The pooled surrogate F-test is included
because the current manuscript retains that historical result.

## Behaviour

```matlab
setup('behaviour');
behaviour_export_inputs('configs/behaviour/export.local.json');
behaviour_group_profiles('configs/behaviour/group_profiles.local.json');
behaviour_run_metric_models('configs/behaviour/metric_models.local.json');
```

The group input directory contains `tfMRI_WM_LR_SubList.mat` (`SubjectList`),
`tfMRI_WM_LR_info.mat` (`Stim_Type_subs`, `Acc`, `RT`) and bilateral
`tfMRI_WM_LR_numPattern_global_{lh,rh}_COF.mat` (`NumSource`, `NumSink`).
The size/OP model additionally requires the corresponding `sizePattern`
arrays and `tfMRI_WM_LR_orderParam_global_LR_COF.mat`
(`order_param_global`). The unrestricted phenotype table must contain unique
`Subject` IDs plus `Age` and `Gender`; the group-profile calculation also
requires its six named cognition/stress columns.

The original inclusion rule, block order `[8 1 6 4 5 7 2 3]`, 404-frame
timeline, maximum 38 frames per block and per-block Box-Cox transformation are
preserved. In particular, Box-Cox is applied to the selected block values
before condition averaging. The original MATLAB `boxcox` dependency must be
available; the repository does not replace its fitted transformation with a
different implementation.

The participant export contains `Subject`, `age`, `sex` and source, sink,
accuracy, Box-Cox accuracy and reaction time for the 0-back and 2-back
conditions. It feeds the Python `behaviour-global`, `behaviour-features` and
`behaviour-local` commands. Family IDs and movement regressors are separate
restricted inputs to those commands. The local analysis requires the original
frozen block-voxel MAT cache; its unavailable producer is not reconstructed by
substituting a newer voxelwise calculation.

`behaviour_group_profiles` retains only the historical fitted vertices,
group assignment and eight descriptive age/sex-adjusted trait means/SEMs.
It verifies the four frozen group input hashes by default.
`behaviour_run_metric_models` retains the original ordinary quadratic PS and
OP models. The broad trait search, partial-correlation/VIF analyses,
Cook-distance refits, and BOLD/FC sensitivity branches are omitted.

## MEG

The reproducible starting point is the archived session-LCMV source-power
and native velocity derivatives used for the manuscript. Raw sensor MEG,
anatomical models, FieldTrip and the HCP MEG event reader remain external.
The original raw-to-LCMV run record refers to a commit plus uncommitted files
copied to the cluster; an exact executable source snapshot was not frozen.
The current raw-LCMV research scripts are therefore not presented as the
verified source of the archived manuscript derivatives.

Archived power directories contain:

- `session_manifest.json`: `source_count`, `total_frames`, `valid_frames`,
  `frame_interval_seconds`, and chunk filenames;
- `session_metadata.mat`: the `metadata` structure;
- chunk MAT files: `source_power` (source by frame), `sample_index`,
  `core_time_mask`, `valid_time_mask`. Chunk cores must be non-overlapping.

Native COF directories contain `cof_manifest.json` (`velocity_chunks`),
`cof_preprocessing.mat` (`preprocessing.mask`), and velocity chunks with
`Ux`, `Uy`, `interval_frame_index`, `interval_time_seconds`.
`meg_native_motif_metrics(cofDir,outputDir)` computes session-wide
divergence/curl normalization before motif thresholding. It does not
normalize each chunk separately.

Set `COF_MEG_SUBJECT_MANIFEST` to the external locked participant list.
The JSON examples in `configs/meg` expose raw, processed, results, FieldTrip,
MEG Connectome and coordinate paths. Original derivative directory names
remain configurable input identifiers. They are not public API names.

| Stage | Entry points |
|---|---|
| Native condition summaries | `meg_native_conditions(configFile,subjectId)`; `meg_native_group(configFile)` |
| Native trial curves | `meg_native_trial_curves(configFile,subjectId)`; `meg_native_group_curves(configFile)` |
| Native participant/group contrast maps | `meg_native_contrast_maps(configFile,subjectId)`; `meg_group_contrast_maps(configFile)` |
| BOLD-like observation and manuscript preprocessing | `meg_bold_subject(configFile,subjectId,scanId,'full')` |
| BOLD-like participant/group contrast maps | `meg_bold_contrast_maps(configFile,subjectId)`; `meg_bold_group_maps(configFile)` |

The BOLD bridge preserves the original HRF convolution, valid-frame mask,
720-ms sampling, high-pass filter, global-signal regression, Wiener operation,
low-pass filter and hemisphere-specific interpolation. Only the configured
manuscript COF parameter set is exposed. Parameter tuning, unused alternate
observation maps and exploratory QC correlation exports are omitted.
Native time curves retain participant means/SEMs; the 600 unreported
timepoint tests are not computed. The native group workflow retains the six
original block contrasts and their participant-bootstrap confidence intervals;
the Python MEG analysis independently verifies their p values and Holm
adjustment. The extra Friedman and four trial-level comparisons are omitted.
Group map entrypoints retain the complete-cohort checks.

## RNN

`python/rnn/train_working_memory.py` contains the original training/model
definitions. A new training run takes `--output-dir`, `--seed`,
`--num_models`, `--num_neuron` and `--model`. Original training seeds were
not recorded; a newly seeded run is a reproducible new realization.

```matlab
setup('rnn');
cfg=jsondecode(fileread('configs/rnn/perturbations.local.json'));
rnn_perturbation_sweep(cfg);
```

Trained-model MAT inputs contain recurrent/input/readout weights, biases,
leak parameters, stimulus features, targets, initial hidden states and
`coordinates_2d_FD` (or the documented legacy aliases accepted by the loader).
`rnn_hidden_metrics` preserves the RNN-specific interpolation, COF and motif
formulas. The sweep contains only input noise, input dropout, hidden noise,
recurrent pruning, recurrent noise and global leak perturbations. Their
order and design-row seed construction match the six-family configuration.
Missing models, malformed trials and metric failures stop the run.

```text
python python/rnn/rnn_scaling.py configs/rnn/scaling.local.json
python python/rnn/rnn_perturbation_summary.py configs/rnn/summary.local.json
```

Scaling consumes hash-pinned `[time, trial, condition, model]` MAT arrays.
The summary consumes the archived hash-pinned perturbation table. These
preserve the historical descriptive reductions. Model-aware inference is
performed by the separate Python analysis module. The scaling source matches
the historical trend, but its exact correspondence to the final plotted
values remains an explicit provenance limitation.

## Verification boundary

The release checks syntax, function/file-name agreement, unique public
names, scope exclusions, synthetic reductions and selected numerical parity
against the original sources. Synthetic comparisons cover the behaviour
loader and ordinary PS/OP regression statistics and the MEG HRF, high-pass,
grid normalization and source-to-BOLD calculation. No restricted participant
dataset or full training run is rerun for release packaging. The local test
environment lacks the original `boxcox` function, so optimization of that
transformation and the positive-path participant export remain unexecuted
locally; the original implementation call is retained and fails explicitly
when its dependency is absent.
