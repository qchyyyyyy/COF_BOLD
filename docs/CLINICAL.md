# Clinical computation and statistics

`setup('clinical')` exposes a single clinical workflow. Author code uses MIT;
third-party software retains its own terms (see `DEPENDENCIES.md`).

## Input and output contracts

1. `clinical_compute(configFile)` reads an explicit scan list, preprocesses
   CIFTI data and estimates OF/PG/COF fields. Set `COF_TEMPLATE_ROOT` through
   `template_root` and install SPM12.
   Inputs are the original fMRIPrep CIFTIs, not arbitrary raw NIfTI files.
2. `clinical_aggregate_pq(configFile)` reads run directories in the declared
   order and concatenates velocity fields **before** derivative standardization
   and motif counting. It writes participant x time `NumSource`/`NumSink`
   arrays plus subject IDs.
   For ds000115 the order is 0-, 1-, 2-back (475 velocity frames per run).
3. `clinical_run(configFile)` reads the original cohort assets and exports
   `group_effects.csv`, `hc_thresholds.csv`, `quadrant_proportions.csv` and,
   when configured, `disorganization_trend.csv`. Use a fresh absolute output
   folder outside this checkout. Optional plot inputs are private derived data.

Configuration templates are under `configs/clinical`. Use `compute_zib.json`
with original TR 0.8 s or `compute_ds000115.json` with original TR 2.5 s. Copy a template outside
the checkout and replace `REPLACE_...` paths. The collector accepts the original
layout: `dataRoot/SCZ_ZIB/{nback_Events_group.mat,HC_Demographics.csv,
Patients_Demographics.csv}`, `resultsRoot/SCZ_ZIB/subjectList.mat`, and
`resultsRoot/SCZ_ZIB/group/task-nback_dir-PA_run-1_numPattern_global_LR_COF.mat`.
For ds000115 it reads `resultsRoot/SCZ_WM/participants.tsv` and
`resultsRoot/SCZ_WM/group/numPattern_global_LR_COF.mat`.

Historical exclusions are explicit: ZIB row 19 of the original ordered cohort;
ds000115 row 31 **after** selecting SCZ/CON. These are provenance-dependent
row exclusions, not rules for a new cohort. Supply the original row order;
the code does not infer an excluded participant from new data. ZIB arrays must
include the excluded row; ds000115 pattern arrays may already exclude it, as
in the original source. Mismatched sizes, duplicate IDs, unknown sex/group
codes and missing required demographic files fail clearly. ZIB retains the
original two-valued numeric sex coding; textual codes must be M/F, MALE/FEMALE
or 1/0. Missing covariate values are excluded from adjusted models only. Pattern-file subject IDs,
when present, must exactly match this order. The original arrays without IDs
require the explicit `allow_historical_positional_input` setting; their order
cannot be independently verified from values alone.

## Preserved calculations

- Velocity parameters are the effective original call values: `rho=5e-4`,
  `tau=1e-8`, `sigma=8`. The old callers named their first two local variables
  in the opposite order; the configuration uses the estimator's argument names.
- Bilateral velocity masks require all four corners of each grid cell to be
  inside the cortical mask. PG retains the original 2 x 2 x 2 valid average.
  The three ds000115 runs are concatenated before pooled standardization;
  standardizing each run separately would change the original PQ estimator.
  Derivatives use sigma 3, pooled spatiotemporal z-scores, threshold +/-2 and
  the original minimum motif size of 10 grid points.
- ZIB condition means use the original event labels without additional edge
  trimming. ds000115 retains original event windows 38:184 and 250:396 within
  each 475-frame run; accuracy averages target/non-target scores with missing
  values omitted.
- Fig. 5 group models use ordinary OLS with age/sex in ZIB and
  age/education/sex in ds000115. The joint score is `z(accuracy) + z(PQ)`;
  it is not a product. The original BH family contains 3 metrics x all loads
  within each dataset (6 or 9 tests). Raw descriptive estimates and adjusted
  models have separate valid-sample counts.
- HC tenth percentiles use participants with both measures finite. Low means
  strictly below the threshold. Missing pairs are explicitly unassigned.
- Fig. 5g uses the original frozen subject table and quadrant ordering, with
  only `z_dis_4grp`. Outcome and ordinal rank are sample-standardized; the
  original two-sided ordinary OLS trend adjusts age/education/sex. The frozen
  labels are checked and mapped to neutral quadrant labels. They are not
  regenerated from the current illustrative HC-reference panel.
- `clinical_load_comparison(pq,events,group)` returns the descriptive means
  and SEM used for ED12, called separately for OF, PG and COF. Inputs must
  already share the retained participant order. No additional test is fitted.

The broad phenotype scans, regime omnibus/targeted tests, normative residual
correlations, HC3 fits and exploratory clinical exports are absent. Current
Methods prose still mentions some of these wider analyses; repository scope
follows the adopted Results and figure outputs, as requested by the authors.
Historical participant inputs are restricted and are not redistributed here.
Archived input values remain unchanged. Fresh computation fails on missing
scans instead of silently creating zero-filled participant rows; it does not
alter the historical cohort or repair its cached values automatically.
