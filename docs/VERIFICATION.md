# Verification

The release is checked against the original numerical sources and with
synthetic inputs. Test data are generated in temporary directories.

## Numerical and integration checks

- MATLAB: original-versus-refactored COF, OF, PG, four derivative modes,
  motif counts, grid conversions, GLS AR1, seeded permutations and clinical
  ordinary OLS trend; known phase-gradient/stationary-flow contracts.
- Clinical: strict HC threshold equality, missing-pair handling, joint z-score
  sums, adjusted group OLS, BH family, independent trend covariance, rank/input
  guards, across-run normalization and PG within-run averaging. File-based
  tests check nested JSON run lists and original/retained subject ordering.
- Python: family CR1 and full family-bootstrap streams, global/local models,
  BH/Holm, RNN weighted curvature/bootstrap, spatial permutations and surface
  interpolation order, MEG sign-flip/bootstrap streams and task extraction are
  compared with the frozen sources. Synthetic command-line tests exercise
  explicit configurations and fresh-output guards. All 20 public Python tests
  pass, covering all 10 command-line entry points. Twelve groups of numerical
  comparisons pass at relative tolerance 1e-12 and absolute tolerance 1e-13.
- Dependencies: CIFTI write/read with values, time and hemisphere checks;
  GIfTI compressed surface read/write, XML parsing, missing-helper closure,
  nanconv and optional raincloud plotting, with a clean MATLAB path.
- Static source checks: Python parsing, MATLAB Code Analyzer, unique author
  MATLAB filenames, data/credential patterns, public layout and SHA-256
  inventory. Archive entries are read back and compared with source hashes.
  Six additional engineering checks cover computation naming/configuration,
  RNN reduction and preservation of the original six MEG bootstrap kernels.

The test environment is MATLAB R2025b and Python 3.12. Exact final Python
package versions are recorded in `TEST_ENVIRONMENT.json`; earlier estimator
checks also ran with NumPy 2.3.5. Original cluster versions may differ.
The GitHub workflow runs synthetic Python tests when the repository is
uploaded; it has not run remotely during local preparation.

## Execution boundaries

Full participant imaging processing, PyTorch training and manuscript cohort
reruns were not performed for this packaging task. The local MATLAB installation
lacks the original `boxcox` dependency (Financial Toolbox); behaviour loader
and model parity are checked using the existing raw-accuracy option. The
production code retains the original per-block Box–Cox optimization and
requires that toolbox. SPM and full MEG preprocessing tools remain external.

MEG workflows start from the archived session-LCMV power/COF inputs. The exact
historical raw-data producer included uncommitted files whose complete version
was not frozen; current upstream files are not presented as that old snapshot.
Other retained provenance gaps include the original simulated Gaussian random
seed, the exact ED4 example input identity/seed, the frozen ED7c correlation
producer and the original local behaviour block-density cache producer.
These workflows accept the historical inputs explicitly. They do not invent
replacement inputs or new manuscript results.

Clinical Fig. 5g uses its frozen historical subject table/order. Original
ID-less clinical arrays require explicit positional-input opt-in. New arrays
with IDs are checked exactly. No clinical broad phenotype results are bundled.

Scientific formulas intentionally preserved include the historical motif
tracking rule, paired ED9 tests and pooled ED4 F-test. Numerical equivalence
on synthetic inputs supports the refactor; it does not resolve the historical
input-provenance gaps above.
