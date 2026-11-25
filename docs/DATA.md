# Data inputs and outputs

The repository contains source code, configuration templates and synthetic
tests. It contains no human scans, participant identifiers, family assignments,
clinical variables, trained model weights or empirical result tables.

Obtain HCP, OpenNeuro and local-cohort inputs under their original data-use
agreements. HCP restricted family/phenotype metadata and local clinical data
are required by some analyses. Dataset access does not imply permission to
redistribute participant-level derived values. Surface templates and annotation
files retain their original terms and are external inputs.

The workflow guides specify input dimensions, ordering and fields. Paths in
configuration templates are examples, not downloadable bundled datasets.
Copy configurations to a private external directory; use a fresh output
directory for each run. MATLAB templates are resolved through
`COF_TEMPLATE_ROOT` (or `COF_DATA_ROOT/templates` for compatible workflows).
Outputs such as `participants_with_fd.csv`, density matrices, quadrant
assignments and optional plot-input MAT files remain private working data.

`.gitignore` and `tools/check_source.py` reject common data/binary payloads.
They supplement review and do not determine dataset sharing permissions.
