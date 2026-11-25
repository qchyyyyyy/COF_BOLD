# Manuscript-to-code map

This map follows the manuscript saved on 25 November 2026. Paths describe
analysis purpose; historical figure numbers do not determine code names.

| Paper component | Computation | Statistics |
| --- | --- | --- |
| Fig. 1: simulations, resting motifs, principal mode | `matlab/workflows/simulations`, `matlab/workflows/rest` | Historical motif summaries; `analysis/spatial` for adopted hierarchy alignment |
| Fig. 2: task reconfiguration and regional/spatial associations | `matlab/workflows/fmri` | `analysis/task`, `analysis/spatial` |
| Fig. 3: native MEG and BOLD-like bridge | `matlab/workflows/meg` | `analysis/meg` |
| Fig. 4: global/local behaviour and descriptive operating groups | `matlab/workflows/behaviour` | `analysis/behaviour`; original descriptive group summaries |
| Fig. 5; ED11–12: clinical comparisons, quadrants, disorganization, method comparison | `matlab/workflows/clinical` | Ordinary OLS and original group correction; Fig. 5g disorganization trend only |
| Fig. 6: RNN scaling and perturbations | `python/rnn`, `matlab/workflows/rnn` | `analysis/rnn` |
| ED4: surrogate comparison | `matlab/workflows/rest`, `matlab/workflows/extended_data` | Original retained tests and summaries |
| ED7: historical method comparisons | `matlab/workflows/extended_data` | Original retained comparisons, with frozen-input requirements |
| ED9: between-load spatial contrasts | `matlab/workflows/extended_data` | Original paired participant tests; no new family adjustment |
| ED10: motif size/order parameter and behaviour | `matlab/workflows/behaviour` | Original adopted models |

The current manuscript still retains the original ED4a pooled F-test values;
the corresponding historical code is retained. No replacement result or new
estimator is introduced during packaging.

Excluded analyses include withdrawn decoding, clinical HC3 alternatives,
clinical broad phenotype/omnibus/residual scans, global behaviour FD diagnostic
correlations, local quadratic/Box–Cox sensitivities, unadopted MEG time-bin tests
and comparisons between old and revised significance tables. No empirical
result files are included. Input arrays needed for historical results remain
external, with their limitations documented in `VERIFICATION.md`.
