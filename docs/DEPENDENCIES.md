# Dependencies

## Bundled MATLAB dependencies

Run the repository setup function, or add these directories individually:

```matlab
addpath(fullfile(repo,'third_party','cifti_matlab'));
addpath(fullfile(repo,'third_party','nanconv'));
addpath(fullfile(repo,'third_party','raincloudplots')); % optional graphics
addpath(fullfile(repo,'tools'));
report = check_dependencies();
```

Do not use `genpath` on `cifti_matlab`: MATLAB resolves its `@gifti`, `@xmltree`
and `private` folders through the library root. Use one CIFTI library on the path
at a time. `check_dependencies` rejects shadowed bundled entry points.

The reader used by the computational code is `ft_read_cifti`. For dense fMRI:

```matlab
c = ft_read_cifti(filename,'readdata',true,'readsurface',false);
```

For a `.dtseries.nii` file, `c.dtseries` contains the vertex-by-time array and
`c.brainstructure` identifies its brain structures. The legacy reader expands
surface brain models to full surface vertex counts, using NaNs for omitted
vertices. Keep the code's existing brain-structure and finite-value masks.
`gifti(filename)` reads surface geometry and supports ASCII, base64 and compressed
base64 data. Compressed GIfTI uses the MATLAB JVM; no MEX compilation is required.
The synthetic test verifies numerical values, time points and hemisphere order.
It is a dependency test, not an end-to-end reconstruction of the manuscript.

When using `ft_write_cifti`, pass an extensionless filename; this archived writer
appends the appropriate `.dtseries.nii`, `.dscalar.nii` or related suffix itself.

The retained CIFTI/GIfTI/XMLTree source is from the original computational
archive. Its FieldTrip `VERSION` records
[`d6ada706fdb09f824c38e4de7655757076be4a21`](https://github.com/fieldtrip/fieldtrip/tree/d6ada706fdb09f824c38e4de7655757076be4a21).
The archive already contains source-level changes, so the release records hashes
of the exact files instead of claiming the entire tree is pristine upstream.
The three missing nested-field helpers were obtained at that pinned commit.
All original `GPL`, `LICENSE.fieldtrip`, `LICENSE.gifti` and `LICENSE.xmltree`
notices are retained in `third_party/cifti_matlab`.

## Optional plotting helpers

The two included [RainCloudPlots](https://github.com/RainCloudPlots/RainCloudPlots)
functions are the original-project copies, with the original MIT licence.
Their default `density_type='ks'` uses `ksdensity` and requires MATLAB's
Statistics and Machine Learning Toolbox. The alternate `'rash'` branch requires
the separately installed [Robust Statistical Toolbox](https://github.com/CPernet/Robust_Statistical_Toolbox);
it is not used by the calculations or tested here. Set the MATLAB random seed
before plotting when reproducible jitter is needed.

The archived `daviolinplot`, `tight_subplot` and `bluewhitered` copies do not
carry sufficiently clear redistribution licences in the available source.
They are omitted. No retained computational entry point requires these helpers.
The archived SCZ plotting scripts and their results are not included.

## External toolboxes

| Software | Required by | Setup |
| --- | --- | --- |
| MATLAB Statistics and Machine Learning Toolbox | Statistical functions, clustering, optional raincloud plots | Install through MATLAB; `check_dependencies` tests optional raincloud plotting when available |
| MATLAB Signal Processing Toolbox | Filtering, Hilbert transforms and related time-series processing | Install through MATLAB |
| MATLAB Image Processing Toolbox | Gaussian kernels (`fspecial`) used by OF/PG/COF and velocity derivatives | Install through MATLAB |
| MATLAB Financial Toolbox | Historical accuracy Box-Cox transformation in the behaviour workflow (`boxcox`) | Install through MATLAB; the [official `boxcox` documentation](https://www.mathworks.com/help/finance/boxcox.html) is in Financial Toolbox |
| SPM12 | Task-fMRI high-pass filtering (`spm_filter`) and haemodynamic response (`spm_hrf`) | Install [official SPM12](https://github.com/spm/spm12); configure `spm_root` / `SPM_ROOT` for the selected workflow |
| Full FieldTrip and HCP MEG Connectome tools | MEG preprocessing and source reconstruction | Configure the workflow's `FIELDTRIP_ROOT` and `MEGCONNECTOME_ROOT`; initialise full FieldTrip only in that workflow |

SPM is an external dependency because `spm_hrf` obtains defaults through
`spm_get_defaults` and `spm_defaults`; extracting a few functions without this
chain can silently change or break its behaviour. A separate installation also
permits use of the original analysis environment. Do not replace the historical
environment with a newer SPM or FieldTrip version when reproducing archived
results. Record the installed versions in each new run's provenance.

`check_dependencies` reports the resolved external paths. If SPM is available,
it also smoke-tests the filter and HRF. The bundled IO and optional graphics
checks use synthetic temporary inputs and delete those inputs afterwards.

The behaviour calculations retain MATLAB's original `boxcox` optimisation.
`check_dependencies` reports whether that entry point is available and its
resolved path. A missing function prevents transformed behaviour analyses;
there is no replacement with a different optimizer or a preselected lambda.
Untransformed input-loader checks do not verify this optimisation stage.
