# Third-party software

These files retain their original copyright notices and licences. The repository's
MIT licence applies to the authors' code; it does not relicense this directory.

| Directory | Included use | Licence and version |
| --- | --- | --- |
| `nanconv/` | NaN-aware convolution used by the velocity-field calculations | Benjamin Kraus, revision 4861 (2013-05-27); BSD-style terms in `license.txt` |
| `cifti_matlab/` | Legacy `ft_read_cifti`, `ft_write_cifti`, GIfTI and XMLTree; source fallbacks included | FieldTrip snapshot identified by `VERSION` as `d6ada706fdb09f824c38e4de7655757076be4a21` (2016-04-16), with archived changes and release additions recorded below; GPL-3.0-or-later for FieldTrip, GPL-2.0-or-later for GIfTI/XMLTree |
| `raincloudplots/` | Optional `raincloud_plot` and `rm_raincloud` graphics helpers | Archived original-project snapshot; MIT, original `LICENSE` retained |

`MANIFEST.json` identifies each shipped dependency file by SHA-256. A snapshot
without recoverable upstream revision metadata is identified by these content
hashes, not assigned a guessed tag or commit.

Only the CIFTI dependency subset is included, not full FieldTrip. Platform MEX
binaries and tutorial datasets are omitted. The bundled XML search uses its
MATLAB source fallback; compressed GIfTI uses the original JVM fallback.

The CIFTI copy contains these explicit dependency repairs:

- `private/issubfield.m`, `private/getsubfield.m`, `private/setsubfield.m` are
  unmodified files from the FieldTrip commit listed above. They close the
  dependencies of the existing warning handler.
- `@gifti/private/xml_findstr.m` is copied from the archived
  `@xmltree/private/xml_findstr.m`, preserving its attribution. This permits GIfTI
  XML parsing without a platform-specific MEX binary.
- `cpm_check_private_helpers.m` is an author-written MIT-licensed integration test.

The original archive's `GPL` text is GPL version 2. `GPL-3.0.txt` additionally
provides the GPL version 3 terms referenced by individual FieldTrip source
headers. Original notices are retained even where the archive-level notice is
older than an individual source header. The included `miniz.c` is version 1.15
and retains Rich Geldreich's public-domain/Unlicense notice. C sources are
included for completeness; the tested MATLAB path uses source/JVM fallbacks.

The original archive mixed a modern `cifti_read` interface with legacy private
functions whose signatures differ. That incompatible mixture is not shipped:
the manuscript's `ft_read_cifti` interface and its matching private functions are
kept together. No scientific calculation in the reader or writer is changed.

See [dependency setup and checks](../docs/DEPENDENCIES.md).
