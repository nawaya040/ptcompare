# Submission comments

This is a new submission. Submitted version: 0.1.0.

## Test environments

* local Windows 11, R-devel (2026-10-26 r90766 ucrt), Rtools45 —
  `R CMD check --as-cran` on 0.1.0
* local Windows 11, R 4.5.2, Rtools45 — `R CMD check --as-cran` on 0.1.0
* win-builder, R-release (R 4.6.1, Windows Server 2022, x86_64) — 0.1.0
* win-builder, R-devel — 0.1.0 and 0.0.1 were uploaded; no result was
  returned for either; see below
* macOS builder (macbuilder.r-project.org) — the run did not check the
  package; see below

0.0.1 and 0.1.0 differ only in the version string, `NEWS.md`, the README, the
heading of the statistical specification and the vignette prose. No R or C++
source file differs between them.

## R CMD check results

### local, R-devel (2026-10-26 r90766 ucrt), `--as-cran` (0.1.0)

No ERRORs, no WARNINGs, 1 NOTE. The checked tarball is the one uploaded to
CRAN. Rcpp and the Suggests packages were installed from CRAN for this R-devel
build. Examples, tests, the vignette and both versions of the manual were
built and run.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Naoki Awaya <nawaya@waseda.jp>'

New submission
```

### local, R 4.5.2, `--as-cran` (0.1.0)

No ERRORs, no WARNINGs, 2 NOTEs. Installation 13s. Examples, tests, the
vignette and both versions of the manual were built and run.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Naoki Awaya <nawaya@waseda.jp>'

New submission

* checking for future file timestamps ... NOTE
unable to verify current time
```

The second is the check being unable to reach the time service from the
machine it ran on. It does not appear on win-builder.

### win-builder, R-release (0.1.0)

Status: 1 NOTE. No ERRORs, no WARNINGs. Installation 22s, check 126s, and the
Windows binary was built. Examples, tests, the vignette and both versions of
the manual were built and run. The run on 0.0.1 gave the same NOTE. The
tarball checked here was built before one section heading was removed from the
vignette; the uploaded tarball differs from it only in that heading.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Naoki Awaya <nawaya@waseda.jp>'

New submission

Possibly misspelled words in DESCRIPTION:
  Awaya (8:67)
  Multiscale (3:17)
  Polya (8:19)
```

All three are intentional. "Awaya" is an author surname, "Multiscale" is used
in the title in its usual statistical sense, and "Polya" is the ASCII
spelling of Pólya, used because the DESCRIPTION file is ASCII. A spell check
of the manual, the vignette and `NEWS.md` reports no other words outside the
dictionary except proper nouns (Mangasarian, Soriano, Wolberg), package names
and the terms "midranks" and "subregions".

### win-builder, R-devel

0.0.1 was uploaded on 2026-10-27 and 0.1.0 on 2026-10-28. Both left the
upload queue, but neither run sent a result. On 2026-10-28 the R-devel queue
held other packages for more than 30 minutes while the R-release queue was
empty. The R-devel check was therefore run locally; see above.

### macOS builder

The submitted run stopped at the dependency step after 22 seconds:

```
* checking package dependencies ... ERROR
Package required but not available: 'Rcpp'
Packages suggested but not available for checking:
  'testthat', 'knitr', 'rmarkdown', 'mclust'
```

None of the dependencies were present in the builder's library, so the
package was never compiled, and no example, test or vignette was run. This is
a property of that run rather than of the package.

## Notes for the reviewer

* The package contains compiled C++ (C++17, via Rcpp). It uses no external
  libraries, writes no files and makes no network connections.
* The vignette uses `mclust` for one example data set. `mclust` is in Suggests
  and every chunk that needs it is guarded with
  `eval = requireNamespace("mclust", quietly = TRUE)`, so the vignette builds
  without it.
* `inst/CITATION` points to the paper the method is taken from: Awaya, N. and
  Ma, L. (2024), Journal of the American Statistical Association, 119(545),
  189--201, <doi:10.1080/01621459.2022.2105223>.
* The README is written in Japanese; the manual and the vignette are in
  English.
