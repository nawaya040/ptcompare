# Submission comments

This is a new submission.

## Test environments

* local Windows 11, R 4.5.2, Rtools45 — `R CMD check --as-cran`
* win-builder, R-devel — pending
* win-builder, R-release — pending
* macOS builder (macbuilder.r-project.org) — pending

## R CMD check results

Local `R CMD check --as-cran`: no ERRORs, no WARNINGs, 2 NOTEs.

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Naoki Awaya <nawaya@waseda.jp>'

New submission

Found the following (possibly) invalid URLs:
  URL: https://github.com/nawaya040/ptcompare
  URL: https://github.com/nawaya040/ptcompare/issues
    Status: 404
```

The first part is expected for a first submission. The 404s are because the
repository is still private; it is made public before submission, after which
both URLs resolve. (If it stays private, the `URL` and `BugReports` fields are
removed from DESCRIPTION instead.)

```
* checking for future file timestamps ... NOTE
unable to verify current time
```

This is the check being unable to reach the time service from the machine it
ran on, not a property of the package.

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
