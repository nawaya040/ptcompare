# Open ptcompare.Rproj before sourcing an example.
if (!file.exists("DESCRIPTION") || !dir.exists(".local-library")) {
  stop("Open ptcompare.Rproj (project root), then run scripts/install-local.R once.")
}
.pt_local_library <- normalizePath(".local-library", winslash = "/", mustWork = TRUE)
.libPaths(c(.pt_local_library, .libPaths()))
if (!requireNamespace("ptcompare", quietly = TRUE)) {
  stop("Run source('scripts/install-local.R') first.")
}
library(ptcompare)
rm(.pt_local_library)
