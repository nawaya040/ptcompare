# Installs this source package into a project-local library.
# Requires Rcpp and, for compilation, Rtools on Windows / a C++ compiler elsewhere.
stopifnot(file.exists("DESCRIPTION"))
if ("ptcompare" %in% loadedNamespaces())
  stop("Restart R before installing (Session > Restart R in RStudio).")
dir.create(".local-library", showWarnings = FALSE)
.libPaths(c(normalizePath(".local-library", winslash = "/"), .libPaths()))
if (!requireNamespace("Rcpp", quietly = TRUE))
  install.packages("Rcpp", repos = "https://cloud.r-project.org")
if (.Platform$OS.type == "windows") {
  rtools <- paste0("C:/rtools", paste0(strsplit(as.character(getRversion()), ".", fixed=TRUE)[[1]][1:2], collapse=""))
  if (dir.exists(rtools)) {
    tools_path <- c(file.path(rtools, "usr/bin"),
                    file.path(rtools, "x86_64-w64-mingw32.static.posix/bin"))
    Sys.setenv(PATH = paste(c(tools_path, Sys.getenv("PATH")), collapse = .Platform$path.sep))
  }
}
rbin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "R.exe" else "R")
status <- system2(rbin, c("CMD", "INSTALL", "--preclean", "--library=.local-library", "."))
if (status != 0) stop("Installation failed; see compiler messages above.")
cat("Installed. Run source('examples/01-local-shift.R').\n")
