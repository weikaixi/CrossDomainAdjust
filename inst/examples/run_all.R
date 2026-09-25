library(CrossDomainAdjust)
folder <- system.file("examples", package = "CrossDomainAdjust")
files <- list.files(folder, pattern = "^[0-9][0-9]_.*[.]R$", full.names = TRUE)
stopifnot(length(files) == 15L)
for (file in files) {
  cat("
Running", basename(file), "
")
  source(file, local = new.env(parent = globalenv()), echo = FALSE)
}
cat("All 15 worked examples completed.
")
