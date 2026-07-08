## Same benchmark as dev/bench_access.R but against the working tree via load_all.
suppressWarnings(suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = FALSE)))
options(neuroim2.bench.loaded = TRUE)
source("dev/bench_access.R", local = FALSE, chdir = FALSE)
