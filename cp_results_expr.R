ds_folder <- "lumB-vs-lumA-johansson-logTransf"

p1 <- "EZH2"
p2 <- "PCNA"

param_data <- readLines(file.path(ds_folder, "parameters.txt"))

result_dt <- read.delim(file.path(ds_folder, "results.txt"))
result_dt$full_relation <- paste0(dt$Source, " ", dt$Relation, " ", dt$Target)


# Rscript cmp_CP_results.R lumB-vs-lumA-johansson lumB-vs-lumA-krug-aggMean cmp
# Rscript cmp_CP_results.R lumB-vs-lumA-johansson-logTransf lumB-vs-lumA-krug-aggMean lumB-vs-lumA-mertins-aggMean cmp
# Rscript cmp_CP_results.R lumB-corr-johansson-logTransf lumB-corr-krug-aggMean lumB-corr-mertins-aggMean corr

require(foreach)

pval_filt <- 0.05

# all_ds <- c("lumB-vs-lumA-johansson", "lumB-vs-lumA-krug-aggMean", "lumB-vs-lumA-mertins-aggMean")
# all_ds <- c("lumB-corr-johansson-logTransf", "lumB-corr-krug-aggMean")

all_ds_withcmp <- commandArgs(trailingOnly = TRUE)

stopifnot(length(all_ds_withcmp) >= 2)

all_ds <- all_ds_withcmp[1:(length(all_ds_withcmp) -1)]
cmpType <- all_ds_withcmp[length(all_ds_withcmp)]
stopifnot(cmpType %in% c("corr", "cmp"))

outFolder <- file.path("CMP_CP_RESULTS", paste0(paste0(all_ds,  collapse = "_vs_"), "_pvalfilt", pval_filt))
dir.create(outFolder, recursive = TRUE)

logFile <- file.path(outFolder, "logFile.txt")

cat(paste0("... relation p-val filt = ",pval_filt, "\n"))
cat(paste0("... relation p-val filt = ",pval_filt, "\n"), file=logFile, append = TRUE)


cat(paste0("# all DS: ", paste0(all_ds, collapse= " vs. "), "\n"))
cat(paste0("# all DS: ", paste0(all_ds, collapse= " vs. "), "\n"), file=logFile, append = TRUE)

# if corr -> Correlation.pval
# Target.change.pval
all_results <- lapply(all_ds, function(x) {
  cat(paste0("> ", x, "\n"))
  cat(paste0("> ", x, "\n"), file=logFile, append = TRUE)
  cat(paste0("... # relations=\t", nrow(dt), "\n"))
  cat(paste0("... # relations=\t", nrow(dt), "\n"), file=logFile, append = TRUE)
  if(cmpType == "corr") {
    dt <- dt[dt$Correlation.pval <= pval_filt,]
  }else if(cmpType == "cmp") {
    dt <- dt[dt$Source.change.pval <= pval_filt & dt$Target.change.pval <= pval_filt,]
  }
  cat(paste0("... # relations after pval filt.=\t", nrow(dt), "\n"))
  cat(paste0("... # relations after pval filt.=\t", nrow(dt), "\n"), file=logFile, append = TRUE)
  dt
})