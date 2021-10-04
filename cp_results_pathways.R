# Rscript cp_results_pathways.R lumB-vs-lumA-johansson-logTransf
# Rscript cp_results_pathways.R lumB-corr-johansson-logTransf
# Rscript cp_results_pathways.R lumB-corr-mertins-aggMean
# Rscript cp_results_pathways.R lumB-corr-krug-aggMean

require(foreach)

outFolder <- "CP_RESULTS_PATHWAYS"
dir.create(outFolder, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)

ds_folder <- "lumB-corr-mertins-aggMean"

ds_folder <- "lumB-vs-lumA-johansson-logTransf"

stopifnot(length(args) == 1)

ds_folder <- args[1]

sif_dt <- read.delim(file.path(ds_folder, "causative.sif"), header=F, col.names = c("p1", "link", "p2", "source", "na_col"), sep="\t")
# don't know why extra column with all na
stopifnot(is.na(sif_dt$na_col))
sif_dt$na_col <- NULL
sif_dt$full_link <- paste0(sif_dt$p1, "_", sif_dt$link, "_", sif_dt$p2) 
sif_dt$source <- gsub("\\s+:\\s+", ":", sif_dt$source) 
sif_dt$source <- gsub("\\s+:", ":", sif_dt$source) 
sif_dt$source <- gsub(":\\s+", ":", sif_dt$source) 


sif_dt$full_link <- paste0(sif_dt$p1, "_", sif_dt$link, "_", sif_dt$p2) 

out_dt <- foreach(i = 1:nrow(sif_dt), .combine='rbind') %do% {
  all_sources <- unlist(strsplit(sif_dt[i, "source"], split=" "))
  data.frame(
    links = rep(sif_dt[i, "full_link"], length(all_sources)), 
    source = all_sources,
    stringsAsFactors = FALSE
  )
}

pw_occ <- table(out_dt$source)
pw_occ <- sort(pw_occ, decreasing = TRUE)

agg_out_dt <- aggregate(links ~source, data=out_dt, function(x) paste0(x, collapse=","))
stopifnot(agg_out_dt$source %in% names(pw_occ))
agg_out_dt$occ <- pw_occ[agg_out_dt$source]
stopifnot(!is.na(agg_out_dt$occ))

agg_out_dt <- agg_out_dt[order(agg_out_dt$occ, decreasing = TRUE),]

outFile <- file.path(outFolder, paste0(ds_folder, "_pathways_occurences.txt"))
write.table(agg_out_dt, file = outFile, col.names=T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))