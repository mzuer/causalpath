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
  dt <- read.delim(file.path(x, "results.txt"))
  cat(paste0("... # relations=\t", nrow(dt), "\n"))
  cat(paste0("... # relations=\t", nrow(dt), "\n"), file=logFile, append = TRUE)
  if(cmpType == "corr") {
    dt <- dt[dt$Correlation.pval <= pval_filt,]
  }else if(cmpType == "cmp") {
    dt <- dt[dt$Source.change.pval <= pval_filt & dt$Target.change.pval <= pval_filt,]
  }
  dt$full_relation <- paste0(dt$Source, " ", dt$Relation, " ", dt$Target)
  cat(paste0("... # relations after pval filt.=\t", nrow(dt), "\n"))
  cat(paste0("... # relations after pval filt.=\t", nrow(dt), "\n"), file=logFile, append = TRUE)
  dt
  })

all_cmbs <- combn(x=1:length(all_results), m=2)

i=1
cmp_dt <- foreach(i = 1:ncol(all_cmbs),.combine='rbind') %do% {
  
  i_ds1 <- all_cmbs[1,i]
  i_ds2 <- all_cmbs[2,i]
  
  ds1 <- all_ds[i_ds1]
  ds2 <- all_ds[i_ds2]
  
  r1 <- all_results[[i_ds1]]
  r2 <- all_results[[i_ds2]]
  
  cat(paste0("... comparison: ", ds1, " vs. ", ds2, "\n"))
  cat(paste0("... comparison: ", ds1, " vs. ", ds2, "\n"), file=logFile, append = TRUE)
  
  stopifnot(!duplicated(r1$full_relation))
  stopifnot(!duplicated(r2$full_relation))
  
  unique_ds1 <- setdiff(r1$full_relation, r2$full_relation)
  cat(paste0("... only ", ds1, "\t= ",length(unique_ds1), "/", nrow(r1), " (",round(length(unique_ds1)/nrow(r1)*100, 2) , "%)\n"))
  cat(paste0("... only ", ds1, "\t= ",length(unique_ds1), "/", nrow(r1), " (",round(length(unique_ds1)/nrow(r1)*100, 2) , "%)\n"), file=logFile, append = TRUE)
  
  unique_ds2 <- setdiff(r2$full_relation, r1$full_relation)
  cat(paste0("... only ", ds2, "\t= ",length(unique_ds2), "/", nrow(r2), " (",round(length(unique_ds2)/nrow(r2)*100, 2) , "%)\n"))
  cat(paste0("... only ", ds2, "\t= ",length(unique_ds2), "/", nrow(r2), " (",round(length(unique_ds2)/nrow(r2)*100, 2) , "%)\n"), file=logFile, append = TRUE)
  
  inter_ds1_ds2 <- intersect(r2$full_relation, r1$full_relation)
  cat(paste0("... common ", ds1, " vs. ", ds2, "\t= ", length(inter_ds1_ds2), "\n"))
  cat(paste0("... common ", ds1, " vs. ", ds2, "\t= ", length(inter_ds1_ds2), "\n"), file=logFile, append = TRUE)
  
  u1 <- data.frame(
    ds_cmp = paste0(ds1, " vs. ", ds2),
    cmp_type = paste0("only_ds1"),
    relation = unique_ds1,
    stringsAsFactors = FALSE
  )
  u2 <- data.frame(
    ds_cmp = paste0(ds1, " vs. ", ds2),
    cmp_type = paste0("only_ds2"),
    relation = unique_ds2,
    stringsAsFactors = FALSE
  )
  if(length(inter_ds1_ds2) == 0) {
    inter12 <- data.frame(
      ds_cmp = paste0(ds1, " vs. ", ds2),
      cmp_type = paste0("common_ds1_ds2"),
      relation = NA,
      stringsAsFactors = FALSE
    )
  } else {
    inter12 <- data.frame(
      ds_cmp = paste0(ds1, " vs. ", ds2),
      cmp_type = paste0("common_ds1_ds2"),
      relation = inter_ds1_ds2,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, list(u1, u2, inter12))
}

inter_all <- Reduce(intersect, lapply(all_results, function(x) x$full_relation))

cat(paste0("... common to all\t= ",length(inter_all), "\n"))
cat(paste0("... common to all\t= ",length(inter_all), "\n"), file=logFile, append = TRUE)

if(length(inter_all) == 0) {
  cmp_dt2 <-data.frame(
    ds_cmp = paste0("all_ds"),
    cmp_type= paste0("common_all"),
    relation= NA,
    stringsAsFactors = FALSE)
}else {
  cmp_dt2 <-data.frame(
    ds_cmp = paste0("all_ds"),
    cmp_type= paste0("common_all"),
    relation= inter_all,
    stringsAsFactors = FALSE)
  
}


out_dt <- rbind(cmp_dt, cmp_dt2)

occurence_dt <- data.frame(table(unlist(lapply(all_results, function(x) x$full_relation))))
colnames(occurence_dt) <- c("relation","frequency")
occurence_dt <- occurence_dt[order(occurence_dt$relation, decreasing = FALSE),]
occurence_dt <- occurence_dt[order(occurence_dt$frequency, decreasing = TRUE),]

outFile <- file.path(outFolder, "relation_frequencies.txt")
write.table(occurence_dt, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=F)
cat(paste0(".... written: ", outFile, "\n"))


outFile <- file.path(outFolder, "relation_intersects.txt")
write.table(out_dt, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=F)
cat(paste0(".... written: ", outFile, "\n"))
cat(paste0(".... written: ", logFile, "\n"))




tmp_dt <- out_dt
tmp_dt$breaklab <-paste0(tmp_dt$ds_cmp, tmp_dt$cmp_type)
i = unique(tmp_dt$breaklab)[1]
for(i in unique(tmp_dt$breaklab)) {
  tmp2_dt <- tmp_dt[tmp_dt$breaklab == i,]
  stopifnot(nrow(tmp2_dt) > 0)
  
  sif_dt <- data.frame(
    col1 = unlist(lapply(tmp2_dt$relation, function(x) unlist(strsplit(x, split=" "))[1])),
    col2 = unlist(lapply(tmp2_dt$relation, function(x) unlist(strsplit(x, split=" "))[2])),
    col3 = unlist(lapply(tmp2_dt$relation, function(x) unlist(strsplit(x, split=" "))[3])), 
    stringsAsFactors=FALSE
  ) 
  
  outFile <- file.path(outFolder, paste0(i, "_relations.sif"))
  write.table(sif_dt, file=outFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=F)
  cat(paste0(".... written: ", outFile, "\n"))
  cat(paste0(".... written: ", logFile, "\n"))
  
}



stop("**** DONE - ok\n")
