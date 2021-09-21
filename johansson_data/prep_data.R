# Rscript prep_data.R

startTime <- Sys.time()

cat(paste0("... start prep_data.R \n"))

outFolder <- file.path("..", "johansson_lumA_lumB_CP_inputs")
dir.create(outFolder, recursive = TRUE)
outFile <- file.path(outFolder, "johansson_lumA_lumB_cp_input.txt")

# manually defined lumA and lumB samples
lumB <- c("OSL.56F", "OSL.405", "OSL.49E","OSL.524", "OSL.443", "OSL.458", "OSL.521", "OSL.46D", "OSL.54D")
lumA <- c("OSL.567", "OSL.4B0", "OSL.485", "OSL.41B", "OSL.4AF", "OSL.46E", "OSL.494", "OSL.457","OSL.48B")
her2 <- c("OSL.43C", "OSL.493","OSL.4D9", "OSL.43A", "OSL.406", "OSL.47C", "OSL.40A", "OSL.40E", "OSL.3FA")
basal <- c("OSL.53E", "OSL.3FF", "OSL.55F", "OSL.46A", "OSL.4D6", "OSL.4B4", "OSL.449", "OSL.44E", "OSL.3EB")
normal <- c("OSL.441", "OSL.430", "OSL.4FA", "OSL.53D","OSL.540", "OSL.42E", "OSL.4BA", "OSL.579", "OSL.57B")

in_dt <- read.delim(file.path("johansson_data_relative_ratios_to_pool.csv"), sep=",",header=TRUE)
nSamp <- 45

stopifnot(length(lumA) == 9)
stopifnot(length(lumB) == 9)
stopifnot(length(her2) == 9)
stopifnot(length(normal) == 9)
stopifnot(length(basal) == 9)

stopifnot(lumA %in% colnames(in_dt))
stopifnot(lumB %in% colnames(in_dt))
stopifnot(her2 %in% colnames(in_dt))
stopifnot(normal %in% colnames(in_dt))
stopifnot(basal %in% colnames(in_dt))

stopifnot(ncol(in_dt) == nSamp+2)


## input for causalpath should have following columns
# 1) ID: A unique text identifier
# 2) HGNC symbol
# 3) Sites: If phosphoprotein measurement, protein sites that are affected  (e.g. S79|S80 S222 -> multiple sites and multiple symbs)
# 4) Effect: If phosphoprotein measurement, effect of phosphorylation on the protein activity (a, i,blank or c)
# 5) Value: numeric value; can be > 1 columns

samp_ids <- colnames(in_dt)[!colnames(in_dt) %in% c("gene_symbol", "ensembl_id")]
stopifnot(length(samp_ids) == nSamp)
stopifnot(!duplicated(samp_ids))

firstColLabs <- c("ID", "Symbols", "Sites", "Effect")

stopifnot(!duplicated(in_dt$ensembl_id))
# -> no one duplicated -> will use it as ID
in_dt$ID <- in_dt$ensembl_id
in_dt$Symbols <- in_dt$gene_symbol
in_dt$Sites <- ""
in_dt$Effect <- ""

input_dt <- in_dt[,c(firstColLabs, lumA, lumB)]

stopifnot(ncol(input_dt) == length(lumA) + length(lumB) + length(firstColLabs))

stopifnot(nrow(input_dt) == nrow(in_dt))


write.table(input_dt, file=outFile, sep="\t", col.names=TRUE, row.names=FALSE)
cat(paste0("... written: ", outFile, "\n"))

cat(paste0("**** Done... \n"))
cat(paste0(startTime, "\n", Sys.time(), "\n"))