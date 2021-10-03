# Rscript cmp_corr_results.R

outFolder <- file.path("CMP_CORR_RESULTS")
dir.create(outFolder, recursive = TRUE)

### results lumB
lumB_dt <- read.delim("lumB-corr-cptac-babur/results.txt")
nrow(lumB_dt)
# 22096
sum(lumB_dt$Correlation.pval <= 0.05)
# 3753
lumB_dt$link <- paste0(lumB_dt$Source, "_", lumB_dt$Relation, "_", lumB_dt$Target)
length(unique(lumB_dt$link))
# 7992
length(unique(lumB_dt$link[lumB_dt$Correlation.pval <= 0.05]))
# 1929

### results lumA
lumA_dt <- read.delim("lumA-corr-cptac-babur/results.txt")
nrow(lumA_dt)
# 22298
sum(lumA_dt$Correlation.pval <= 0.05)
# 4283
lumA_dt$link <- paste0(lumA_dt$Source, "_", lumA_dt$Relation, "_", lumA_dt$Target)
length(unique(lumA_dt$link))
# 7915
length(unique(lumA_dt$link[lumA_dt$Correlation.pval <= 0.05]))
# 2009

### results Her2
her2_dt <- read.delim("her2-corr-cptac-babur/results.txt")
nrow(her2_dt)
# 17523
sum(her2_dt$Correlation.pval <= 0.05)
# 2090
her2_dt$link <- paste0(her2_dt$Source, "_", her2_dt$Relation, "_", her2_dt$Target)
length(unique(her2_dt$link))
# 7247
length(unique(her2_dt$link[her2_dt$Correlation.pval <= 0.05]))
# 1362


length(intersect(unique(her2_dt$link), unique(lumB_dt$link)))
# 5547

length(intersect(unique(her2_dt$link), unique(lumA_dt$link)))
# 5630

length(intersect(unique(lumB_dt$link), unique(lumA_dt$link)))
# 6011

interLinks_lumB_her2 <- intersect(unique(her2_dt$link), unique(lumB_dt$link))
diffLinks_lumBnotLumA <- setdiff(unique(lumB_dt$link), unique(lumA_dt$link))

interLinks_lumB_her2_lumBsignif <- intersect(unique(her2_dt$link), unique(lumB_dt$link[lumB_dt$Correlation.pval <= 0.05]))
diffLinks_lumBnotLumA_lumBsignif <- setdiff(unique(lumB_dt$link[lumB_dt$Correlation.pval <= 0.05]), 
                                            unique(lumA_dt$link))


diffLinks_her2notLumA <- setdiff(unique(her2_dt$link), unique(lumA_dt$link))

lumB_her2_notLumA <- intersect(diffLinks_her2notLumA,diffLinks_lumBnotLumA )
lumB_her2_notLumA_lumBsignif <- intersect(diffLinks_her2notLumA,diffLinks_lumBnotLumA_lumBsignif )


######################################################################
# take the sif file, and retrieve only the links that are common

lumB_sif_dt <- read.delim("lumB-corr-cptac-babur/causative.sif", header=F)
colnames(lumB_sif_dt)[1:3] <- c("Source", "Relation", "Target") 
lumB_sif_dt$link <- paste0(lumB_sif_dt$Source, "_", lumB_sif_dt$Relation, "_", lumB_sif_dt$Target)
stopifnot(lumB_sif_dt$link %in% lumB_dt$link)

interHe2_lumB_sif_dt <- lumB_sif_dt[lumB_sif_dt$link %in% interLinks_lumB_her2,]
stopifnot(length(interLinks_lumB_her2) == nrow(interHe2_lumB_sif_dt))

stopifnot(diffLinks_lumBnotLumA %in% lumB_sif_dt$link)
notLumA_sif_dt <- lumB_sif_dt[lumB_sif_dt$link %in% diffLinks_lumBnotLumA,]
stopifnot(length(diffLinks_lumBnotLumA) == nrow(notLumA_sif_dt))

interHe2_lumB_sif_dt$link <- notLumA_sif_dt$link<- NULL

outFile <- file.path(outFolder, "lumB_corr_cptac_babur_interHe2_lumB_causative.sif")
write.table(interHe2_lumB_sif_dt, file=outFile, col.names = F, row.names=F, sep="\t", quote=FALSE)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "lumB_corr_cptac_babur_notLumA_lumB_causative.sif")
write.table(notLumA_sif_dt, file=outFile, col.names = F, row.names=F, sep="\t", quote=FALSE)
cat(paste0("... written: ", outFile, "\n"))


signifOnly_interHe2_notLumA_lumB_sif_dt <- lumB_sif_dt[lumB_sif_dt$link %in% lumB_her2_notLumA_lumBsignif,]
stopifnot(length(lumB_her2_notLumA_lumBsignif) == nrow(signifOnly_interHe2_notLumA_lumB_sif_dt))

outFile <- file.path(outFolder, "lumB_corr_cptac_babur_interHe2notLumA_lumB_signifOnly_causative.sif")
write.table(signifOnly_interHe2_notLumA_lumB_sif_dt, file=outFile, col.names = F, row.names=F, sep="\t", quote=FALSE)
cat(paste0("... written: ", outFile, "\n"))

interHe2_notLumA_lumB_sif_dt <- lumB_sif_dt[lumB_sif_dt$link %in% lumB_her2_notLumA,]
stopifnot(length(lumB_her2_notLumA) == nrow(interHe2_notLumA_lumB_sif_dt))

outFile <- file.path(outFolder, "lumB_corr_cptac_babur_interHe2notLumA_lumB_causative.sif")
write.table(interHe2_notLumA_lumB_sif_dt, file=outFile, col.names = F, row.names=F, sep="\t", quote=FALSE)


signifOnly_interHe2_lumB_sif_dt <- lumB_sif_dt[lumB_sif_dt$link %in% interLinks_lumB_her2_lumBsignif,]
stopifnot(length(interLinks_lumB_her2_lumBsignif) == nrow(signifOnly_interHe2_lumB_sif_dt))

signifOnly_notLumA_sif_dt <- lumB_sif_dt[lumB_sif_dt$link %in% diffLinks_lumBnotLumA_lumBsignif,]
stopifnot(length(diffLinks_lumBnotLumA_lumBsignif) == nrow(signifOnly_notLumA_sif_dt))

signifOnly_notLumA_sif_dt$link <- signifOnly_interHe2_lumB_sif_dt$link<- NULL

outFile <- file.path(outFolder, "lumB_corr_cptac_babur_interHe2_lumB_signifOnly_causative.sif")
write.table(signifOnly_interHe2_lumB_sif_dt, file=outFile, col.names = F, row.names=F, sep="\t", quote=FALSE)
cat(paste0("... written: ", outFile, "\n"))

outFile <- file.path(outFolder, "lumB_corr_cptac_babur_notLumA_lumB_signifOnly_causative.sif")
write.table(signifOnly_notLumA_sif_dt, file=outFile, col.names = F, row.names=F, sep="\t", quote=FALSE)
cat(paste0("... written: ", outFile, "\n"))


selected_symbs <- c( "DTX3", "MRPS23", "EIF2S2", "EIF6", "SLC2A10")
any(selected_symbs %in% lumB_dt$Source)
lumB_dt[lumB_dt$Source %in% selected_symbs |lumB_dt$Target %in% selected_symbs ,]

interHe2_notLumA_lumB_sif_dt[interHe2_notLumA_lumB_sif_dt$Source %in% selected_symbs |
                               interHe2_notLumA_lumB_sif_dt$Target %in% selected_symbs ,]



prolif_dt <- read.delim("CELL_PROLIFERATION_GO_0008283.txt", skip=2, header=F)

interHe2_notLumA_lumB_sif_dt[interHe2_notLumA_lumB_sif_dt$Source %in% prolif_dt[,1] |
                               interHe2_notLumA_lumB_sif_dt$Target %in% prolif_dt[,1] ,]


interHe2_notLumA_lumB_sif_dt[interHe2_notLumA_lumB_sif_dt$Source %in% prolif_dt[,1] &
                               interHe2_notLumA_lumB_sif_dt$Target %in% prolif_dt[,1] ,]

signifOnly_interHe2_notLumA_lumB_sif_dt[signifOnly_interHe2_notLumA_lumB_sif_dt$Source %in% prolif_dt[,1] |
                                          signifOnly_interHe2_notLumA_lumB_sif_dt$Target %in% prolif_dt[,1] ,]

signifOnly_interHe2_notLumA_lumB_sif_dt[signifOnly_interHe2_notLumA_lumB_sif_dt$Source %in% prolif_dt[,1] &
                                          signifOnly_interHe2_notLumA_lumB_sif_dt$Target %in% prolif_dt[,1] ,]

