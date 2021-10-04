# Rscript cp_results_expr.R lumB-vs-lumA-johansson-logTransf EZH2 PCNA cmp
# Rscript cp_results_expr.R lumB-corr-johansson-logTransf ETS1 FLI1 corr

require(ggplot2)

args <- commandArgs(trailingOnly = TRUE)


ds_folder <- "lumB-vs-lumA-johansson-logTransf"
p1 <- "EZH2"
p2 <- "PCNA"
cmpType <- "cmp"

ds_folder <- "lumB-corr-johansson-logTransf"
p1 <- "ETS1"
p2 <- "FLI1"
cmpType <- "corr"

stopifnot(length(args) == 4)

ds_folder <- args[1]
p1 <- args[2]
p2 <- args[3]
cmpType <- args[4]


plotType <- "png"
myWidthGG <- 6
myHeightGG <- 6


stopifnot(cmpType %in% c("cmp", "corr"))


outFolder <- file.path("CP_RESULTS_EXPR")
dir.create(outFolder, recursive = TRUE)

param_data <- readLines(file.path(ds_folder, "parameters.txt"))

result_dt <- read.delim(file.path(ds_folder, "results.txt"))

result_dt$full_relation <- paste0(result_dt$Source, " ", result_dt$Relation, " ", result_dt$Target)


logTransf_all <- param_data[grepl("do-log-transform", param_data)]

logTransf_param <- Filter(function(x) {strsplit(x, split="")[1] != "#"}, trimws(logTransf_all))
if(length(logTransf_param) == 0) {
  logTransf_bool <- FALSE
} else {
  stopifnot(length(logTransf_param) == 1)
  logTransf_bool <- unlist(lapply(logTransf_param, function(x) trimws(unlist(strsplit(x, split="="))[[2]])))
  stopifnot(logTransf_bool == "true" | logTransf_bool == "false")
  logTransf_bool <- toupper(logTransf_bool)
}


expr_file_all <- param_data[grepl("proteomics-values-file", param_data)]
expr_file <- Filter(function(x) {strsplit(x, split="")[1] != "#"}, trimws(expr_file_all))
stopifnot(length(expr_file) == 1)
exprFile <- unlist(lapply(expr_file, function(x) trimws(unlist(strsplit(x, split="="))[[2]])))
stopifnot(file.exists(exprFile))
expr_dt <- read.delim(exprFile, header = TRUE)
### !!! NB: next line will work only if pass in the same file (can be passed in proteomics-platform-file)
stopifnot(c("ID", "Symbols", "Sites", "Effect") %in% colnames(expr_dt))
stopifnot(p1 %in% expr_dt$Symbols)
stopifnot(p2 %in% expr_dt$Symbols)

value_cols_all <- param_data[grepl("value-column", param_data)]
# remove lines with comment:
value_cols <- Filter(function(x) {strsplit(x, split="")[1] != "#"}, trimws(value_cols_all))

samptypes <- unlist(lapply(value_cols, function(x) trimws(unlist(strsplit(x, split="="))[[1]])))

sampnames <- unlist(lapply(value_cols, function(x) trimws(unlist(strsplit(x, split="="))[[2]])))
stopifnot(sampnames %in% colnames(expr_dt))

sampLabs <- setNames(gsub("-value-column", "", samptypes), sampnames)

plot_dt <- expr_dt[expr_dt$Symbols == p1 | expr_dt$Symbols == p2, colnames(expr_dt) %in% c("Symbols", sampnames)]
stopifnot(nrow(plot_dt) == 2)
## with some datasets might not always be true -> might need to be adapted to ID instead
stopifnot(!duplicated(plot_dt$Symbols))
rownames(plot_dt) <- plot_dt$Symbols
plot_dt$Symbols <- NULL
if(logTransf_bool){
  plot_dt <- log2(plot_dt+1)  
}
plot_dt <- as.data.frame(t(plot_dt))
stopifnot(p1 %in% colnames(plot_dt) & p2 %in% colnames(plot_dt))
stopifnot(sampnames %in% rownames(plot_dt))
plot_dt$sampType <- sampLabs[rownames(plot_dt)]
stopifnot(!is.na(plot_dt$sampType))

mytheme <-   theme(plot.title = element_text(hjust=0.5, face = "bold"),
                   plot.subtitle = element_text(hjust=0.5),
                   strip.background = element_rect(fill = "grey95", color = "grey30", size = 1),
                   legend.text = element_text(size = 10),
                   strip.text = element_text(size=12),
                   axis.text = element_text(size = 10),
                   axis.title = element_text(size = 12),
                   legend.title = element_text(size = 10))


mytit <- paste0(ds_folder)
mysub_ <- basename(expr_file)
mysubFile <- basename(expr_file)
if(logTransf_bool) {
  mysub <- paste0(mysub_, " [log2(.+1)]")
  mysubFile <- paste0(mysub_, "_log2")
}


if(cmpType=="corr"){
  stopifnot(samptypes == "value-column")
  p <- ggplot(plot_dt, aes_string(x=paste0(p1), y=paste0(p2), color="sampType")) +
    stat_smooth(se = FALSE, method = "lm", size = 1,
                color = 'grey30') +
    geom_point(size=3) +
    scale_color_brewer(palette="Set1")+
    guides(color=FALSE)
  
} else {
  stopifnot(samptypes %in% c("test-value-column", "control-value-column"))
  plot_dt_nocond <- plot_dt
  plot_dt_nocond$sampType <- NULL
  p <- ggplot(plot_dt, aes_string(x=paste0(p1), y=paste0(p2))) +
    
    geom_point(data = plot_dt_nocond, fill = "grey80", color = "grey80", size = 3) +
    geom_point(aes(color = sampType),  size = 3) +
    labs(color="")+
    facet_wrap(~sampType)+
    stat_smooth(se = FALSE, method = "lm", size = 1,
                color = 'grey30') +
    scale_color_brewer(palette="Set1") 
  myWidthGG <- myWidthGG * 1.8
}


p <- p + 
  ggtitle(paste0(mytit), subtitle=paste0(mysub))+
  ylab(paste0(p2))+
  xlab(paste0(p1))+
  mytheme 


outFile <- file.path(outFolder, paste0(p1, "_vs_", p2, "_", mysubFile, ".", plotType))
ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))



if(cmpType == "cmp") {
  
  stopifnot(samptypes %in% c("test-value-column", "control-value-column"))

    p <- ggplot(plot_dt, aes_string(x=paste0(p1), y=paste0(p2), color="sampType")) +
    geom_point(size = 3) +
    labs(color="")+
    facet_wrap(~sampType, scales = "free")+
    stat_smooth(se = FALSE, method = "lm", size = 1,
                color = 'grey30') +
    scale_color_brewer(palette="Set1") + 
    ggtitle(paste0(mytit), subtitle=paste0(mysub))+
    ylab(paste0(p2))+
    xlab(paste0(p1))+
    mytheme 

    outFile <- file.path(outFolder, paste0(p1, "_vs_", p2, "_", mysubFile, "_v2.", plotType))
    ggsave(p, filename = outFile, height=myHeightGG, width=myWidthGG)
    cat(paste0("... written: ", outFile, "\n"))
    
    
}



stop("-ok\n")

require(paxtoolsr)

getURI()


