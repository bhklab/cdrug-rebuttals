########################
## Benjamin Haibe-Kains
## All rights Reserved
## April 8, 2014
########################



## This script runs the following scripts:
##  1. 'normalization_cgp.R' performing curation, annotation and normalization of CGP data
##  2. 'normalization_ccle.R' performing curation, annotation and normalization of CCLE data
##  3. 'cdrug2_format.R' performing additional curation to identify common cell lines, tissue types and drugs investigated in CGP and CCLE
##  4. 'cdrug2_analysis.R' performing the analyses presented in our response to Gelleher et al.
##  5. 'cdrug2_pipeline.R' computing all the correlations and generating all the tables/figures for the paper


## remove all existing objects from the workspace
rm(list=ls(all=TRUE))

require(parallel) || stop("Library parallel is not available")

## define functions required for the analysis pipeline
require(MetaGx) || stop("Library MetaGx is not available")
require(PharmacoGx) || stop("Library PharmacoGx is not available")
require(genefu) || stop("Library genefu is not available")
require(gdata) || stop("Library gdata is not available")

## additional functions
source(file.path("code", "cdrug2_foo.R"))

########################
## global parameters

## set method for downloading
# options(download.file.method="auto")
options(download.file.method="wget")
## change to curl, wget or internal depending on your system

## prevent strings to be converted into factors
options(stringsAsFactors=FALSE)

## set random seed to ensuer reproducibility of the resuls
set.seed(54321)

## number of cpu cores available for the analysis pipeline
## set to 'NULL' if all the available cores should be used
nbcore <- 16
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)

## tissue type specific analyses
tissue.specific <- FALSE

## list of characters to be removed from row and column names
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

## directory where all the analysis results will be stored
saveres <- "saveres"
if(!file.exists(saveres)) { dir.create(saveres, showWarnings=FALSE, recursive=TRUE) }

## maximum number of CEL files to normalize at once
max.celfiles <- 300

## minimum number of samples to compute correlation
minsample <- 10

## minimum number of samples in each category for the adaptive MCC computation
min.cat <- 3

## max fdr, threshold used to identify genes with significant concordance index
myfdr <- 0.20

## method to estimate the association gene-drug, controlled for tissue type
# genedrugm <- c("lm", "cindex")
genedrugm <- "lm"

## consistency between associations
concordance.method <- "spearman"
fdr.cuts <- c(1, 0.5, 0.2, 0.05, 0.01)
quantile.cuts <- c(1, 0.5, 0.2, 0.05, 0.01)

## list of genes to consider for concordance
## listg <- c("all", "l1000")
listg <- "all"

## Broad landmark genes
dir.create(file.path("data", "L1000", "dwl"), recursive=TRUE, showWarnings=FALSE)
myfn <- file.path("saveres", "l1000_genes.RData")
if (!file.exists(myfn)) {
  message("Download landmark genes from L1000")
  dwl.status <- download.file(url="http://www.lincscloud.org/l1000/example_files/Landmark_Genes_n978.xlsx", destfile=file.path("data", "L1000", "dwl", "Landmark_Genes_n978.xlsx"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path("data", "L1000", "dwl", "Landmark_Genes_n978.xlsx"), to=file.path("data", "L1000", "l1000_genes.xlsx"))
  ## read xls into data frame
  l1000.genes <- gdata::read.xls(xls=file.path("data", "L1000", "l1000_genes.xlsx"), sheet=1)
  l1000.genes <- l1000.genes[!is.na(l1000.genes[ , "Entrez.Gene.ID"]) & !duplicated(l1000.genes[ , "Entrez.Gene.ID"]), , drop=FALSE]
  rownames(l1000.genes) <- paste("geneid", as.character(l1000.genes[ , "Entrez.Gene.ID"]), sep="_")
  save(list=c("l1000.genes"), compress=TRUE, file=myfn)
} else { load(myfn) }


if(!file.exists("cdrug2_log.txt")) {
  steps <- c("cdrug2_analysis_birtwistle", "cdrug2_normalization_cgp", "cdrug2_normalization_ccle", "cdrug2_format", "cdrug2_analysis_huang")
  progress.log <- cbind(steps, rep("...", length(steps)))
  dimnames(progress.log) <- list(paste("step", 0:(length(steps)-1), sep="."), c("script", "progress"))
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
} else {
    progress.log <- read.table(file=file.path("cdrug2_log.txt"), sep="\t", header=TRUE, stringsAsFactor=FALSE)
}

########################
## re-analysis of drug sensitivitry data

message("\n-----------------------------\n| Normalization of CGP data |\n-----------------------------")
if (progress.log["step.0", "progress"] != "done") {
  progress.log["step.0", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
  source(file.path("code", "cdrug2_analysis_birtwistle.R"))
  progress.log["step.0", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
}
message("\t-> DONE")

########################
## curation, annotation and normalization of CGP data

message("\n-----------------------------\n| Normalization of CGP data |\n-----------------------------")
if (progress.log["step.1", "progress"] != "done") {
  progress.log["step.1", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
  source(file.path("code", "cdrug2_normalization_cgp.R"))
  progress.log["step.1", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
}
message("\t-> DONE")

#######################
## curation, annotation and normalization of CCLE data

message("\n------------------------------\n| Normalization of CCLE data |\n------------------------------")
if (progress.log["step.2", "progress"] != "done") {
  progress.log["step.2", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
  source(file.path("code", "cdrug2_normalization_ccle.R"))
  progress.log["step.2", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
}
message("\t-> DONE")

#######################
## intersection of CGP and CCLE

message("\n-------------------------------------\n| Intersection between GGP and CCLE |\n-------------------------------------")
if (progress.log["step.3", "progress"] != "done") {
  progress.log["step.3", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
  source(file.path("code", "cdrug2_format.R"))
  progress.log["step.3", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
}
message("\t-> DONE")

########################
## script performing the correlation analyses at the level of gene expressions

message("\n-------------------------------------------------------------------\n| Additional analysis regarding consistency between GGP and CCLE |\n-------------------------------------------------------------------")
if (progress.log["step.4", "progress"] != "done") {
  progress.log["step.4", "progress"] <- "in progress"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
  source(file.path("code", "cdrug2_analysis_huang.R"))
  progress.log["step.4", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
}
message("\t-> DONE")


## save session info
write(toLatex(sessionInfo(), locale = FALSE), file="sessionInfoR.tex", append=FALSE)

