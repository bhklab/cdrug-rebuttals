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
nbcore <- 30
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

## max fdr, threshold used to identify genes with significant concordance index
myfdr <- 0.20

## method to estimate the association gene-drug, controlled for tissue type
# genedrugm <- c("lm", "cindex")
genedrugm <- "lm"

## consistency between associations
concordance.method <- "spearman"
fdr.cuts <- c(1, 0.5, 0.2, 0.05, 0.01)
quantile.cuts <- c(1, 0.5, 0.2, 0.05, 0.01)

## additional functions
source(file.path("code", "cdrug2_foo.R"))

if(!file.exists("cdrug2_log.txt")) {
  steps <- c("cdrug2_normalization_cgp", "cdrug2_normalization_ccle", "cdrug2_format", "cdrug2_analysis")
  progress.log <- cbind(steps, rep("...", length(steps)))
  dimnames(progress.log) <- list(paste("step", 1:length(steps), sep="."), c("script", "progress"))
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
} else {
    progress.log <- read.table(file=file.path("cdrug2_log.txt"), sep="\t", header=TRUE, stringsAsFactor=FALSE)
}

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

# #######################
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
  source(file.path("code", "cdrug2_analysis.R"))
  progress.log["step.4", "progress"] <- "done"
  write.table(progress.log, sep="\t", row.names=TRUE, col.names=TRUE, file=file.path("cdrug2_log.txt"), quote=FALSE)
}
message("\t-> DONE")


## save session info
write(toLatex(sessionInfo(), locale = FALSE), file="sessionInfoR.tex", append=FALSE)

