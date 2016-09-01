options(stringsAsFactors=FALSE)

## load functions
source(file.path("code", "cdrug_foo.R"))

saveres <- file.path("output_birtwistle")

### install all the libraries at once
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("VennDiagram", "Hmisc", "xtable", "RColorBrewer", "pROC", "Biobase", "genefu", "PharmacoGx", "xlsx"))
library(VennDiagram)
library(Hmisc)
library(xtable)
library(RColorBrewer)
library(pROC)
library(Biobase)
library(genefu)
library(xlsx)
# install the latest devel version of the PharmacoGx package
# library(devtools)
# devtools::install_github("bhklab/PharmacoGx", ref="master", lib="/mnt/work1/users/bhklab/Rlib/")
library(PharmacoGx)s

### cell lines used in Haibe-Kains et al, Nature 2013
nature2013.common.cellines <- read.csv(file=file.path("data", "HaibeKains_Nature_2013_common_cellines.csv"), stringsAsFactors=FALSE)[ , 1]

### llist of known targets and biomarkers
known.biomarkers <- read.csv(file.path("data", "known_biomarkers.csv"), stringsAsFactors=FALSE)

### global parameters

nbcore <- 4

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

confine.analyses.to.nature.common.cell.lines <- FALSE

cor.method <- "pearson"


#################################################
## read irtwistle's results
#################################################



### read drug identifiers
drug.names <- read.xlsx(file=file.path("data", "Birtwistle", "DrugNames.xlsx"), sheetIndex=1, stringsAsFactors=FALSE)
rownames(drug.names) <- drug.names[ , "PharmacoGx.Name"]

### read cell line-drug identifiers
expn <- read.xlsx(file=file.path("data", "Birtwistle", "CelllineNamesInds.xlsx"), sheetIndex=1, stringsAsFactors=FALSE)

celline.names <- cbind("Cell.Line.Name"=sort(unique(expn[ , "Cell.Line.Name"])), "PharmacoGx.Name"=NA)
pgx.names <- sort(unique(expn[ , "Cell.Line.Name"])) ### UPDATE!!!
celline.names[ , "PharmacoGx.Name"] <- pgx.names
rownames(celline.names) <- celline.names[ , "PharmacoGx.Name"]

### update experiment names
exp.names <- cbind("Experiment.ID"=expn[ , "Cell.Line..Drug.Pair.ID"], t(sapply(expn[ , "Cell.Line..Drug.Pair.ID"], function (x, exps, cells, drugs) {
        iix <- which(exps[ , "Cell.Line..Drug.Pair.ID"] == x)[1]
        cc <- exps[iix, "Cell.Line.Name"]
        cc <- cells[which(cells[ , "Cell.Line.Name"] == cc)[1], "PharmacoGx.Name"]
        dd <- exps[iix, "Drug.ID.number"]
        dd <- drugs[which(drugs[ , "Drug.ID"] == dd)[1], "PharmacoGx.Name"]
        return (c("Cell.Line.Name"=cc, "Drug.Name"=dd))       
    }, exps=expn, cells=celline.names, drugs=drug.names)))
    

### read manual curations and map them to experiments
study.id <- c("CCLE"=1, "CGP"=2)
mc <- dir(path=file.path("data", "Birtwistle"), pattern="IP", full.names=TRUE)
mc.all <- NULL
for (i in 1:length(mc)) {
    tt <- read.csv(mc[i], header=FALSE, stringsAsFactors=FALSE)
    colnames(tt) <- c("Study.ID", "Experiment.ID", "Classification")
    ### Study names
    ss <- factor(tt[ , "Study.ID"], levels=c(1, 2))
    levels(ss) <- names(study.id)
    ss <- as.character(ss)
    ### drug names
    ee <- tt[ , "Experiment.ID"]
   
   
   
    oo <- factor(tt[ , "Classification"], levels=c(1, 2))
    levels(oo) <- c("sensitive", "insensitive")
    oo <- as.character(oo)
    
    mc.cgp <- cbind("Experiment.ID"=tt[ss == "CGP", "Experiment.ID"], ee[ss == "CGP", , drop=FALSE], "Classification"=oo[ss == "CGP"])
    mc.ccle <- cbind("Experiment.ID"=tt[ss == "CCLE", "Experiment.ID"], ee[ss == "CCLE", , drop=FALSE], oo[ss == "CCLE"])
    
    mc.all <- c(mc.all, list(cbind(ss, ))
}


#################################################
## get pharmacogenomic datasets
#################################################


myfn <- file.path(saveres, "data_cgp_ccle_strict.RData")
if (!file.exists(myfn)) {
  ### download curated pharmacogenomic data from CGP and CCLE
  CGP <- PharmacoGx::downloadPSet("GDSC", saveDir=file.path(saveres, "PSets"))
  CGP@annotation$name <- "CGP"
  CCLE <- PharmacoGx::downloadPSet("CCLE", saveDir=file.path(saveres, "PSets")) 
  ### common concentration range
  if(confine.analyses.to.nature.common.cell.lines) {
    common <- PharmacoGx::intersectPSet(pSets = list("CCLE"=CCLE, "CGP"=CGP), intersectOn = c("cell.lines", "drugs"), cells=nature2013.common.cellines, strictIntersect=TRUE)  
  } else { 
    common <- PharmacoGx::intersectPSet(pSets = list("CCLE"=CCLE, "CGP"=CGP), intersectOn = c("cell.lines", "drugs"), strictIntersect=TRUE)
  }
  ### CGP
  cgp.auc <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CGP, sensitivity.measure="auc_recomputed", summary.stat="median")
  cgp.ic50 <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CGP, sensitivity.measure="ic50_recomputed", summary.stat="median")
  cgp.slope <- apply(X=common$CGP@sensitivity$raw, MARGIN=1, FUN=function (x) {
    return (PharmacoGx::computeSlope(concentration=x[ , "Dose"], viability=x[ , "Viability"]))
  })
  common$CGP@sensitivity$profiles <- cbind(common$CGP@sensitivity$profiles, "slope_recomputed"=cgp.slope)
  cgp.slope <- summarizeSensitivityProfiles(pSet=common$CGP, sensitivity.measure="slope_recomputed", summary.stat="median")
  ### CCLE
  ccle.auc <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="auc_recomputed", summary.stat="median")
  ccle.ic50 <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="ic50_recomputed", summary.stat="median")
  ccle.slope <- apply(X=common$CCLE@sensitivity$raw, MARGIN=1, FUN=function (x) {
    return (PharmacoGx::computeSlope(concentration=x[ , "Dose"], viability=x[ , "Viability"]))
  })
  common$CCLE@sensitivity$profiles <- cbind(common$CCLE@sensitivity$profiles, "slope_recomputed"=ccle.slope)
  ccle.slope <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="slope_recomputed", summary.stat="median")
  ### note that ic50 and auc recomputed using a unified pipeline could be selected by using ic50_recomputed or auc_recomputed
  save(list=c("CGP", "CCLE", "common", "cgp.auc", "cgp.ic50", "cgp.slope", "ccle.auc", "ccle.ic50", "ccle.slope"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}


#################################################
### Supplementary Figure 1
### scatter plot with red points for missed cell lines in nature study
#################################################

pdf(file.path(saveres, "cgp_ccle_scatterplot_auc_missed_red_points.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
i <- 0
for(drugn in drugNames(common$CGP)) {
  i <- i + 1
  xx <- cgp.auc[drugn, ]
  yy <- ccle.auc[drugn, ]
  xxlim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
  yylim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
  nnn <- sum(complete.cases(xx, yy))
  par(mar=c(4, 4, 3, 1) + 0.1)
  mycol <- rep(blues9[7], length(xx))
  names(mycol) <- names(xx)
  mycol[!is.element(names(mycol), nature2013.common.cellines)] <- "red3"
  myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 11, "AUC (CGP)", ""), ylab=ifelse((i %% 4) == 1, "AUC (CCLE)", ""), main=drugn, xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75, col=mycol)
  abline(a=0, b=1, col="black")
}
dev.off()

#################################################
### xrug-dose response curves for PACLITAXEL:JVM−3 and NILOTINIB:EM−2
#################################################