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
library(psych)
# install the latest devel version of the PharmacoGx package
# library(devtools)
# devtools::install_github("bhklab/PharmacoGx", ref="master")
library(PharmacoGx)

### global parameters

nbcore <- 4
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

cor.method <- "pearson"


#################################################
## get pharmacogenomic datasets
#################################################

### download curated pharmacogenomic data for CGP and CCLE
CGP <- PharmacoGx::downloadPSet("GDSC", saveDir=file.path(saveres, "PSets"))
CGP@annotation$name <- "CGP"
CCLE <- PharmacoGx::downloadPSet("CCLE", saveDir=file.path(saveres, "PSets")) 

myfn <- file.path(saveres, "data_cgp_ccle.RData")
if (!file.exists(myfn)) {
  ### full concentration range
  common <- PharmacoGx::intersectPSet(pSets = list("CCLE"=CCLE, "CGP"=CGP), intersectOn = c("cell.lines", "drugs"), strictIntersect=TRUE, nthread=nbcore)
  ### CGP
  cgp.auc <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CGP, sensitivity.measure="auc_recomputed", summary.stat="median")
  # cgp.ic50 <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CGP, sensitivity.measure="ic50_recomputed", summary.stat="median")
  cgp.slope <- apply(X=common$CGP@sensitivity$raw, MARGIN=1, FUN=function (x) {
    return (PharmacoGx::computeSlope(concentration=x[ , "Dose"], viability=x[ , "Viability"], verbose=FALSE))
  })
  common$CGP@sensitivity$profiles <- cbind(common$CGP@sensitivity$profiles, "slope_recomputed"=cgp.slope)
  cgp.slope <- summarizeSensitivityProfiles(pSet=common$CGP, sensitivity.measure="slope_recomputed", summary.stat="median")
  ### CCLE
  ccle.auc <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="auc_recomputed", summary.stat="median")
  # ccle.ic50 <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="ic50_recomputed", summary.stat="median")
  ccle.slope <- apply(X=common$CCLE@sensitivity$raw, MARGIN=1, FUN=function (x) {
    return (PharmacoGx::computeSlope(concentration=x[ , "Dose"], viability=x[ , "Viability"], verbose=FALSE))
  })
  common$CCLE@sensitivity$profiles <- cbind(common$CCLE@sensitivity$profiles, "slope_recomputed"=ccle.slope)
  ccle.slope <- PharmacoGx::summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="slope_recomputed", summary.stat="median")
  ### note that ic50 and auc recomputed using a unified pipeline could be selected by using ic50_recomputed or auc_recomputed
  save(list=c("common", "cgp.auc", "cgp.slope", "ccle.auc", "ccle.slope"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}

myfn <- file.path(saveres, "data_cgp_ccle_strict.RData")
if (!file.exists(myfn)) {
  ### common concentration range
  commons <- PharmacoGx::intersectPSet(pSets = list("CCLE"=CCLE, "CGP"=CGP), intersectOn = c("cell.lines", "drugs", "concentrations"), strictIntersect=TRUE, nthread=nbcore)
  ### CGP
  cgp.aucs <- PharmacoGx::summarizeSensitivityProfiles(pSet=commons$CGP, sensitivity.measure="auc_recomputed_star", summary.stat="median")
  # cgp.ic50s <- PharmacoGx::summarizeSensitivityProfiles(pSet=commons$CGP, sensitivity.measure="ic50_recomputed", summary.stat="median")
  cgp.slopes <- apply(X=commons$CGP@sensitivity$raw, MARGIN=1, FUN=function (x) {
    return (PharmacoGx::computeSlope(concentration=x[ , "Dose"], viability=x[ , "Viability"], verbose=FALSE))
  })
  commons$CGP@sensitivity$profiles <- cbind(commons$CGP@sensitivity$profiles, "slope_recomputed_star"=cgp.slopes)
  cgp.slopes <- summarizeSensitivityProfiles(pSet=commons$CGP, sensitivity.measure="slope_recomputed_star", summary.stat="median")
  ### CCLE
  ccle.aucs <- PharmacoGx::summarizeSensitivityProfiles(pSet=commons$CCLE, sensitivity.measure="auc_recomputed", summary.stat="median")
  # ccle.ic50s <- PharmacoGx::summarizeSensitivityProfiles(pSet=commons$CCLE, sensitivity.measure="ic50_recomputed", summary.stat="median")
  ccle.slopes <- apply(X=commons$CCLE@sensitivity$raw, MARGIN=1, FUN=function (x) {
    return (PharmacoGx::computeSlope(conc=x[ , "Dose"], viability=x[ , "Viability"], verbose=FALSE))
  })
  commons$CCLE@sensitivity$profiles <- cbind(commons$CCLE@sensitivity$profiles, "slope_recomputed_star"=ccle.slopes)
  ccle.slopes <- summarizeSensitivityProfiles(pSet=commons$CCLE, sensitivity.measure="slope_recomputed_star", summary.stat="median")
  ### note that ic50 and auc recomputed using a unified pipeline could be selected by using ic50_recomputed or auc_recomputed
  save(list=c("commons", "cgp.aucs", "cgp.slopes", "ccle.aucs", "ccle.slopes"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}

drugn <- rownames(drugInfo(common$CCLE))

#################################################
## read irtwistle's results
#################################################

myfn <- file.path(saveres, "data_birtwistle.RData")
if (!file.exists(myfn)) {
  ### read drug identifiers
  drug.names <- read.xlsx(file=file.path("data", "Birtwistle", "DrugNames.xlsx"), sheetIndex=1, stringsAsFactors=FALSE)
  rownames(drug.names) <- drug.names[ , "PharmacoGx.Name"]

  ### read cell line-drug identifiers
  expn <- read.xlsx(file=file.path("data", "Birtwistle", "CelllineNamesInds.xlsx"), sheetIndex=1, stringsAsFactors=FALSE)

  ## read cell line names
  celline.names <- read.xlsx(file=file.path("data", "Birtwistle", "CellLineNames.xlsx"), sheetIndex=1, stringsAsFactors=FALSE)
  rownames(celline.names) <- celline.names[ , "PharmacoGx.Name"]

  ### update experiment names
  exp.names <- cbind("Experiment.ID"=paste("exp", expn[ , "Cell.Line..Drug.Pair.ID"], sep="_"), t(sapply(expn[ , "Cell.Line..Drug.Pair.ID"], function (x, exps, cells, drugs) {
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
  mc.cgp <- mc.ccle <- data.frame(matrix(NA, nrow=nrow(exp.names), ncol=ncol(exp.names) + length(mc), dimnames=list(exp.names[ , "Experiment.ID"], c(colnames(exp.names), paste("IP", 1:length(mc), sep="")))))
  mc.ccle[ , colnames(exp.names)] <- mc.cgp[ , colnames(exp.names)] <- exp.names
  ## each curator
  for (i in 1:length(mc)) {
      tt <- read.csv(mc[i], header=FALSE, stringsAsFactors=FALSE)
      colnames(tt) <- c("Study.ID", "Experiment.ID", "Classification")
      ### CGP
      oo <- tt[!is.na(tt[ , "Study.ID"]) & tt[ , "Study.ID"] == study.id["CGP"], "Classification"]
      names(oo) <- paste("exp", tt[!is.na(tt[ , "Study.ID"]) & tt[ , "Study.ID"] == study.id["CGP"], "Experiment.ID"], sep="_")
      mc.cgp[ , paste("IP", i, sep="")] <- oo[rownames(mc.cgp)]
      ### CCLE
      oo <- tt[!is.na(tt[ , "Study.ID"]) & tt[ , "Study.ID"] == study.id["CCLE"], "Classification"]
      names(oo) <- paste("exp", tt[!is.na(tt[ , "Study.ID"]) & tt[ , "Study.ID"] == study.id["CCLE"], "Experiment.ID"], sep="_")
      mc.ccle[ , paste("IP", i, sep="")] <- oo[rownames(mc.cgp)]
  }
  ### consensus of curators
  ### CGP
  tt <- apply(mc.cgp[ , paste("IP", 1:length(mc), sep="")], 1, function (x) {
    tt <- table(x)
    if(length(tt) == 2) {
      tt <- tt[c("1", "2")]
    }
    tt <- names(tt)[which.max(tt)]
    return (tt)
  })
  mc.cgp <- data.frame(mc.cgp, "IPALL"=tt)
  ### CCLE
  tt <- apply(mc.ccle[ , paste("IP", 1:length(mc), sep="")], 1, function (x) {
    tt <- table(x)
    if(length(tt) == 2) {
      tt <- tt[c("1", "2")]
    }
    tt <- names(tt)[which.max(tt)]
    return (tt)
  })
  mc.ccle <- data.frame(mc.ccle, "IPALL"=tt)
  save(list=c("mc.cgp", "mc.ccle"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}

### transform list of experiments in matrix format
### CGP
cgp.mcall <- matrix(NA, nrow=nrow(cgp.slope), ncol=ncol(cgp.slope), dimnames=dimnames(cgp.slope))
for (i in 1:nrow(mc.cgp)) {
  cgp.mcall[mc.cgp[i, "Drug.Name"], mc.cgp[i, "Cell.Line.Name"]] <- mc.cgp[i, "IPALL"]
}
### CCLE
ccle.mcall <- matrix(NA, nrow=nrow(ccle.slope), ncol=ncol(ccle.slope), dimnames=dimnames(ccle.slope))
for (i in 1:nrow(mc.ccle)) {
  ccle.mcall[mc.ccle[i, "Drug.Name"], mc.ccle[i, "Cell.Line.Name"]] <- mc.ccle[i, "IPALL"]
}
### 1 = sensitive; 2 = insensitive


#################################################
### scatter plot with red points for missed cell lines in nature study
#################################################

pdf(file.path(saveres, "cgp_ccle_pool_drugs_cor.pdf"), height=7, width=7)
par(mfrow=c(2, 2), cex=0.8, las=1)
### AUC
xx <- as.numeric(cgp.auc)
yy <- as.numeric(ccle.auc)
xxlim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
yylim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
nnn <- sum(complete.cases(xx, yy))
par(mar=c(4, 4, 3, 1) + 0.1)
mycol <- rep(blues9[7], length(xx))
names(mycol) <- names(xx)
myScatterPlot(x=xx, y=yy, xlab="AUC (CGP)", ylab="AUC (CCLE)", main="", xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.25, col=mycol)
abline(a=0, b=1, col="black")
cc.auc <- cor.test(xx, yy, method=cor.method, use="complete.obs")
nn.auc <- sum(complete.cases(xx, yy))
legend("bottomright", legend=sprintf("rho = %.2g", cc.auc$estimate), bty="n", cex=1.25)

### AUC_s
xx <- as.numeric(cgp.aucs)
yy <- as.numeric(ccle.aucs)
xxlim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
yylim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
nnn <- sum(complete.cases(xx, yy))
par(mar=c(4, 4, 3, 1) + 0.1)
mycol <- rep(blues9[7], length(xx))
names(mycol) <- names(xx)
myScatterPlot(x=xx, y=yy, xlab="AUC_s (CGP)", ylab="AUC_s (CCLE)", main="", xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.25, col=mycol)
abline(a=0, b=1, col="black")
cc.aucs <- cor.test(xx, yy, method=cor.method, use="complete.obs")
nn.aucs <- sum(complete.cases(xx, yy))
legend("bottomright", legend=sprintf("rho = %.2g", cc.aucs$estimate), bty="n", cex=1.25)
### test for difference in correlation
cor.diff <- psych::r.test(n=nn.auc, n2=nn.aucs, r12=cc.auc$estimate, r34=cc.aucs$estimate)
message(sprintf("rho AUC_s = %.2g vs rho AUC = %.2g -> p = %.1E", cc.aucs$estimate, cc.auc$estimate, cor.diff$p))

### m
xx <- as.numeric(cgp.slope)
yy <- as.numeric(ccle.slope)
xxlim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
yylim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
nnn <- sum(complete.cases(xx, yy))
par(mar=c(4, 4, 3, 1) + 0.1)
mycol <- rep(blues9[7], length(xx))
names(mycol) <- names(xx)
myScatterPlot(x=xx, y=yy, xlab="m (CGP)", ylab="m (CCLE)", main="", xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.25, col=mycol)
abline(a=0, b=1, col="black")
cc.m <- cor.test(xx, yy, method=cor.method, use="complete.obs")
nn.m <- sum(complete.cases(xx, yy))
legend("bottomright", legend=sprintf("rho = %.2g", cc.m$estimate), bty="n", cex=1.25)
### m_s
xx <- as.numeric(cgp.slopes)
yy <- as.numeric(ccle.slopes)
xxlim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
yylim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
nnn <- sum(complete.cases(xx, yy))
par(mar=c(4, 4, 3, 1) + 0.1)
mycol <- rep(blues9[7], length(xx))
names(mycol) <- names(xx)
myScatterPlot(x=xx, y=yy, xlab="m_s (CGP)", ylab="m_s (CCLE)", main="", xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.25, col=mycol)
abline(a=0, b=1, col="black")
cc.ms <- cor.test(xx, yy, method=cor.method, use="complete.obs")
nn.ms <- sum(complete.cases(xx, yy))
legend("bottomright", legend=sprintf("rho = %.2g", cc.ms$estimate), bty="n", cex=1.25)
### test for difference in correlation
cor.diff <- psych::r.test(n=nn.m, n2=nn.ms, r12=cc.m$estimate, r34=cc.ms$estimate)
message(sprintf("rho m_s = %.2g vs rho m = %.2g -> p = %.1E", cc.ms$estimate, cc.m$estimate, cor.diff$p))
dev.off()


pdf(file.path(saveres, "cgp_ccle_each_drug_cor_auc.pdf"), height=5, width=10)
par(mar=c(7, 4, 1, 2) + 0.1)
cc.auc.all <- matrix(NA, nrow=2, ncol=length(drugn), dimnames=list(c("AUC", "AUC_s"), drugn))
for (i in 1:length(drugn)) {
  ### AUC
  xx <- cgp.auc[drugn[i], ]
  yy <- ccle.auc[drugn[i], ]
  cc <- cor(xx, yy, method=cor.method, use="complete.obs")
  cc.auc.all["AUC", drugn[i]] <- cc
  ### AUC_s
  xx <- cgp.aucs[drugn[i], ]
  yy <- ccle.aucs[drugn[i], ]
  cc <- cor(xx, yy, method=cor.method, use="complete.obs")
  cc.auc.all["AUC_s", drugn[i]] <- cc
}
bp <- barplot(cc.auc.all, las=2, beside=TRUE, space=c(0.1, 1), ylim=c(0, 1), ylab=expression(rho), cex.lab=1.5, col=c("grey20", "grey70"))
legend("topright", legend=c("AUC", "AUC_s"), col=c("grey20", "grey70"), bty="n", pch=15, cex=1.25)
dev.off()
wilcox.test(cc.auc.all["AUC_s", ], cc.auc.all["AUC", ], paired=TRUE, alternative="greater", exact=FALSE)

pdf(file.path(saveres, "cgp_ccle_each_drug_cor_slope.pdf"), height=5, width=10)
par(mar=c(7, 4, 1, 2) + 0.1)
cc.slope.all <- matrix(NA, nrow=2, ncol=length(drugn), dimnames=list(c("m", "m_s"), drugn))
for (i in 1:length(drugn)) {
  ### AUC
  xx <- cgp.slope[drugn[i], ]
  yy <- ccle.slope[drugn[i], ]
  cc <- cor(xx, yy, method=cor.method, use="complete.obs")
  cc.slope.all["m", drugn[i]] <- cc
  ### AUC_s
  xx <- cgp.slopes[drugn[i], ]
  yy <- ccle.slopes[drugn[i], ]
  cc <- cor(xx, yy, method=cor.method, use="complete.obs")
  cc.slope.all["m_s", drugn[i]] <- cc
}
bp <- barplot(cc.slope.all, las=2, beside=TRUE, space=c(0.1, 1), ylim=c(0, 1), ylab=expression(rho), cex.lab=1.5, col=c("grey20", "grey70"))
legend("topright", legend=c("m", "m_s"), col=c("grey20", "grey70"), bty="n", pch=15, cex=1.25)
dev.off()
wilcox.test(cc.slope.all["m_s", ], cc.slope.all["m", ], paired=TRUE, alternative="greater", exact=FALSE)


### AUC
pdf(file.path(saveres, "cgp_ccle_scatterplot_aucs.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
i <- 0
for(drugn in drugNames(common$CGP)) {
  i <- i + 1
  xx <- cgp.aucs[drugn, ]
  yy <- ccle.aucs[drugn, ]
  xxlim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
  yylim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
  nnn <- sum(complete.cases(xx, yy))
  par(mar=c(4, 4, 3, 1) + 0.1)
  mycol <- rep(blues9[7], length(xx))
  names(mycol) <- names(xx)
  myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 11, "AUC_s (CGP)", ""), ylab=ifelse((i %% 4) == 1, "AUC_s (CCLE)", ""), main=drugn, xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.5, col=mycol)
  abline(a=0, b=1, col="black")
}
### pooled drugs
i <- i + 1
xx <- as.numeric(cgp.aucs)
yy <- as.numeric(ccle.aucs)
xxlim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
yylim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
nnn <- sum(complete.cases(xx, yy))
par(mar=c(4, 4, 3, 1) + 0.1)
mycol <- rep(blues9[7], length(xx))
names(mycol) <- names(xx)
myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 11, "AUC_s (CGP)", ""), ylab=ifelse((i %% 4) == 1, "AUC_s (CCLE)", ""), main="All drugs", xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.25, col=mycol)
abline(a=0, b=1, col="black")
dev.off()

### Slope
pdf(file.path(saveres, "cgp_ccle_scatterplot_ms.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
i <- 0
for(drugn in drugNames(common$CGP)) {
  i <- i + 1
  xx <- cgp.slopes[drugn, ]
  yy <- ccle.slopes[drugn, ]
  xxlim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
  yylim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
  nnn <- sum(complete.cases(xx, yy))
  par(mar=c(4, 4, 3, 1) + 0.1)
  mycol <- rep(blues9[7], length(xx))
  names(mycol) <- names(xx)
  myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 11, "m_s (CGP)", ""), ylab=ifelse((i %% 4) == 1, "m_s (CCLE)", ""), main=drugn, xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.5, col=mycol)
  abline(a=0, b=1, col="black")
}
### pooled drugs
i <- i + 1
xx <- as.numeric(cgp.slopes)
yy <- as.numeric(ccle.slopes)
xxlim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
yylim <- c(0, ceiling(max(xx, yy, na.rm=TRUE) * 10) / 10)
nnn <- sum(complete.cases(xx, yy))
par(mar=c(4, 4, 3, 1) + 0.1)
mycol <- rep(blues9[7], length(xx))
names(mycol) <- names(xx)
myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 11, "m_s (CGP)", ""), ylab=ifelse((i %% 4) == 1, "m_s (CCLE)", ""), main="All drugs", xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.25, col=mycol)
abline(a=0, b=1, col="black")
dev.off()

#################################################
### MCC overall and per drug
#################################################

mycol <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00')

### discreatization of AUC and Slope using 0.2 as cutoff
### CGP
cgp.auc2 <- apply(cgp.auc, 2, function (x, y) {
  x <- factor(x > y, levels=c(TRUE, FALSE))
  levels(x) <- c("1", "2")
  return(x)
}, y=0.2)
cgp.aucs2 <- apply(cgp.aucs, 2, function (x, y) {
  x <- factor(x > y, levels=c(TRUE, FALSE))
  levels(x) <- c("1", "2")
  return(x)
}, y=0.2)
cgp.slope2 <- apply(cgp.slope, 2, function (x, y) {
  x <- factor(x > y, levels=c(TRUE, FALSE))
  levels(x) <- c("1", "2")
  return(x)
}, y=0.2)
cgp.slopes2 <- apply(cgp.slopes, 2, function (x, y) {
  x <- factor(x > y, levels=c(TRUE, FALSE))
  levels(x) <- c("1", "2")
  return(x)
}, y=0.2)
### CCLE
ccle.auc2 <- apply(ccle.auc, 2, function (x, y) {
  x <- factor(x > y, levels=c(TRUE, FALSE))
  levels(x) <- c("1", "2")
  return(x)
}, y=0.2)
ccle.aucs2 <- apply(ccle.aucs, 2, function (x, y) {
  x <- factor(x > y, levels=c(TRUE, FALSE))
  levels(x) <- c("1", "2")
  return(x)
}, y=0.2)
ccle.slope2 <- apply(ccle.slope, 2, function (x, y) {
  x <- factor(x > y, levels=c(TRUE, FALSE))
  levels(x) <- c("1", "2")
  return(x)
}, y=0.2)
ccle.slopes2 <- apply(ccle.slopes, 2, function (x, y) {
  x <- factor(x > y, levels=c(TRUE, FALSE))
  levels(x) <- c("1", "2")
  return(x)
}, y=0.2)

### example of why MCC is a more approriate mesure of consistency
# i <- 2; print(drugn[i]); tt <- table("CGP"=cgp.auc2[drugn[i],], "CCLE"=ccle.auc2[drugn[i],]) ; print(tt); PharmacoGx:::.mcc(ct=tt); sum(diag(tt)) / sum(tt)

### pooled drugs
mcc.pool <- c(
  "AUC"=mcc(x=factor(as.character(cgp.auc2), levels=c("1", "2")), y=factor(as.character(ccle.auc2), , levels=c("1", "2")), nperm=0)$estimate,
  "AUC_s"=mcc(x=factor(as.character(cgp.aucs2), levels=c("1", "2")), y=factor(as.character(ccle.aucs2), , levels=c("1", "2")), nperm=0)$estimate,
  "m"=mcc(x=factor(as.character(cgp.slope2), levels=c("1", "2")), y=factor(as.character(ccle.slope2), , levels=c("1", "2")), nperm=0)$estimate,
  "m_s"=mcc(x=factor(as.character(cgp.slopes2), levels=c("1", "2")), y=factor(as.character(ccle.slopes2), , levels=c("1", "2")), nperm=0)$estimate,
  "Manual curation"=mcc(x=factor(as.character(cgp.mcall), levels=c("1", "2")), y=factor(as.character(ccle.mcall), , levels=c("1", "2")), nperm=0)$estimate
)

pdf(file.path(saveres, "mcc_all_drugs.pdf"), width=7, height=5)
par(mar=c(7, 4, 1, 2) + 0.1)
barplot(mcc.pool, col=mycol, ylim=c(0, 1), las=2, space=0.5, main="", ylab="MCC")
dev.off()
# sapply(mcc.pool, function (x) { return (x$p.value)})

### for each drug separately
mcc.drugs <- NULL
for (i in 1:length(drugn)) {
  mcc.drug <- c(
    "AUC"=mcc(x=factor(as.character(cgp.auc2[drugn[i], ]), levels=c("1", "2")), y=factor(as.character(ccle.auc2[drugn[i], ]), , levels=c("1", "2")), nperm=0)$estimate,
    "AUC_s"=mcc(x=factor(as.character(cgp.aucs2[drugn[i], ]), levels=c("1", "2")), y=factor(as.character(ccle.aucs2[drugn[i], ]), , levels=c("1", "2")), nperm=0)$estimate,
    "m"=mcc(x=factor(as.character(cgp.slope2[drugn[i], ]), levels=c("1", "2")), y=factor(as.character(ccle.slope2[drugn[i], ]), , levels=c("1", "2")), nperm=0)$estimate,
    "m_s"=mcc(x=factor(as.character(cgp.slopes2[drugn[i], ]), levels=c("1", "2")), y=factor(as.character(ccle.slopes2[drugn[i], ]), , levels=c("1", "2")), nperm=0)$estimate,
    "Manual curation"=mcc(x=factor(as.character(cgp.mcall[drugn[i], ]), levels=c("1", "2")), y=factor(as.character(ccle.mcall[drugn[i], ]), , levels=c("1", "2")), nperm=0)$estimate
  )
  mcc.drugs <- rbind(mcc.drugs, mcc.drug)
}
rownames(mcc.drugs) <- drugn

pdf(file.path(saveres, "mcc_each_drug.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, mar=c(7, 4, 1, 2) + 0.1)
i <- 0
for(i in 1:length(drugn)) {
  barplot(mcc.drugs[drugn[i], ], col=mycol, ylim=c(0, 1), space=0.5, main=drugn[i], ylab="MCC")
}
barplot(mcc.pool, col=mycol, ylim=c(0, 1), space=0.5, main="All drugs", ylab="MCC", las=2)

# axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
# text(x=apply(mp, 2, mean) + 1.45, y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, font=2)
dev.off()

### Does manual curation yield higher MCC than AUC_s
wilcox.test(x=mcc.drugs[ , "Manual curation"], y=mcc.drugs[ , "AUC_s"], alternative="greater", paired=TRUE, exact=FALSE)

### Does manual curation yield higher MCC than m_s
wilcox.test(x=mcc.drugs[ , "Manual curation"], y=mcc.drugs[ , "m_s"], alternative="greater", paired=TRUE, exact=FALSE)



#################################################
### drug-dose response curves for PACLITAXEL:JVM−3 and NILOTINIB:EM−2
#################################################

drugDoseResponseCurve(drug="paclitaxel", cellline="JVM-3", pSets=list(CGP, CCLE))
drugDoseResponseCurve(drug="Nilotinib", cellline="EM-2", pSets=list(CGP, CCLE))

## save session info
write(toLatex(sessionInfo(), locale=FALSE), file=file.path(saveres, "sessionInfoR.tex"), append=FALSE)

### end
