options(stringsAsFactors=FALSE)

## load functions
source(file.path("code", "cdrug_foo.R"))

saveres <- file.path("output_aittokallio")

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
# library(xlsx)
library(abind)
library(psych)
# install the latest devel version of the PharmacoGx package
# library(devtools)
# devtools::install_github("bhklab/PharmacoGx", ref="master")
library(PharmacoGx)

### global parameters

nbcore <- 16
availcore <- parallel::detectCores()
if (is.null(nbcore) || nbcore > availcore) { nbcore <- availcore }
options("mc.cores"=nbcore)

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

cor.method <- "pearson"


#################################################
## get pharmacogenomic datasets
#################################################

### download curated pharmacogenomic data for CGP and CCLE

### CGP
CGP <- PharmacoGx::downloadPSet("GDSC", saveDir=file.path(saveres, "PSets"))
CGP@annotation$name <- "CGP"

### CCLE
CCLE <- PharmacoGx::downloadPSet("CCLE", saveDir=file.path(saveres, "PSets")) 

### FIMM
load("/mnt/work1/users/bhklab/Projects/PharmacoGxTest/PSets/FIMM.RData")

message("Intersection of PSets")
myfn <- file.path(saveres, "data_common_cgp_ccle_fimm.RData")
if (!file.exists(myfn)) {
 
  ### full concentration range
  message("")
  message("\tFull concentration range")
 
  message("\t\tCCLE vs FIMM")
  common.ccle.fimm <- PharmacoGx::intersectPSet(pSets = list("CCLE"=CCLE, "FIMM"=FIMM), intersectOn = c("cell.lines", "drugs"), strictIntersect=TRUE, nthread=nbcore)
 
  message("\t\tCGP vs FIMM")
  common.cgp.fimm <- PharmacoGx::intersectPSet(pSets = list("CGP"=CGP, "FIMM"=FIMM), intersectOn = c("cell.lines", "drugs"), strictIntersect=TRUE, nthread=nbcore)
  
  message("\t\tCGP vs CCLE")
  common.cgp.ccle <- PharmacoGx::intersectPSet(pSets = list("CGP"=CGP, "CCLE"=CCLE), intersectOn = c("cell.lines", "drugs"), strictIntersect=TRUE, nthread=nbcore)
  
  ### common concentration range
  message("")
  message("\tCommon concentration range")
  
  message("\t\tCCLE vs FIMM")
  ppl <- PharmacoGx::intersectPSet(pSets = list("CCLE"=CCLE, "FIMM"=FIMM), intersectOn = c("cell.lines", "drugs", "concentrations"), strictIntersect=TRUE, nthread=nbcore)
  for (i in 1:length(ppl)) {
    ppl[[i]]@sensitivity$profiles[ , grep("*_recomputed", colnames(ppl[[i]]@sensitivity$profiles))] <- NA
    ### recompute auc on common concentration range
    aa <- apply(X=ppl[[i]]@sensitivity$raw, MARGIN=1, FUN=function (x) {
        return (PharmacoGx::computeAUC(concentration=x[ , "Dose"], viability=x[ , "Viability"], verbose=FALSE))
      })
    ppl[[i]]@sensitivity$profiles <- cbind(ppl[[i]]@sensitivity$profiles, "auc_recomputed"=aa / 100)
  }
  commons.ccle.fimm <- ppl
  rm(ppl)
  
  message("\t\tCGP vs FIMM")
  ppl <- PharmacoGx::intersectPSet(pSets = list("CGP"=CGP, "FIMM"=FIMM), intersectOn = c("cell.lines", "drugs", "concentrations"), strictIntersect=TRUE, nthread=nbcore)
  for (i in 1:length(ppl)) {
    ppl[[i]]@sensitivity$profiles[ , grep("*_recomputed", colnames(ppl[[i]]@sensitivity$profiles))] <- NA
    ### recompute auc on common concentration range
    aa <- apply(X=ppl[[i]]@sensitivity$raw, MARGIN=1, FUN=function (x) {
        return (PharmacoGx::computeAUC(concentration=x[ , "Dose"], viability=x[ , "Viability"], verbose=FALSE))
      })
    ppl[[i]]@sensitivity$profiles <- cbind(ppl[[i]]@sensitivity$profiles, "auc_recomputed"=aa / 100)
  }
  commons.cgp.fimm <- ppl
  rm(ppl)
  
  message("\t\tCGP vs CCLE")
  ppl <- PharmacoGx::intersectPSet(pSets = list("CGP"=CGP, "CCLE"=CCLE), intersectOn = c("cell.lines", "drugs", "concentrations"), strictIntersect=TRUE, nthread=nbcore)
  for (i in 1:length(ppl)) {
    ppl[[i]]@sensitivity$profiles[ , grep("*_recomputed", colnames(ppl[[i]]@sensitivity$profiles))] <- NA
    ### recompute auc on common concentration range
    aa <- apply(X=ppl[[i]]@sensitivity$raw, MARGIN=1, FUN=function (x) {
        return (PharmacoGx::computeAUC(concentration=x[ , "Dose"], viability=x[ , "Viability"], verbose=FALSE))
      })
    ppl[[i]]@sensitivity$profiles <- cbind(ppl[[i]]@sensitivity$profiles, "auc_recomputed"=aa / 100)
  }
  commons.cgp.ccle <- ppl
  rm(ppl)
  
  comp <- c("CCLE.vs.FIMM"="common.ccle.fimm", "CGP.vs.FIMM"="common.cgp.fimm", "CGP.vs.CCLE"="common.cgp.ccle")
  comps <- c("CCLE.vs.FIMM"="commons.ccle.fimm", "CGP.vs.FIMM"="commons.cgp.fimm", "CGP.vs.CCLE"="commons.cgp.ccle"
  save(list=c(comp, comps, c("comp", "comps")), compress=TRUE, file=myfn)
} else {
  load(myfn)
}


#################################################
### boxplot for AUC unharmonised vs harmonised
#################################################

across.ccl <- between.ccl <- across.ccl.rank <- between.ccl.rank <-NULL

for (i in 1:length(comp)) {
  
  ## full concentration range
  auc1 <- PharmacoGx::summarizeSensitivityProfiles(pSet=get(comp[i])[[1]], sensitivity.measure="auc_recomputed", summary.stat="median")
  auc2 <- PharmacoGx::summarizeSensitivityProfiles(pSet=get(comp[i])[[2]], sensitivity.measure="auc_recomputed", summary.stat="median")
  auc2 <- auc2[rownames(auc1), colnames(auc1)]
  auc <- abind(auc1, auc2, along=3)

  ### pearson
  tt <- apply(auc, 1, function (x, y) {
    if (sum(complete.cases(x)) < 10) {
      rr <- NA
    } else {
      rr <- cor(x=x[ , 1], y=x[ , 2], method=y, use="complete.obs")
    }
    return(rr)
  }, y="pearson")
  tt2 <- apply(auc, 2, function (x, y) {
    if (sum(complete.cases(x)) < 10) {
      rr <- NA
    } else {
      rr <- cor(x=x[ , 1], y=x[ , 2], method=y, use="complete.obs")
    }
    return(rr)
  }, y="pearson")
  
  ### spearman
  tt.rank <- apply(auc, 1, function (x, y) {
    if (sum(complete.cases(x)) < 10) {
      rr <- NA
    } else {
      rr <- cor(x=x[ , 1], y=x[ , 2], method=y, use="complete.obs")
    }
    return(rr)
  }, y="spearman")
  tt2.rank <- apply(auc, 2, function (x, y) {
    if (sum(complete.cases(x)) < 10) {
      rr <- NA
    } else {
      rr <- cor(x=x[ , 1], y=x[ , 2], method=y, use="complete.obs")
    }
    return(rr)
  }, y="spearman")
  
  ## common concentration range
  auc1 <- PharmacoGx::summarizeSensitivityProfiles(pSet=get(comps[i])[[1]], sensitivity.measure="auc_recomputed", summary.stat="median")
  auc2 <- PharmacoGx::summarizeSensitivityProfiles(pSet=get(comps[i])[[2]], sensitivity.measure="auc_recomputed", summary.stat="median")
  auc2 <- auc2[rownames(auc1), colnames(auc1)]
  auc <- abind(auc1, auc2, along=3)
  
  ### pearson
  tts <- apply(auc, 1, function (x, y) {
    if (sum(complete.cases(x)) < 10) {
      rr <- NA
    } else {
      rr <- cor(x=x[ , 1], y=x[ , 2], method=y, use="complete.obs")
    }
    return(rr)
  }, y="pearson")
  tts2 <- apply(auc, 2, function (x, y) {
    if (sum(complete.cases(x)) < 10) {
      rr <- NA
    } else {
      rr <- cor(x=x[ , 1], y=x[ , 2], method=y, use="complete.obs")
    }
    return(rr)
  }, y="pearson")
  
  ### spearman
  tts.rank <- apply(auc, 1, function (x, y) {
    if (sum(complete.cases(x)) < 10) {
      rr <- NA
    } else {
      rr <- cor(x=x[ , 1], y=x[ , 2], method=y, use="complete.obs")
    }
    return(rr)
  }, y="spearman")
  tts2.rank <- apply(auc, 2, function (x, y) {
    if (sum(complete.cases(x)) < 10) {
      rr <- NA
    } else {
      rr <- cor(x=x[ , 1], y=x[ , 2], method=y, use="complete.obs")
    }
    return(rr)
  }, y="spearman")
  
  across.ccl <- c(across.ccl, list(tt), list(tts))
  across.ccl.rank <- c(across.ccl.rank, list(tt.rank), list(tts.rank))
  if (i < length(comp)) {
    across.ccl <- c(across.ccl, list(NA))
    across.ccl.rank <- c(across.ccl.rank, list(NA))
  }
    
  between.ccl <- c(between.ccl, list(tt2), list(tts2))
  between.ccl.rank <- c(between.ccl.rank, list(tt2.rank), list(tts2.rank))
  if (i < length(comps)) {
    between.ccl <- c(between.ccl, list(NA))
    between.ccl.rank <- c(between.ccl.rank, list(NA))
  }
}

transparency <- 0.4
mycol <- c(rgb(red=1, green=0, blue=0, alpha=transparency, maxColorValue=1), rgb(red=0, green=0, blue=1, alpha=transparency, maxColorValue=1), rgb(red=1, green=1, blue=1, alpha=transparency, maxColorValue=1))
yylim <- c(-0.3, 1)

pdf(file.path(saveres, "boxplot_auc.pdf"), width=10, height=10)
par(mfrow=c(2,2))

par(las=2, mar=c(4, 4, 4, 2) + 0.1, xaxt="n")
bp <- boxplot(across.ccl, outline=FALSE, ylim=yylim, ylab="Pearson correlation", main="Consistency across cell lines")
axis(1, at=0.5 + seq(1, length(bp$names), by=3), tick=TRUE, labels=TRUE)
text(x=1.5 + seq(1, length(bp$names), by=3), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=names(comp), srt=0, xpd=NA)
stripchart(across.ccl, vertical=TRUE, method = "jitter", jitter=0.2, add=TRUE, pch=20, col=mycol, cex=2)
legend("bottomleft", legend=c("Harmonised", "Un-harmonised"), col=mycol[1:2], pch=20, pt.cex=2, bty="n")

par(las=2, mar=c(4, 4, 4, 2) + 0.1, xaxt="n")
bp <- boxplot(between.ccl, outline=FALSE, ylim=yylim, ylab="Pearson correlation", main="Consistency between cell lines")
axis(1, at=0.5 + seq(1, length(bp$names), by=3), tick=TRUE, labels=TRUE)
text(x=1.5 + seq(1, length(bp$names), by=3), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=names(comp), srt=0, xpd=NA)
stripchart(between.ccl, vertical=TRUE, method = "jitter", jitter=0.2, add=TRUE, pch=20, col=mycol, cex=2)
legend("bottomleft", legend=c("Harmonised", "Un-harmonised"), col=mycol[1:2], pch=20, pt.cex=2, bty="n")

par(las=2, mar=c(4, 4, 4, 2) + 0.1, xaxt="n")
bp <- boxplot(across.ccl.rank, outline=FALSE, ylim=yylim, ylab="Spearman correlation", main="Consistency across cell lines")
axis(1, at=1.5 + seq(1, length(bp$names), by=3), tick=TRUE, labels=TRUE)
text(x=1.5 + seq(1, length(bp$names), by=3), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=names(comp), srt=0, xpd=NA)
stripchart(across.ccl.rank, vertical=TRUE, method = "jitter", jitter=0.2, add=TRUE, pch=20, col=mycol, cex=2)
legend("bottomleft", legend=c("Harmonised", "Un-harmonised"), col=mycol[1:2], pch=20, pt.cex=2, bty="n")

par(las=2, mar=c(4, 4, 4, 2) + 0.1, xaxt="n")
bp <- boxplot(between.ccl.rank, outline=FALSE, ylim=yylim, ylab="Spearman correlation", main="Consistency between cell lines")
axis(1, at=1.5 + seq(1, length(bp$names), by=3), tick=TRUE, labels=TRUE)
text(x=1.5 + seq(1, length(bp$names), by=3), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=names(comp), srt=0, xpd=NA)
stripchart(between.ccl.rank, vertical=TRUE, method = "jitter", jitter=0.2, add=TRUE, pch=20, col=mycol, cex=2)
legend("bottomleft", legend=c("Harmonised", "Un-harmonised"), col=mycol[1:2], pch=20, pt.cex=2, bty="n")

dev.off()


#################################################
### box plots, correlation, across vs between
#################################################


myfn <- file.path(saveres, "data_cgp_ccle.RData")
if (!file.exists(myfn)) {
  
  ### common affymetrix probes
  common.features <- intersect(rownames(featureInfo(common.cgp.ccle$CCLE, "rna"))[featureInfo(common.cgp.ccle$CCLE, "rna")[,"BEST"]==T], rownames(featureInfo(common.cgp.ccle$CGP, "rna"))[featureInfo(common.cgp.ccle$CGP, "rna")[,"BEST"]==T])
  ### CGP
  cgp.ge <- PharmacoGx::summarizeMolecularProfiles(pSet=common.cgp.ccle$CGP, mDataType="rna", summary.stat="median")
  cgp.ge <- cgp.ge[common.features,]
  cgp.auc <- PharmacoGx::summarizeSensitivityProfiles(pSet=common.cgp.ccle$CGP, sensitivity.measure="auc_published", summary.stat="median")
  cgp.ic50 <- PharmacoGx::summarizeSensitivityProfiles(pSet=common.cgp.ccle$CGP, sensitivity.measure="ic50_published", summary.stat="median")
  ### CCLE
  ccle.ge <- PharmacoGx::summarizeMolecularProfiles(pSet=common.cgp.ccle$CCLE, mDataType="rna", summary.stat="median")
  ccle.ge <- ccle.ge[common.features,]
  ccle.auc <- PharmacoGx::summarizeSensitivityProfiles(pSet=common.cgp.ccle$CCLE, sensitivity.measure="auc_published", summary.stat="median")
  ccle.ic50 <- PharmacoGx::summarizeSensitivityProfiles(pSet=common.cgp.ccle$CCLE, sensitivity.measure="ic50_published", summary.stat="median")
  ### note that ic50 and auc recomputed using a unified pipeline could be selected by using ic50_recomputed or auc_recomputed
  save(list=c("common.features", "cgp.ge", "cgp.auc", "cgp.ic50", "ccle.ge", "ccle.auc", "ccle.ic50"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}


### correlation between cell lines
### spearman
ge.between <- sapply(1:length(cellNames(common.cgp.ccle$CCLE)), function(x){ cor(exprs(ccle.ge)[,x], exprs(cgp.ge)[,x], method="spearman", use="pairwise.complete.obs") })
auc.between <- cor(ccle.auc, cgp.auc, method="spearman", use="pairwise.complete.obs")
ic50.between <- cor(ccle.ic50, cgp.ic50, method="spearman", use="pairwise.complete.obs")
w1 <- wilcox.test(x=ge.between, y=diag(auc.between), conf.int=TRUE)
w2 <- wilcox.test(x=ge.between, y=diag(ic50.between), conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E", w1$p.value, w2$p.value)
pdf(file.path(saveres, "cgp_ccle_spearman_between_cellines_boxplot.pdf"), height=7, width=7)
boxplot(list("GE"=ge.between, "AUC"=diag(auc.between), "IC50"=diag(ic50.between)), main="Concordance between cell lines\nSpearman", ylab=expression(rho), sub=ss, ylim=yylim, col="lightgrey", pch=20, border="black")
dev.off()
### pearson
ge.between <- sapply(1:length(cellNames(common.cgp.ccle$CCLE)), function(x){ cor(exprs(ccle.ge)[,x], exprs(cgp.ge)[,x], method="pearson", use="pairwise.complete.obs") })
auc.between <- cor(ccle.auc, cgp.auc, method="pearson", use="pairwise.complete.obs")
ic50.between <- cor(ccle.ic50, cgp.ic50, method="pearson", use="pairwise.complete.obs")
w1 <- wilcox.test(x=ge.between, y=diag(auc.between), conf.int=TRUE)
w2 <- wilcox.test(x=ge.between, y=diag(ic50.between), conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E", w1$p.value, w2$p.value)
pdf(file.path(saveres, "cgp_ccle_pearson_between_cellines_boxplot.pdf"), height=7, width=7)
boxplot(list("GE"=ge.between, "AUC"=diag(auc.between), "IC50"=diag(ic50.between)), main="Concordance between cell lines\nPearson", ylab=expression(rho), sub=ss, ylim=yylim, col="lightgrey", pch=20, border="black")
dev.off()

### correlation across cell lines
### spearman
ge.across <- sapply(1:length(common.features), function(x){ cor(exprs(ccle.ge)[x,], exprs(cgp.ge)[x,], method="spearman", use="pairwise.complete.obs") })
auc.across <- cor(t(ccle.auc), t(cgp.auc), method="spearman", use="pairwise.complete.obs")
ic50.across <- cor(t(ccle.ic50), t(cgp.ic50), method="spearman", use="pairwise.complete.obs")
w1 <- wilcox.test(x=ge.across, y=diag(auc.across), conf.int=TRUE)
w2 <- wilcox.test(x=ge.across, y=diag(ic50.across), conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E", w1$p.value, w2$p.value)
pdf(file.path(saveres, "cgp_ccle_spearman_across_cellines_boxplot.pdf"), height=7, width=7)
boxplot(list("GE"=ge.across, "AUC"=diag(auc.across), "IC50"=diag(ic50.across)), main="Concordance across cell lines\nSpearman", ylab=expression(rho), sub=ss, ylim=yylim, col="lightgrey", pch=20, border="black")
dev.off()
### pearson
ge.across <- sapply(1:length(common.features), function(x){ cor(exprs(ccle.ge)[x,], exprs(cgp.ge)[x,], method="pearson", use="pairwise.complete.obs") })
auc.across <- cor(t(ccle.auc), t(cgp.auc), method="pearson", use="pairwise.complete.obs")
ic50.across <- cor(t(ccle.ic50), t(cgp.ic50), method="pearson", use="pairwise.complete.obs")
w1 <- wilcox.test(x=ge.across, y=diag(auc.across), conf.int=TRUE)
w2 <- wilcox.test(x=ge.across, y=diag(ic50.across), conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E", w1$p.value, w2$p.value)
pdf(file.path(saveres, "cgp_ccle_pearson_across_cellines_boxplot.pdf"), height=7, width=7)
boxplot(list("GE"=ge.across, "AUC"=diag(auc.across), "IC50"=diag(ic50.across)), main="Concordance across cell lines\nPearson", ylab=expression(rho), sub=ss, ylim=yylim, col="lightgrey", pch=20, border="black")
dev.off()

#################################################
### save session info
#################################################

write(toLatex(sessionInfo(), locale=FALSE), file=file.path(saveres, "sessionInfoR.tex"), append=FALSE)

### end
