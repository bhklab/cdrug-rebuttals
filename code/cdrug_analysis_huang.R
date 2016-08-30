# args <- commandArgs(trailingOnly = TRUE)
# if(length(nbcore) == 0){ nbcore <- as.integer(args[1]) }

## load functions
source(file.path("code", "cdrug_foo.R"))

saveres <- file.path("output")

### install all the libraries at once
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("VennDiagram", "Hmisc", "xtable", "RColorBrewer", "pROC", "Biobase", "genefu", "PharmacoGx"))
library(VennDiagram)
library(Hmisc)
library(xtable)
library(RColorBrewer)
library(pROC)
library(Biobase)
library(genefu)
# install the latest devel version of the PharmacoGx package
# library(devtools)
# devtools::install_github("bhklab/PharmacoGx", ref="master", lib="/mnt/work1/users/bhklab/Rlib/")
library(PharmacoGx)

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
## get pharmacogenomic datasets
#################################################


myfn <- file.path(saveres, "data_cgp_ccle.RData")
if (!file.exists(myfn)) {
  ### download curated pharmacogenomic data from CGP and CCLE
  CGP <- PharmacoGx::downloadPSet("GDSC_2013", saveDir=file.path(saveres, "PSets")) 
  CCLE <- PharmacoGx::downloadPSet("CCLE_2013", saveDir=file.path(saveres, "PSets")) 
  if(confine.analyses.to.nature.common.cell.lines) {
    common <- intersectPSet(pSets = list("CCLE"=CCLE, "CGP"=CGP), intersectOn = c("cell.lines", "drugs"), cells=nature2013.common.cellines)  
  } else { 
    common <- intersectPSet(pSets = list("CCLE"=CCLE, "CGP"=CGP), intersectOn = c("cell.lines", "drugs"))
  }
  ### 512 cell lines and 15 drugs in comon

  ### reorder drugs to put nioltinib first
  drugix <- rownames(common$CGP@drug)
  drugix <- c(which(drugix == "Nilotinib"), setdiff(1:length(drugix), which(drugix == "Nilotinib")))

  ### common affymetrix probes
  common.features <- intersect(rownames(featureInfo(common$CCLE, "rna"))[featureInfo(common$CCLE, "rna")[,"BEST"]==T], rownames(featureInfo(common$CGP, "rna"))[featureInfo(common$CGP, "rna")[,"BEST"]==T])
  ### CGP
  cgp.ge <- summarizeMolecularProfiles(pSet=common$CGP, mDataType="rna", summary.stat="median")
  cgp.ge <- cgp.ge[common.features,]
  cgp.auc <- summarizeSensitivityProfiles(pSet=common$CGP, sensitivity.measure="auc_published", summary.stat="median")[drugix, , drop=FALSE]
  cgp.ic50 <- summarizeSensitivityProfiles(pSet=common$CGP, sensitivity.measure="ic50_published", summary.stat="median")[drugix, , drop=FALSE]
  ### CCLE
  ccle.ge <- summarizeMolecularProfiles(pSet=common$CCLE, mDataType="rna", summary.stat="median")
  ccle.ge <- ccle.ge[common.features,]
  ccle.auc <- summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="auc_published", summary.stat="median")[drugix, , drop=FALSE]
  ccle.ic50 <- summarizeSensitivityProfiles(pSet=common$CCLE, sensitivity.measure="ic50_published", summary.stat="median")[drugix, , drop=FALSE]
  ### note that ic50 and auc recomputed using a unified pipeline could be selected by using ic50_recomputed or auc_recomputed

  save(list=c("CGP", "CCLE", "common", "common.features", "cgp.ge", "cgp.auc", "cgp.ic50", "ccle.ge", "ccle.auc", "ccle.ic50"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}

mycol <- RColorBrewer::brewer.pal(n=7, name="Set1")
### venn diagram of common cell lines between studies
pdf(file.path(saveres, "celllines.pdf"), height=4, width=4)
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@cell), area2=nrow(CGP@cell), cross.area=nrow(common$CCLE@cell), fill=c(mycol[1], mycol[2]), lty="blank",cex=1.5, cat.cex=1, cat.col = c("black", "black"))
dev.off()
### venn diagram of common drugs between studies
pdf(file.path(saveres, "drugs.pdf"), height=4, width=4)
venn.plot <- VennDiagram::draw.pairwise.venn(area1=nrow(CCLE@drug), area2=nrow(CGP@drug), cross.area=nrow(common$CCLE@drug), fill=c(mycol[1], mycol[2]), lty="blank",cex=1.5, cat.cex=1, cat.col = c("black", "black"))
dev.off()

### breaking down common cell lines based on tissue types
tt <- table(common$CCLE@cell[, "tissueid"])
mm <- cbind("Tissue type"= Hmisc::capitalize(gsub("_", " ", names(tt))), "Number of cell lines"=tt)
mm <- mm[mm[ , 2] != 0, , drop=FALSE]
mm <- mm[order(as.numeric(mm[ , 2]), decreasing=TRUE), , drop=FALSE]
xtable::print.xtable(xtable::xtable(mm), include.rownames=FALSE, floating=FALSE, file=file.path(saveres, "tissue_type.tex"), append=FALSE)

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
## Supplementary Figure 2 : box plots, correlation, across vs between
#################################################

### correlation between cell lines
### spearman
ge.between <- sapply(1:length(cellNames(common$CCLE)), function(x){ cor(exprs(ccle.ge)[,x], exprs(cgp.ge)[,x], method="spearman", use="pairwise.complete.obs") })
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
ge.between <- sapply(1:length(cellNames(common$CCLE)), function(x){ cor(exprs(ccle.ge)[,x], exprs(cgp.ge)[,x], method="pearson", use="pairwise.complete.obs") })
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
## Supplementary Figure 3: AMCC for drug AUC across cell lines
#################################################

min.cat <- 3
myfn <- file.path(saveres, "auc_cgp_ccle_amcc_across.RData")
if (!file.exists(myfn)) {
  pdf(file.path(saveres, "auc_cgp_ccle_amcc_across.pdf"), height=5, width=9)
  mcc.auc <- NULL
  for(drugn in drugNames(common$CGP)) {
    # message(sprintf("compute AMCC for %s", drugn))
    xx <- cgp.auc[drugn , ]
    yy <- ccle.auc[drugn , ]
    ccix <- complete.cases(xx, yy)
    mm <- PharmacoGx::amcc(x=xx, y=yy, min.cat=min.cat, nperm=10^3, nthread=nbcore)
    mcc.auc <- rbind(mcc.auc, mm$amcc)
    mm <- mm$mcc
    rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
    mccix <- max(which(mm[-rmix, "estimate"] == max(mm[-rmix, "estimate"], na.rm=TRUE))) + (min.cat - 1)
    par(mfrow=c(1, 2))
    ## scatterplot with sperman correlation
    myScatterPlot(xx, yy, main=drugn, xlab="AUC (CGP)", ylab="AUC (CCLE)")
    rs <- cor.test(x=xx, y=yy, method=cor.method, use="pairwise.complete.obs")
    oo <- sort(xx[ccix], decreasing=TRUE)
    abline(v=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
    oo <- sort(yy[ccix], decreasing=TRUE)
    abline(h=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
    legend("topleft", legend=c(sprintf("rho=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# cell lines=%i", sum(ccix))), col="white", pch=0, bty="n")
    ## mcc plot
    plot(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "estimate"], ylim=c(min(mm[-rmix, "estimate"]), 1), xlab="# sensitive cell lines", ylab="Matthews Correlation coefficient", pch=20, col="lightgrey", main="")
    lines(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "estimate"], col="darkgrey")
    abline(v=(1:nrow(mm[-rmix, ]))[which.max(mm[-rmix, "estimate"])], col="red", lty=2)
    legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "estimate"], mm[mccix, "p.value"]), sprintf("# sensitive cell lines=%i", mccix)), col="white", pch=0, bty="n", text.font=1)
  }
  rownames(mcc.auc) <- drugNames(common$CGP)
  dev.off()
  save(list=c("mcc.auc"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}
drug.amcc <- mcc.auc[, "mcc"]

#################################################
### Supplementary Figure 4
#################################################
## reproducibility between different screening sites
## camptothecin was screened at MGH (drug id 195) and WTSI (drug id 1003)
## data only available in the supplementary infomration of the Nature website

myfn <- file.path(saveres, "nature_supplementary_information.xls")
if (!file.exists(myfn)) {
  dwl.status <- download.file(url="http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip", destfile=file.path(saveres, "nature11005-s2.zip"), quiet=TRUE)
  ff <- as.character(unzip(zipfile=file.path(saveres, "nature11005-s2.zip"), list=TRUE)[1, 1])
  unzip(zipfile=file.path(saveres, "nature11005-s2.zip"), exdir=saveres)
  file.copy(from=file.path(saveres, ff), to=myfn)
}
myfn2 <- file.path(saveres, "nature_supplinfo_drugpheno_cgp.RData")
if(!file.exists(myfn2)) {
  drugpheno.nature <- gdata::read.xls(xls=file.path(saveres, "nature_supplementary_information.xls"), sheet=2)
  drugpheno.nature[drugpheno.nature == "" | drugpheno.nature == " "] <- NA
  save(list="drugpheno.nature", compress=TRUE, file=myfn2)
} else {
  load(myfn2)
}
## format column names
coln2 <- gsub(" ", "", sapply(drugpheno.nature[1,], as.character))
coln2[coln2 == ""] <- NA
drugpheno.nature <- drugpheno.nature[-1, ,drop=FALSE]
coln <- colnames(drugpheno.nature)
coln2[is.na(coln2)] <- coln[is.na(coln2)]
coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
colnames(drugpheno.nature) <- coln2
rownames(drugpheno.nature) <- as.character(drugpheno.nature[ , "Cell.Line"])
myx <- sapply(strsplit(colnames(drugpheno.nature), "_"), function(x) { return(all(x[c(length(x)-1, length(x))] == c("IC", "50"))) })
ic50 <- drugpheno.nature[ , myx, drop=FALSE]
nn <- dimnames(ic50)
nn[[2]] <- gsub("_IC_50", "", nn[[2]])
ic50 <- apply(ic50, 2, as.numeric)
dimnames(ic50) <- nn
ic50 <- exp(ic50) / 10^6
cgp.camptothecin.wtsi <- ic50[ , "drugid_195", drop=FALSE]
cgp.camptothecin.mgh <- ic50[ , "drugid_1003", drop=FALSE]

drugn <- "CAMPTOTHECIN"
xx <- -log10(cgp.camptothecin.wtsi * 10^6)
yy <- -log10(cgp.camptothecin.mgh * 10^6)
ccix <- complete.cases(xx, yy)
mm <- amcc(x=xx, y=yy, min.cat=min.cat, nperm=10^3, nthread=nbcore)
mcc.auc <- mm$amcc
mm <- mm$mcc
rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
mccix <- max(which(mm[-rmix, "estimate"] == max(mm[-rmix, "estimate"], na.rm=TRUE))) + (min.cat - 1)

pdf(file.path(saveres, "auc_cgp_campthotecin_amcc_across.pdf"), height=5, width=9)
par(mfrow=c(1, 2))
## scatterplot with sperman correlation
myScatterPlot(xx, yy, main=drugn, xlab="-Log10 IC50 (WTSI)", ylab="-Log10 IC50 (MGH)")
rs <- cor.test(x=xx, y=yy, method=cor.method, use="pairwise.complete.obs", exact=FALSE)
oo <- sort(xx[ccix], decreasing=TRUE)
abline(v=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
oo <- sort(yy[ccix], decreasing=TRUE)
abline(h=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
legend("topleft", legend=c(sprintf("rho=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# cell lines=%i", sum(ccix))), col="white", pch=0, bty="n")
## mcc plot
plot(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "estimate"], ylim=c(min(mm[-rmix, "estimate"]), 1), xlab="# sensitive cell lines", ylab="Matthews Correlation coefficient", pch=20, col="lightgrey", main="")
lines(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "estimate"], col="darkgrey")
abline(v=(1:nrow(mm[-rmix, ]))[which.max(mm[-rmix, "estimate"])], col="red", lty=2)
legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "estimate"], mm[mccix, "p.value"]), sprintf("# sensitive cell lines=%i", mccix)), col="white", pch=0, bty="n", text.font=1)
dev.off()


drugn <- "AZD6482"
ids <- unlist(strsplit(CGP@drug$drugid[which(drugNames(CGP)==drugn)], split="///"))
tt <- cbind(CGP@sensitivity$info, "drugid2"=sapply(strsplit(rownames(CGP@sensitivity$info), split="_"), function(x){x[2]}))
cells <- intersect(tt[grep(ids[1], rownames(tt)),"cellid"], tt[grep(ids[2], rownames(tt)),"cellid"])
id1 <- rownames(tt)[which(tt$drugid2 == ids[1] & tt$cellid %in% cells)]
id2 <- rownames(tt)[which(tt$drugid2 == ids[2] & tt$cellid %in% cells)]
message(sprintf("compute AMCC for %s", drugn))
xx <- -log10(CGP@sensitivity$profiles[id1, "ic50_published"])
yy <- -log10(CGP@sensitivity$profiles[id2, "ic50_published"])
ccix <- complete.cases(xx, yy)
mm <- amcc(x=xx, y=yy, min.cat=min.cat, nperm=10^3, nthread=nbcore)
mcc.auc <- mm$amcc
mm <- mm$mcc
rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
mccix <- max(which(mm[-rmix, "estimate"] == max(mm[-rmix, "estimate"], na.rm=TRUE))) + (min.cat - 1)

pdf(file.path(saveres, "auc_cgp_AZD6482_amcc_across.pdf"), height=5, width=9)
par(mfrow=c(1, 2))
## scatterplot with sperman correlation
myScatterPlot(xx, yy, main=drugn, xlab="-Log10 IC50 (WTSI)", ylab="-Log10 IC50 (MGH)")
rs <- cor.test(x=xx, y=yy, method=cor.method, use="pairwise.complete.obs", exact=FALSE)
oo <- sort(xx[ccix], decreasing=TRUE)
abline(v=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
oo <- sort(yy[ccix], decreasing=TRUE)
abline(h=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
legend("topleft", legend=c(sprintf("rho=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# cell lines=%i", sum(ccix))), col="white", pch=0, bty="n")
## mcc plot
plot(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "estimate"], ylim=c(min(mm[-rmix, "estimate"]), 1), xlab="# sensitive cell lines", ylab="Matthews Correlation coefficient", pch=20, col="lightgrey", main="")
lines(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "estimate"], col="darkgrey")
abline(v=(1:nrow(mm[-rmix, ]))[which.max(mm[-rmix, "estimate"])], col="red", lty=2)
legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "estimate"], mm[mccix, "p.value"]), sprintf("# sensitive cell lines=%i", mccix)), col="white", pch=0, bty="n", text.font=1)
dev.off()

#################################################
### Supplementary Figure 5
#################################################
## AMCC for gene expression across cell lines
myfn <- file.path(saveres, "ge_cgp_ccle_amcc_across.RData")
if (!file.exists(myfn)) {
  splitix <- parallel::splitIndices(nx=length(common.features), ncl=nbcore)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, xx, yy, min.cat) {
    res <- t(sapply(x, function(x, xx, yy, min.cat) {
      xx <- xx[ , x]
      yy <- yy[ , x]
      res <- amcc(x=xx, y=yy, step.prct=0, min.cat=min.cat, nperm=0, nthread=1)$amcc
      return (res)
    }, xx=xx, yy=yy, min.cat=min.cat))
    return (res)
  }, xx=t(exprs(ccle.ge)), yy=t(exprs(cgp.ge)), min.cat=min.cat)
  mcc.ge <- do.call(rbind, mcres)
  # rownames(mcc.ge) <- rownames(annot.ge)
  dimnames(mcc.ge) <- list(common.features, c("mcc", "p", "n1", "n2", "n"))
  save(list=c("mcc.ge"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}

myfn <- file.path(saveres, "ic50_cgp_ccle_amcc_across.RData")
if (!file.exists(myfn)) {
  pdf(file.path(saveres, "ic50_cgp_ccle_amcc_across.pdf"), height=5, width=9)
  mcc.ic50 <- NULL
  for(drugn in drugNames(common$CGP)) {
    # message(sprintf("compute AMCC for %s", drugn))
    xx <- cgp.ic50[drugn , ]
    yy <- ccle.ic50[drugn , ]
    ccix <- complete.cases(xx, yy)
    mm <- amcc(x=xx, y=yy, min.cat=min.cat, nperm=10^3, nthread=nbcore)
    mcc.ic50 <- rbind(mcc.ic50, mm$amcc)
    mm <- mm$mcc
    rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
    mccix <- max(which(mm[-rmix, "estimate"] == max(mm[-rmix, "estimate"], na.rm=TRUE))) + (min.cat - 1)
    par(mfrow=c(1, 2))
    ## scatterplot with sperman correlation
    myScatterPlot(xx, yy, main=drugn, xlab="AUC (CGP)", ylab="AUC (CCLE)")
    rs <- cor.test(x=xx, y=yy, method=cor.method, use="pairwise.complete.obs")
    oo <- sort(xx[ccix], decreasing=TRUE)
    abline(v=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
    oo <- sort(yy[ccix], decreasing=TRUE)
    abline(h=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
    legend("topleft", legend=c(sprintf("rho=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# cell lines=%i", sum(ccix))), col="white", pch=0, bty="n")
    ## mcc plot
    plot(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "estimate"], ylim=c(min(mm[-rmix, "estimate"]), 1), xlab="# sensitive cell lines", ylab="Matthews Correlation coefficient", pch=20, col="lightgrey", main="")
    lines(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "estimate"], col="darkgrey")
    abline(v=(1:nrow(mm[-rmix, ]))[which.max(mm[-rmix, "estimate"])], col="red", lty=2)
    legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "estimate"], mm[mccix, "p.value"]), sprintf("# sensitive cell lines=%i", mccix)), col="white", pch=0, bty="n", text.font=1)
  }
  rownames(mcc.ic50) <- drugNames(common$CGP)
  dev.off()
  save(list=c("mcc.ic50"), compress=TRUE, file=myfn)
} else {
  load(myfn)
}

## boxplot of AMCC
pdf(file.path(saveres, "boxplot_amcc_across.pdf"), width=6, height=6)
ll <- list("GE"=mcc.ge[ , "mcc"], "AUC"=drug.amcc, "IC50"=mcc.ic50[ , "mcc"])
kt <- kruskal.test(x=ll)
wt1 <- wilcox.test(x=ll$GE, y=ll$AUC)
wt2 <- wilcox.test(x=ll$GE, y=ll$IC50)
boxplot(ll, ylab="AMCC", main="Concordance across cell lines", pch=20, col="lightgrey", sub=sprintf("GE vs. AUC=%.1E\nGE vs. IC50=%.1E", wt1$p.value, wt2$p.value), border="black")
text(x=2.25, y=1, "nilotinib")
dev.off()



#################################################
### Supplementary Figure 6
#################################################
cgp.auc.all <- summarizeSensitivityProfiles(pSet=CGP, sensitivity.measure="auc_published", summary.stat="median")
ccle.auc.all <- summarizeSensitivityProfiles(pSet=CCLE, sensitivity.measure="auc_published", summary.stat="median")

### median abslute deviation of drug sensitivities
### all data
drug.mad.cgp.all <- apply(cgp.auc.all, 1, mad, na.rm=TRUE)
drug.mad.ccle.all <- apply(ccle.auc.all, 1, mad, na.rm=TRUE)
### common data
drug.mad.cgp <- apply(cgp.auc, 1, mad, na.rm=TRUE)
drug.mad.ccle <- apply(ccle.auc, 1, mad, na.rm=TRUE)

## MAD of cytotoxic drugs vs the rest in full studies
iix.cgp <- factor(!is.na(drugInfo(CGP)[ , "Drug.class.II"]) & drugInfo(CGP)[ , "Drug.class.II"] == "Cytotoxic", levels=c(TRUE, FALSE))
levels(iix.cgp) <- c("Cytotoxic", "Targeted")
iix.ccle <- factor(!is.na(drugInfo(CCLE)[ , "Class"]) & drugInfo(CCLE)[ , "Class"] == "Cytotoxic", levels=c(TRUE, FALSE))
levels(iix.ccle) <- c("Cytotoxic", "Targeted")

## find optimal cutoff to descriminate cytotoxic vs targeted drugs based on AUC MAD
rroc.cgp <- pROC::roc(formula=drug.type ~ drug.mad, data=data.frame("drug.mad"=drug.mad.cgp.all, "drug.type"=iix.cgp))
threshold.cgp <- pROC::coords(roc=rroc.cgp, x="best", best.method="youden")[1]
rroc.ccle <- pROC::roc(formula=drug.type ~ drug.mad, data=data.frame("drug.mad"=drug.mad.ccle.all, "drug.type"=iix.ccle))
threshold.ccle <- pROC::coords(roc=rroc.ccle, x="best", best.method="youden")[1]
mad.cytotoxic.cutoff <- round(mean(threshold.cgp, threshold.ccle) * 100) / 100

pdf(file.path(saveres, "cgp_auc_mad_vs_cytotoxic_full.pdf"), width=10, height=5)
par(mfrow=c(1, 2), mar=c(5, 4, 3, 2) + 0.1)
yylim <- c(0, ceiling(max(c(drug.mad.cgp.all, drug.mad.ccle.all), na.rm=TRUE) * 10) / 10)
w1 <- wilcox.test(x=drug.mad.cgp.all[which(iix.cgp == "Cytotoxic")], y=drug.mad.cgp.all[which(iix.cgp == "Targeted")])
ss <- sprintf("Wilcoxon rank sum test = %.1E", w1$p.value)
boxplot(drug.mad.cgp.all ~ iix.cgp, ylim=yylim, col="lightgrey", border="black", pars=list(outcol="black", outpch=20, outcex=0.5), ylab="MAD of AUC", sub=ss, main="Distribution of drug sensitivity(AUC) in CGP")
abline(h=mad.cytotoxic.cutoff, lty=2, col="red")
mtext(side=3, at=par("usr")[1] - 0.33, text=LETTERS[1], line=2, font=2, cex=0.8)
w1 <- wilcox.test(x=drug.mad.ccle.all[which(iix.ccle == "Cytotoxic")], y=drug.mad.ccle.all[which(iix.ccle == "Targeted")])
ss <- sprintf("Wilcoxon rank sum test = %.1E", w1$p.value)
boxplot(drug.mad.ccle.all ~ iix.ccle, ylim=yylim, col="lightgrey", border="black", pars=list(outcol="black", outpch=20, outcex=0.5), ylab="MAD of AUC", sub=ss, main="Distribution of drug sensitivity(AUC) in CCLE")
abline(h=mad.cytotoxic.cutoff, lty=2, col="red")
mtext(side=3, at=par("usr")[1] - 0.33, text=LETTERS[2], line=2, font=2, cex=0.8)
dev.off()

pdf(file.path(saveres, "cgp_ccle_auc_mad_vs_mad.pdf"), width=5, height=5)
par(mar=c(5, 4, 1, 2) + 0.1, cex=0.7)
xxlim <- c(floor(min(drug.mad.cgp, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.cgp, na.rm=TRUE) * 1200) / 1000)
yylim <- c(floor(min(drug.mad.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.ccle, na.rm=TRUE) * 1200) / 1000)
llim <- c(min(xxlim[1], yylim[1]), max(xxlim[2], yylim[2]))
plot(x=drug.mad.cgp, y=drug.mad.ccle, xlim=llim, ylim=llim, pch=20, col=blues9[7], xlab="MAD of AUC in CGP", ylab="MAD of AUC in CCLE")
abline(a=0, b=1, col="darkgrey", lwd=0.5)
abline(h=threshold.ccle, col="red", lty=2, lwd=0.5)
abline(v=threshold.cgp, col="red", lty=2, lwd=0.5)
text(x=drug.mad.cgp, y=drug.mad.ccle, labels=drugNames(common$CGP), cex=0.7, font=1, srt=45, pos=4)
dev.off()


#################################################
### Supplementary Figure 7
#################################################

pdf(file.path(saveres, "cgp_ccle_auc_cor_vs_mad.pdf"), width=10, height=5)
par(mfrow=c(1, 2), mar=c(5, 4, 1, 2) + 0.1, cex=0.7)
xxlim <- c(floor(min(diag(auc.across), na.rm=TRUE) * 1000) / 1000, ceiling(max(diag(auc.across), na.rm=TRUE) * 1200) / 1000)
## variability in cgp
cc <- cor.test(drug.mad.cgp, diag(auc.across), method=cor.method, use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.cgp, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.cgp, na.rm=TRUE) * 1200) / 1000)
plot(x=diag(auc.across), y=drug.mad.cgp, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="Correlation of AUC between CCLE and CGP", ylab="MAD of AUC in CGP")
text(x=diag(auc.across), y=drug.mad.cgp, labels=drugNames(common$CGP), cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("rho=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=1, cex=1)
## variability in ccle
nnn <- sum(complete.cases(drug.mad.ccle, diag(auc.across)))
cc <- cor.test(drug.mad.ccle, diag(auc.across), method=cor.method, use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.ccle, na.rm=TRUE) * 1200) / 1000)
plot(x=diag(auc.across), y=drug.mad.ccle, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="Correlation of AUC between CCLE and CGP", ylab="MAD of AUC in CCLE")
text(x=diag(auc.across), y=drug.mad.ccle, labels=drugNames(common$CGP), cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("rho=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=1, cex=1)
dev.off()

pdf(file.path(saveres, "cgp_ccle_auc_amcc_vs_mad.pdf"), width=10, height=5)
par(mfrow=c(1, 2), mar=c(5, 4, 1, 2) + 0.1, cex=0.7)
xxlim <- c(floor(min(drug.amcc, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.amcc, na.rm=TRUE) * 1200) / 1000)
## variability in cgp
cc <- cor.test(drug.mad.cgp, drug.amcc, method=cor.method, use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.cgp, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.cgp, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.amcc, y=drug.mad.cgp, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="AMCC of AUC between CCLE and CGP", ylab="MAD of AUC in CGP")
text(x=drug.amcc, y=drug.mad.cgp, labels=drugNames(common$CGP), cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("rho=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=1, cex=1)
## variability in ccle
nnn <- sum(complete.cases(drug.mad.ccle, drug.amcc))
cc <- cor.test(drug.mad.ccle, drug.amcc, method=cor.method, use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.ccle, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.amcc, y=drug.mad.ccle, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="AMCC of AUC between CCLE and CGP", ylab="MAD of AUC in CCLE")
text(x=drug.amcc, y=drug.mad.ccle, labels=drugNames(common$CGP), cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("rho=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=1, cex=1)
dev.off()


#################################################
### known biomarkers
#################################################

celline.bcrabl <- c("K-562", "KYO-1", "EM-3", "AR230", "KCL22", "BV-173", "CML-T1", "EM-2", "KU812", "LAMA-84", "MEG-01")
### update CGP PSet
myfn <- file.path(saveres, "gdsc2013_mutation.csv")
if (!file.exists(myfn)) {
  dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub4/cancerrxgene/releases/release-2.0/gdsc_mutation_w2.csv", destfile=myfn, quiet=TRUE)
}
fusion.cgp <- read.csv(myfn, stringsAsFactors=FALSE)
rownames(fusion.cgp) <- rownames(cellInfo(CGP))[match(fusion.cgp[ , "Cell.Line"], cellInfo(CGP)[ , "GDSC.cellid"])]
fusion.cgp[!is.na(fusion.cgp) & (fusion.cgp == "" | fusion.cgp == " " | fusion.cgp == "na")] <- NA
fusion.cgp <- t(fusion.cgp[ , c("BCR_ABL", "EWS_FLI1", "MLL_AFF1")])
mypdata <- data.frame(cbind("batchid"=NA, "cellid"=colnames(fusion.cgp)), stringsAsFactors=FALSE)
rownames(mypdata) <- colnames(fusion.cgp)
myfdata <- data.frame("Symbol"=rownames(fusion.cgp), stringsAsFactors=FALSE)
rownames(myfdata) <- rownames(fusion.cgp)
CGP@molecularProfiles$fusion <- ExpressionSet(assayData=fusion.cgp, phenoData=AnnotatedDataFrame(mypdata), featureData=AnnotatedDataFrame(myfdata))
CGP@molecularProfiles$mutation <- CGP@molecularProfiles$mutation[!is.element(rownames(fData(CGP@molecularProfiles$mutation)), c("BCR_ABL", "EWS_FLI1", "MLL_AFF1")), ]
exprs(CGP@molecularProfiles$fusion)["BCR_ABL", colnames(exprs(CGP@molecularProfiles$fusion)) %in% celline.bcrabl] <- "BCR Exon_13 to ABL Exon_2"
### update CCLE PSet
colnames(fData(CCLE@molecularProfiles$rna))[colnames(fData(CCLE@molecularProfiles$rna)) == "symbol"] <- "Symbol"
fusion.ccle <- CCLE@molecularProfiles$mutation
exprs(fusion.ccle) <- matrix(NA, nrow=1, ncol=ncol(exprs(fusion.ccle)), dimnames=list("BCR_ABL", colnames(exprs(fusion.ccle))))
exprs(fusion.ccle)["BCR_ABL", colnames(exprs(fusion.ccle)) %in% celline.bcrabl] <- "BCR Exon_13 to ABL Exon_2"
myfdata <- data.frame("Symbol"=rownames(fusion.ccle), stringsAsFactors=FALSE)
rownames(myfdata) <- rownames(fusion.ccle)
fData(fusion.ccle) <- myfdata
CCLE@molecularProfiles$fusion <- fusion.ccle

# gene.cgp <- list("rna"=fData(CGP@molecularProfiles$rna)[ , "Symbol"], "mutation"=fData(CGP@molecularProfiles$mutation)[ , "Symbol"], "fusion"=rownames(fData(CGP@molecularProfiles$fusion)))
# known.biomarkers[!is.element(known.biomarkers[ , "gene"], sort(unique(unlist(gene.cgp)))), ]
# gene.ccle <- list("rna"=fData(CCLE@molecularProfiles$rna)[ , "Symbol"], "mutation"=fData(CCLE@molecularProfiles$mutation)[ , "Symbol"], "fusion"=fData(CCLE@molecularProfiles$fusion)[ , "Symbol"])
# known.biomarkers[!is.element(known.biomarkers[ , "gene"], sort(unique(unlist(gene.ccle)))), ]

### for each known biomarker, estimate gene-drug association for mutation, fusion and expression
mut.biomarkers <- drugSensitivitySig(pSet=CGP, mDataType="mutation", drugs=known.biomarkers[, "drug"], features=known.biomarkers[ , "gene"], sensitivity.measure="auc_published", molecular.summary.stat="or")


