########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## April 8, 2014
## March 1, 2015
########################

require(amap) || stop("Library amap is not available!")
require(vcd) || stop("Library vcd is not available!")
require(gplots) || stop("Library gplots is not available!")
require(WriteXLS) || stop("Library WriteXLS is not available!")
require(xtable) || stop("Library xtable is not available!")
require(epibasix) || stop("Library gplots is not available!")
require(OptimalCutpoints) || stop("Library OptimalCutpoints is not available!")

## list of genes
switch (listg,
	"all"={
		## select all the genes with an overall jetset score > 0.20 to remove genes with probes of low coverage, specificty and robustness
		list.genes <- annot.ge[annot.ge[ , "overall"] > 0.20, , drop=FALSE]
	},
	"l1000"={
		## Broad's landmark genes
		list.genes <- l1000.genes[intersect(rownames(annot.ge), rownames(l1000.genes)), , drop=FALSE]
	}
)

## load data shared between CGP and CCLE
load(file.path(saveres, "cdrug2_cgp_ccle_common.RData"))
tissue.ccle <- sampleinfo.ccle[ , "tissue.type"]
names(tissue.ccle) <- rownames(sampleinfo.ccle)
tissue.cgp <- sampleinfo.cgp[ , "tissue.type"]
names(tissue.cgp) <- rownames(sampleinfo.cgp)
drugsn <- gsub("drugid_", "", rownames(druginfo))
drugs.color <- rainbow(length(drugsn), v=0.9)
names(drugs.color) <- drugsn
nature2013.common.cellines <- read.csv(file=file.path("code", "HaibeKains_Nature_2013_common_cellines.csv"))[ , 1]
sampleinfo <- sampleinfo.cgp
tissue <- sampleinfo[ , "tissue.type"]
names(tissue) <- rownames(sampleinfo)
annot.ge <- annot.ge[rownames(list.genes), , drop=FALSE]
data.ge.ccle <- data.ge.ccle[ , rownames(list.genes), drop=FALSE]
data.ge.cgp <- data.ge.cgp[ , rownames(list.genes), drop=FALSE]

########################
## sample size per tissue types
########################

ccelln <- intersect(rownames(sampleinfo.cgp), rownames(sampleinfo.ccle))
tt <- table(sampleinfo.cgp[ccelln, "tissue.type"], sampleinfo.cgp[ccelln, "tissue.type"])
if(any(tt[upper.tri(tt)] > 0) || any(tt[lower.tri(tt)] > 0)) { warning("Discrepancies in tissue type classification between CGP and CCLE") }
mm <- cbind("Tissue type"=names(diag(tt)), "Number of cell lines"=diag(tt))
mm <- mm[mm[ , 2] != 0, , drop=FALSE]
mm <- mm[order(as.numeric(mm[ , 2]), decreasing=TRUE), , drop=FALSE]
xtable::print.xtable(xtable::xtable(mm), include.rownames=FALSE, floating=FALSE, file=file.path(saveres, "tissue_type_cgp_ccle_paper.tex"), append=FALSE)

utissue <- table(as.character(tissue))
coltissue <- coltissuet <- factor(names(utissue))
# levels(coltissue) <- gplots::rich.colors(length(utissue))
levels(coltissue) <- rainbow(n=length(utissue), s=0.5, v=1)
levels(coltissuet) <- rainbow(n=length(utissue), s=0.75, v=1)
coltissue <- as.character(coltissue)
coltissuet <- as.character(coltissuet)
names(coltissue) <- names(coltissuet) <- names(utissue)

########################
## concordance based on AMCC
########################

## AMCC for drug IC50 across cell lines
myfn <- file.path(saveres, "ic50_cgp_ccle_amcc_across.RData")
if (!file.exists(myfn)) {
  pdf(file.path(saveres, "ic50_cgp_ccle_amcc_across.pdf"), height=5, width=9)
  mcc.ic50 <- NULL
  for(i in 1:nrow(druginfo)) {
    drugn <- gsub("drugid_", "", rownames(druginfo))[i]
    # message(sprintf("compute AMCC for %s", drugn))
    xx <- drugpheno.cgp$IC50[ , i]
    yy <- drugpheno.ccle$IC50[ , i]
	 ccix <- complete.cases(xx, yy)
    mm <- amcc(x=xx, y=yy, min.cat=min.cat, nperm=10^3, nthread=nbcore)
	 mcc.ic50 <- rbind(mcc.ic50, mm$amcc)
	 mm <- mm$mcc
	 rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
	 mccix <- max(which(mm[-rmix, "mcc"] == max(mm[-rmix, "mcc"], na.rm=TRUE))) + (min.cat - 1)
	 par(mfrow=c(1, 2))
    ## scatterplot with sperman correlation
    myScatterPlot(xx, yy, main=drugn, xlab="IC50 (CGP)", ylab="IC50 (CCLE)")
    rs <- cor.test(x=xx, y=yy, method=concordance.method, use="pairwise.complete.obs")
	 oo <- sort(xx[ccix], decreasing=TRUE)
	 abline(v=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
	 oo <- sort(yy[ccix], decreasing=TRUE)
	 abline(h=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
    legend("topleft", legend=c(sprintf("Rs=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# cell lines=%i", sum(ccix))), col="white", pch=0, bty="n")
    ## mcc plot
    plot(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], ylim=c(min(mm[-rmix, "mcc"]), 1), xlab="# sensitive cell lines", ylab="Matthews Correlation coefficient", pch=20, col="lightgrey", main="")
    lines(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], col="darkgrey")
	 abline(v=(1:nrow(mm[-rmix, ]))[which.max(mm[-rmix, "mcc"])], col="red", lty=2)
    legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "mcc"], mm[mccix, "p"]), sprintf("# sensitive cell lines=%i", mccix)), col="white", pch=0, bty="n", text.font=1)
  }
  rownames(mcc.ic50) <- rownames(druginfo)
  dev.off()
  save(list=c("mcc.ic50"), compress=TRUE, file=myfn)
} else { load(myfn) }
mcc.ic50.across <- mcc.ic50

## AMCC for drug AUC across cell lines
myfn <- file.path(saveres, "auc_cgp_ccle_amcc_across.RData")
if (!file.exists(myfn)) {
  pdf(file.path(saveres, "auc_cgp_ccle_amcc_across.pdf"), height=5, width=9)
  mcc.auc <- NULL
  for(i in 1:nrow(druginfo)) {
    drugn <- gsub("drugid_", "", rownames(druginfo))[i]
    # message(sprintf("compute AMCC for %s", drugn))
    xx <- drugpheno.cgp$AUC[ , i]
    yy <- drugpheno.ccle$AUC[ , i]
	 ccix <- complete.cases(xx, yy)
    mm <- amcc(x=xx, y=yy, min.cat=min.cat, nperm=10^3, nthread=nbcore)
	 mcc.auc <- rbind(mcc.auc, mm$amcc)
	 mm <- mm$mcc
	 rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
	 mccix <- max(which(mm[-rmix, "mcc"] == max(mm[-rmix, "mcc"], na.rm=TRUE))) + (min.cat - 1)
	 par(mfrow=c(1, 2))
    ## scatterplot with sperman correlation
    myScatterPlot(xx, yy, main=drugn, xlab="AUC (CGP)", ylab="AUC (CCLE)")
    rs <- cor.test(x=xx, y=yy, method=concordance.method, use="pairwise.complete.obs")
	 oo <- sort(xx[ccix], decreasing=TRUE)
	 abline(v=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
	 oo <- sort(yy[ccix], decreasing=TRUE)
	 abline(h=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
    legend("topleft", legend=c(sprintf("Rs=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# cell lines=%i", sum(ccix))), col="white", pch=0, bty="n")
    ## mcc plot
    plot(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], ylim=c(min(mm[-rmix, "mcc"]), 1), xlab="# sensitive cell lines", ylab="Matthews Correlation coefficient", pch=20, col="lightgrey", main="")
    lines(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], col="darkgrey")
	 abline(v=(1:nrow(mm[-rmix, ]))[which.max(mm[-rmix, "mcc"])], col="red", lty=2)
    legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "mcc"], mm[mccix, "p"]), sprintf("# sensitive cell lines=%i", mccix)), col="white", pch=0, bty="n", text.font=1)
  }
  rownames(mcc.auc) <- rownames(druginfo)
  dev.off()
  save(list=c("mcc.auc"), compress=TRUE, file=myfn)
} else { load(myfn) }
mcc.auc.across <- mcc.auc

## AMCC for drug AUC across cell lines
myfn <- file.path(saveres, "auc_cgp_campthotecin_amcc_across.RData")
load(file.path(saveres, "cgp_camptothecin_mgh_wtsi.RData"))
if (!file.exists(myfn)) {
	pdf(file.path(saveres, "auc_cgp_campthotecin_amcc_across.pdf"), height=5, width=9)

	drugn <- "CAMPTOTHECIN"
	# message(sprintf("compute AMCC for %s", drugn))
	xx <- -log10(cgp.camptothecin.wtsi)
	yy <- -log10(cgp.camptothecin.mgh)
	ccix <- complete.cases(xx, yy)
	mm <- amcc(x=xx, y=yy, min.cat=min.cat, nperm=10^3, nthread=nbcore)
	mcc.auc <- mm$amcc
	mm <- mm$mcc
	rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
	mccix <- max(which(mm[-rmix, "mcc"] == max(mm[-rmix, "mcc"], na.rm=TRUE))) + (min.cat - 1)
	par(mfrow=c(1, 2))
	## scatterplot with sperman correlation
	myScatterPlot(xx, yy, main=drugn, xlab="-Log10 IC50 (WTSI)", ylab="-Log10 IC50 (MGH)")
	rs <- cor.test(x=xx, y=yy, method=concordance.method, use="pairwise.complete.obs")
	oo <- sort(xx[ccix], decreasing=TRUE)
	abline(v=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
	oo <- sort(yy[ccix], decreasing=TRUE)
	abline(h=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
	legend("topleft", legend=c(sprintf("Rs=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# cell lines=%i", sum(ccix))), col="white", pch=0, bty="n")
	## mcc plot
	plot(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], ylim=c(min(mm[-rmix, "mcc"]), 1), xlab="# sensitive cell lines", ylab="Matthews Correlation coefficient", pch=20, col="lightgrey", main="")
	lines(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], col="darkgrey")
	abline(v=(1:nrow(mm[-rmix, ]))[which.max(mm[-rmix, "mcc"])], col="red", lty=2)
	legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "mcc"], mm[mccix, "p"]), sprintf("# sensitive cell lines=%i", mccix)), col="white", pch=0, bty="n", text.font=1)

  dev.off()
  save(list=c("mcc.auc"), compress=TRUE, file=myfn)
} else { load(myfn) }

## AMCC for gene expression across cell lines
myfn <- file.path(saveres, sprintf("ge_cgp_ccle_amcc_across_%s.RData", listg))
if (!file.exists(myfn)) {
  splitix <- parallel::splitIndices(nx=ncol(data.ge.cgp), ncl=nbcore)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, xx, yy, min.cat) {
    res <- t(sapply(x, function(x, xx, yy, min.cat) {
      xx <- xx[ , x]
      yy <- yy[ , x]
      res <- amcc(x=xx, y=yy, step.prct=0, min.cat=min.cat, nperm=0, nthread=1)$amcc
      return (res)
    }, xx=xx, yy=yy, min.cat=min.cat))
    return (res)
  }, xx=data.ge.cgp, yy=data.ge.ccle, min.cat=min.cat)
  mcc.ge <- do.call(rbind, mcres)
  # rownames(mcc.ge) <- rownames(annot.ge)
  dimnames(mcc.ge) <- list(colnames(data.ge.cgp), c("mcc", "p", "n1", "n2", "n"))
  save(list=c("mcc.ge"), compress=TRUE, file=myfn)
} else { load(myfn) }
mcc.ge.across <- mcc.ge

## AMCC for drug IC50 between cell lines
myfn <- file.path(saveres, "ic50_cgp_ccle_amcc_between.RData")
if (!file.exists(myfn)) {
	# pdf(file.path(saveres, "ic50_cgp_ccle_amcc_between.pdf"), height=5, width=9)
	mcc.ic50 <- NULL
	for(i in 1:nrow(sampleinfo)) {
	celln <- rownames(sampleinfo)[i]
		# message(sprintf("compute AMCC for %s", celln))
		xx <- -log10(drugpheno.cgp$IC50[i, ] / 10^6)
		yy <- -log10(drugpheno.ccle$IC50[i, ] / 10^6)
		mcc.ic50 <- rbind(mcc.ic50, amcc(x=xx, y=yy, min.cat=min.cat, nperm=10^3, nthread=nbcore)$amcc)
  }
  rownames(mcc.ic50) <- rownames(sampleinfo)
  save(list=c("mcc.ic50"), compress=TRUE, file=myfn)
} else { load(myfn) }
mcc.ic50.between <- mcc.ic50

## AMCC for drug AUC between cell lines
myfn <- file.path(saveres, "auc_cgp_ccle_amcc_between.RData")
if (!file.exists(myfn)) {
  # pdf(file.path(saveres, "auc_cgp_ccle_amcc_between.pdf"), height=5, width=9)
  mcc.auc <- NULL
  for(i in 1:nrow(sampleinfo)) {
    celln <- rownames(sampleinfo)[i]
    # message(sprintf("compute AMCC for %s", celln))
    xx <- drugpheno.cgp$AUC[i, ]
    yy <- drugpheno.ccle$AUC[i, ]
	 mcc.auc <- rbind(mcc.auc, amcc(x=xx, y=yy, min.cat=min.cat, nperm=10^3, nthread=nbcore)$amcc)
  }
  rownames(mcc.auc) <- rownames(sampleinfo)
  save(list=c("mcc.auc"), compress=TRUE, file=myfn)
} else { load(myfn) }
mcc.auc.between <- mcc.auc

## AMCC for gene expression between cell lines
myfn <- file.path(saveres, sprintf("ge_cgp_ccle_amcc_between_%s.RData", listg))
if (!file.exists(myfn)) {
   splitix <- parallel::splitIndices(nx=nrow(sampleinfo), ncl=nbcore)
   splitix <- splitix[sapply(splitix, length) > 0]
   mcres <- parallel::mclapply(splitix, function(x, ss, xx, yy, list.genes, min.cat) {
      mcc.ge <- NULL
		for (i in 1:length(x)) {
			celln <- rownames(ss)[x[i]]
	      # message(sprintf("compute AMCC for %s", celln))
	      xx2 <- xx[celln, rownames(list.genes)]
	      yy2 <- yy[celln, rownames(list.genes)]
			res <- amcc(x=xx2, y=yy2, step.prct=0, min.cat=min.cat, nperm=0, nthread=1)$amcc
			mcc.ge <- rbind(mcc.ge, res)
		}
		return (mcc.ge)
	}, ss=sampleinfo, xx=data.ge.cgp, yy=data.ge.ccle, list.genes=list.genes, min.cat=min.cat)
	mcc.ge <- do.call(rbind, mcres)
	dimnames(mcc.ge) <- list(rownames(sampleinfo), c("mcc", "p", "n1", "n2", "n"))
	save(list=c("mcc.ge"), compress=TRUE, file=myfn)
} else { load(myfn) }
mcc.ge.between <- mcc.ge


## boxplot of AMCC
pdf(file.path(saveres, "boxplot_amcc_across.pdf"), width=6, height=6)
ll <- list("GE"=mcc.ge.across[ , "mcc"], "AUC"=mcc.auc.across[ , "mcc"], "IC50"=mcc.ic50.across[ , "mcc"])
kt <- kruskal.test(x=ll)
wt1 <- wilcox.test(x=ll$GE, y=ll$AUC)
wt2 <- wilcox.test(x=ll$GE, y=ll$IC50)
boxplot(ll, ylab="AMCC", main="Concordance across cell lines", pch=20, col="lightgrey", sub=sprintf("GE vs. AUC=%.1E\nGE vs. IC50=%.1E", wt1$p.value, wt2$p.value), border="black")
text(x=2.25, y=1, "nilotinib")
dev.off()

pdf(file.path(saveres, "boxplot_amcc_between.pdf"), width=6, height=6)
ll <- list("GE"=mcc.ge.between[ , "mcc"], "AUC"=mcc.auc.between[ , "mcc"], "IC50"=mcc.ic50.between[ , "mcc"])
kt <- kruskal.test(x=ll)
wt1 <- wilcox.test(x=ll$GE, y=ll$AUC)
wt2 <- wilcox.test(x=ll$GE, y=ll$IC50)
mp <- boxplot(ll, ylab="AMCC", main="Concordance between cell lines", pch=20, col="lightgrey", sub=sprintf("GE vs. AUC=%.1E\nGE vs. IC50=%.1E", wt1$p.value, wt2$p.value), border="black")
dev.off()

## create tables summarizing all the correlations
tt <- matrix(NA, nrow=length(drugsn), ncol=2, dimnames=list(drugsn, c("drug.sensitivity", "gene.drug")))
correlations <- list("ic50"=tt, "ic50.call"=tt, "auc"=tt, "auc.call"=tt)
## correlation statistics
tt <- matrix(NA, nrow=length(drugsn), ncol=5, dimnames=list(drugsn, c("rho", "lower", "upper", "p", "n")))
correlations.stats <- list("ic50"=tt, "ic50.call"=tt, "auc"=tt, "auc.call"=tt)

## consistency between IC50 with spearman correlation
pdf(file.path(saveres, "cgp_ccle_scatterplot_ic50_pres.pdf"), height=9, width=16)
par(mfrow=c(3, 5), cex=0.8, las=1)
for(i in 1:nrow(druginfo)) {
  drugn <- gsub("drugid_", "", colnames(drugpheno.ccle$IC50)[i])
  xx <- -log10(drugpheno.cgp$IC50[, i] / 10^6)
  yy <- -log10(drugpheno.ccle$IC50[, i] / 10^6)
  xxlim <- c(floor(min(xx, na.rm=TRUE) * 10) / 10, ceiling(max(xx, na.rm=TRUE) * 10) / 10)
  yylim <- c(floor(min(yy, na.rm=TRUE) * 10) / 10, ceiling(max(yy, na.rm=TRUE) * 10) / 10)
  yylim[2] <- yylim[2] + (yylim[2] * 0.1)
  nnn <- sum(complete.cases(xx, yy))
  if(nnn >= minsample) {
    cc <- cor.test(xx, yy, method=concordance.method, use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  par(mar=c(4, 4, 3, 1) + 0.1)
  mycol <- rep(blues9[7], length(xx))
  names(mycol) <- names(xx)
  mycol[!is.element(names(mycol), nature2013.common.cellines)] <- "red3"
  myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 10, "-log10 IC50 (CGP)", ""), ylab=ifelse((i %% 5) == 1, "-log10 IC50 (CCLE)", ""), main=drugn, xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75, col=mycol)
  abline(a=0, b=1, col="black")
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2, cex=1)
  correlations[["ic50"]][i, "drug.sensitivity"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  correlations.stats[["ic50"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
}
dev.off()
## scatter plot
pdf(file.path(saveres, "cgp_ccle_scatterplot_ic50_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:nrow(druginfo)) {
  drugn <- gsub("drugid_", "", colnames(drugpheno.ccle$IC50)[i])
  xx <- -log10(drugpheno.cgp$IC50[, i] / 10^6)
  yy <- -log10(drugpheno.ccle$IC50[, i] / 10^6)
  xxlim <- c(floor(min(xx, na.rm=TRUE) * 10) / 10, ceiling(max(xx, na.rm=TRUE) * 10) / 10)
  yylim <- c(floor(min(yy, na.rm=TRUE) * 10) / 10, ceiling(max(yy, na.rm=TRUE) * 10) / 10)
  nnn <- sum(complete.cases(xx, yy))
  par(mar=c(4, 4, 3, 1) + 0.1)
  mycol <- rep(blues9[7], length(xx))
  names(mycol) <- names(xx)
  mycol[!is.element(names(mycol), nature2013.common.cellines)] <- "red3"
  myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 11, "-log10 IC50 (CGP)", ""), ylab=ifelse((i %% 4) == 1, "-log10 IC50 (CCLE)", ""), main=drugn, xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75, col=mycol)
  abline(a=0, b=1, col="black")
}
dev.off()

## consistency between AUC with spearman correlation
pdf(file.path(saveres, "cgp_ccle_scatterplot_auc_pres.pdf"), height=9, width=16)
par(mfrow=c(3, 5), cex=0.8, las=1)
for(i in 1:nrow(druginfo)) {
  drugn <- gsub("drugid_", "", colnames(drugpheno.ccle$IC50)[i])
  xx <- drugpheno.cgp$AUC[, i]
  yy <- drugpheno.ccle$AUC[, i]
  xxlim <- c(0, ceiling(max(xx, na.rm=TRUE) * 10) / 10)
  yylim <- c(0, ceiling(max(yy, na.rm=TRUE) * 10) / 10)
  yylim[2] <- yylim[2] + 0.1
  nnn <- sum(complete.cases(xx, yy))
  if(nnn >= minsample) {
    cc <- cor.test(xx, yy, method=concordance.method, use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  par(mar=c(4, 4, 3, 1) + 0.1)
  mycol <- rep(blues9[7], length(xx))
  names(mycol) <- names(xx)
  mycol[!is.element(names(mycol), nature2013.common.cellines)] <- "red3"
  myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 10, "AUC (CGP)", ""), ylab=ifelse((i %% 5) == 1, "AUC (CCLE)", ""), main=drugn, xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75, col=mycol)
  abline(a=0, b=1, col="black")
  # legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2, cex=1)
  correlations[["auc"]][i, "drug.sensitivity"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  correlations.stats[["auc"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
}
dev.off()

## scatter plot with barplot
pdf(file.path(saveres, "cgp_ccle_scatterplot_auc_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), cex=0.8, las=1)
for(i in 1:nrow(druginfo)) {
  drugn <- gsub("drugid_", "", colnames(drugpheno.ccle$IC50)[i])
  xx <- drugpheno.cgp$AUC[, i]
  yy <- drugpheno.ccle$AUC[, i]
  xxlim <- c(0, ceiling(max(xx, na.rm=TRUE) * 10) / 10)
  yylim <- c(0, ceiling(max(yy, na.rm=TRUE) * 10) / 10)
  nnn <- sum(complete.cases(xx, yy))
  par(mar=c(4, 4, 3, 1) + 0.1)
  mycol <- rep(blues9[7], length(xx))
  names(mycol) <- names(xx)
  mycol[!is.element(names(mycol), nature2013.common.cellines)] <- "red3"
  myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 11, "AUC (CGP)", ""), ylab=ifelse((i %% 4) == 1, "AUC (CCLE)", ""), main=drugn, xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75, col=mycol)
  abline(a=0, b=1, col="black")
}
dev.off()

########################
## correlation across cell lines
########################

myfn <- file.path(saveres, sprintf("cgp_ccle_spearman_across_cellines_%s.RData", listg))
if (!file.exists(myfn)) {
  ## correlation for gene expression across cell lines
  # gg <- rownames(annot.ge)
  gg <- rownames(list.genes)
  ge.cor <- sapply(gg, function (x, d1, d2) {
    return (cor(d1[ , x], d2[ , x], method=concordance.method, use="pairwise.complete.obs"))
  }, d1=data.ge.cgp[ , rownames(list.genes), drop=FALSE], d2=data.ge.ccle[ , rownames(list.genes), drop=FALSE])
  ## correlation for ic50 across cell lines
  dd <- rownames(druginfo)
  ic50.cor <- sapply(dd, function (x, d1, d2) {
    return (cor(d1[ , x], d2[ , x], method=concordance.method, use="pairwise.complete.obs"))
  }, d1=drugpheno.cgp$IC50, d2=drugpheno.ccle$IC50)
  ## correlation for auc across cell lines
  dd <- rownames(druginfo)
  auc.cor <- sapply(dd, function (x, d1, d2) {
    return (cor(d1[ , x], d2[ , x], method=concordance.method, use="pairwise.complete.obs"))
  }, d1=drugpheno.cgp$AUC, d2=drugpheno.ccle$AUC)
  ## correlation for mutation across cell lines
  gg <- colnames(mutation.cgp)
  mut.kappa <- sapply(gg, function (x, d1, d2) {
    tt <- table(d1[ , x] != "wt", d2[ , x] != "wt")
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
      rr <- NA
    } else {rr <- rr$kappa }
    return (rr)
  }, d1=mutation.cgp, d2=mutation.ccle)
  ## correlation for IC50 ternary calls
  dd <- rownames(druginfo)
  ic50.call.kappa <- sapply(dd, function (x, d1, d2) {
    tt <- table(unlist(d1[ , x]), unlist(d2[ , x]))
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
      rr <- NA
    } else {
      rr <- rr$kappa
    }
    return (rr)
  }, d1=drugpheno.cgp$IC50.CALL3, d2=drugpheno.ccle$IC50.CALL3)
  ## correlation for AUC ternary calls
  dd <- rownames(druginfo)
  auc.call.kappa <- sapply(dd, function (x, d1, d2) {
    tt <- table(unlist(d1[ , x]), unlist(d2[ , x]))
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
      rr <- NA
    } else {
      rr <- rr$kappa
    }
    return (rr)
  }, d1=drugpheno.cgp$AUC.CALL3, d2=drugpheno.ccle$AUC.CALL3)
  save(list=c("ge.cor", "ic50.cor", "auc.cor", "mut.kappa", "ic50.call.kappa", "auc.call.kappa"), compress=TRUE, file=myfn)
} else { load(myfn) }

pdf(file.path(saveres, sprintf("cgp_ccle_spearman_across_cellines_boxplot_%s.pdf", listg)))
## test significance of the difference between genomic and drug sensitivity data
w1 <- wilcox.test(x=ge.cor, y=auc.cor, conf.int=TRUE)
w2 <- wilcox.test(x=ge.cor, y=ic50.cor, conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E", w1$p.value, w2$p.value)
boxplot(list("GE"=ge.cor, "AUC"=auc.cor, "IC50"=ic50.cor), main="Concordance across cell lines", ylab=expression(R[s]), sub=ss, ylim=yylim, col="lightgrey", pch=20, border="black")
dev.off()

pdf(file.path(saveres, sprintf("cgp_ccle_kappa_across_cellines_boxplot_%s.pdf", listg)))
## test significance of the difference between genomic and drug sensitivity data
w1 <- wilcox.test(x=mut.kappa, y=ic50.call.kappa, conf.int=TRUE)
w2 <- wilcox.test(x=mut.kappa, y=auc.call.kappa, conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("Mutation vs. IC50 calls = %.1E\nMutation vs. AUC calls = %.1E", w1$p.value, w2$p.value)
boxplot(list("Mutations"=mut.kappa, "IC50 calls"=ic50.call.kappa, "AUC calls"=auc.call.kappa), main="Concordance across cell lines", ylab=expression(kappa), sub=ss, ylim=yylim, col="lightgrey", pch=20, border="black")
dev.off()

########################
## correlation between cell lines
########################

myfn <- file.path(saveres, sprintf("cgp_ccle_spearman_between_cellines_%s.RData", listg))
if (!file.exists(myfn)) {
  cellid <- rownames(data.ge.cgp)
  ## correlation for gene expression across cell lines
  ge.cor <- sapply(cellid, function (x, d1, d2) {
    return (cor(d1[x, ], d2[x, ], method=concordance.method, use="pairwise.complete.obs"))
  }, d1=data.ge.cgp[ , rownames(list.genes), drop=FALSE], d2=data.ge.ccle[ , rownames(list.genes), drop=FALSE])
  ## correlation for ic50 across cell lines
  ic50.cor <- sapply(cellid, function (x, d1, d2) {
    return (cor(d1[x, ], d2[x, ], method=concordance.method, use="pairwise.complete.obs"))
  }, d1=drugpheno.cgp$IC50, d2=drugpheno.ccle$IC50)
  ## correlation for auc across cell lines
  auc.cor <- sapply(cellid, function (x, d1, d2) {
    return (cor(d1[x, ], d2[x, ], method=concordance.method, use="pairwise.complete.obs"))
  }, d1=drugpheno.cgp$AUC, d2=drugpheno.ccle$AUC)
  ## correlation for mutation across cell lines
  mut.kappa <- sapply(cellid, function (x, d1, d2) {
    tt <- table(d1[x, ] != "wt", d2[x, ] != "wt")
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
      rr <- NA
    } else {
      rr <- rr$kappa
    }
    return (rr)
  }, d1=mutation.cgp, d2=mutation.ccle)
  ## correlation for IC50 ternary calls
  ic50.call.kappa <- sapply(cellid, function (x, d1, d2) {
    tt <- table(unlist(d1[x, ]), unlist(d2[x, ]))
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
      rr <- NA
    } else {
      rr <- rr$kappa
    }
    return (rr)
  }, d1=drugpheno.cgp$IC50.CALL3, d2=drugpheno.ccle$IC50.CALL3)
  ## correlation for AUC ternary calls
  auc.call.kappa <- sapply(cellid, function (x, d1, d2) {
    tt <- table(unlist(d1[x, ]), unlist(d2[x, ]))
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
      rr <- NA
    } else {
      rr <- rr$kappa
    }
    return (rr)
  }, d1=drugpheno.cgp$AUC.CALL3, d2=drugpheno.ccle$AUC.CALL3)
  save(list=c("ge.cor", "ic50.cor", "auc.cor", "mut.kappa", "ic50.call.kappa", "auc.call.kappa"), compress=TRUE, file=myfn)
} else { load(myfn) }

pdf(file.path(saveres, sprintf("cgp_ccle_spearman_between_cellines_boxplot_%s.pdf", listg)))
## test significance of the difference between genomic and drug sensitivity data
w1 <- wilcox.test(x=ge.cor, y=auc.cor, conf.int=TRUE)
w2 <- wilcox.test(x=ge.cor, y=ic50.cor, conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E", w1$p.value, w2$p.value)
boxplot(list("GE"=ge.cor, "AUC"=auc.cor, "IC50"=ic50.cor), main="Concordance between cell lines", ylab=expression(R[s]), sub=ss, ylim=yylim, col="lightgrey", pch=20, border="black")
dev.off()

pdf(file.path(saveres, sprintf("cgp_ccle_kappa_between_cellines_boxplot_%s.pdf", listg)))
## test significance of the difference between genomic and drug sensitivity data
w1 <- wilcox.test(x=mut.kappa, y=ic50.call.kappa, conf.int=TRUE)
w2 <- wilcox.test(x=mut.kappa, y=auc.call.kappa, conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("Mutation vs. IC50 calls = %.1E\nMutation vs. AUC calls = %.1E", w1$p.value, w2$p.value)
boxplot(list("Mutations"=mut.kappa, "IC50 calls"=ic50.call.kappa, "AUC calls"=auc.call.kappa), main="Concordance between cell lines", ylab=expression(kappa), sub=ss, ylim=yylim, col="lightgrey", pch=20, border="black")
dev.off()


########################
message("Correlation of tuned AUC sensitivity calls:")

## new sensitivity call
auc.call3.cgp <- drugpheno.cgp$AUC.CALL3[rownames(data.ge.cgp), rownames(druginfo), drop=FALSE]
auc.call3.ccle <- drugpheno.ccle$AUC.CALL3[rownames(data.ge.ccle), rownames(druginfo), drop=FALSE]

pdf(file.path(saveres, "cgp_ccle_auc_call3_paper.pdf"), width=8, height=9)
par(mfrow=c(5, 3), mar=c(5, 4, 4, 2) + 0.1)
mykap <- mykap2 <- NULL
for (i in 1:ncol(auc.call3.ccle)) {
  ff1 <- factor(x=auc.call3.ccle[ , i], levels=c("resistant", "intermediate", "sensitive"))
  ff2 <- factor(x=auc.call3.cgp[ , i], levels=c("resistant", "intermediate", "sensitive"))
  ## with intermediate
  levels(ff1) <- c("res", "inter", "sens")
  levels(ff2) <- c("res", "inter", "sens")
  nnn <- sum(complete.cases(ff1, ff2))
  tt <- table("CCLE"=ff1, "CGP"=ff2)
  tts <- vcd::assocstats(tt)
  ## check if all elements are on the diagonal
  ttt <- tt
  diag(ttt) <- 0
  if(sum(ttt) == 0) {
    rr <- rr.l <- rr.u <- 1
  } else {
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
		 rr <- rr.l <- rr.u <- NA
	 } else {
		rr.l <- rr$CIL
		rr.u <- rr$CIU
		rr <- rr$kappa
	}
  }
  mykap <- c(mykap, rr)
  ## without intermediate
  levels(ff1)[2] <- NA
  levels(ff2)[2] <- NA
  nnn2 <- sum(complete.cases(ff1, ff2))
  tt2 <- table("CCLE"=ff1, "CGP"=ff2)
  tts2 <- vcd::assocstats(tt2)
  ## check if all elements are on the diagonal
  ttt <- tt
  diag(ttt) <- 0
  if(sum(ttt) == 0) {
    rr2 <- rr2.l <- rr2.u <- 1
  } else {
    err <- try(rr2 <- epibasix::epiKappa(tt2, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
		 rr2 <- rr2.l <- rr2.u <- NA
	 } else {
		rr2.l <- rr2$CIL
		rr2.u <- rr2$CIU
		rr2 <- rr2$kappa
	}
  }
  mykap2 <- c(mykap2, rr2)
  plot.new()
  plotrix::addtable2plot(x=par("usr")[1] + 0.1, y=par("usr")[4] - 0.4, xjust=0, yjust=0, table=as.matrix(tt), bty="o", display.rownames=TRUE, display.colnames=TRUE, hlines=FALSE, vlines=FALSE, title="CCLE   vs   CGP", cex=1.25)
  title(main=sprintf("AUC sensitivity calling\n%s", gsub("drugid_", "", colnames(auc.call3.ccle)[i])), sub=sprintf("Kappa=%.2g, 95%%CI [%.2g,%.2g], p=%.1E\nKappa2=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", rr, rr.l, rr.u, tts$chisq_tests[1,3], rr2, rr2.l, rr2.u, tts2$chisq_tests[1,3]))
  correlations[["auc.call"]][i, "drug.sensitivity"] <- rr
  ## statistics
  correlations.stats[["auc.call"]][i, ] <- c(rr, rr.l, rr.u, tts$chisq_tests[1,3], nnn)
}
names(mykap) <- names(mykap2) <- colnames(auc.call3.ccle)
dev.off()

########################
message("Correlation of tuned AUC sensitivity calls (15 most sensitive vs 55 most resistant for common cell lines):")

## new sensitivity call
auc.call3.cgp <- drugpheno.cgp$AUC[rownames(data.ge.cgp), rownames(druginfo), drop=FALSE]
auc.call3.ccle <- drugpheno.ccle$AUC[rownames(data.ge.ccle), rownames(druginfo), drop=FALSE]
## select the cell lines that have been tested in both studies
myx <- !is.na(auc.call3.cgp) & !is.na(auc.call3.ccle)
auc.call3.cgp[!myx] <- NA
auc.call3.ccle[!myx] <- NA

myfn <- function(x, n1=30, n2=10) {
	ss <- sort(x, na.last=NA)
	rr <- c(ss[n1], ss[length(ss) - n2 + 1])
	ss <- cut2(x=x, cuts=rr)
	if(length(levels(ss)) == 3) { levels(ss) <- c("resistant", "intermediate", "sensitive") }
	if(length(levels(ss)) == 2) { levels(ss) <- c("resistant", "sensitive") }
	if(length(levels(ss)) == 1) { levels(ss) <- c("resistant") }
	return (ss)
}
auc.call3.cgp <- apply(auc.call3.cgp, 2, myfn)
dimnames(auc.call3.cgp) <- list(rownames(data.ge.cgp), rownames(druginfo))
auc.call3.ccle <- apply(auc.call3.ccle, 2, myfn)
dimnames(auc.call3.ccle) <- list(rownames(data.ge.ccle), rownames(druginfo))

pdf(file.path(saveres, "cgp_ccle_auc_call_extreme_common_paper.pdf"), width=8, height=9)
par(mfrow=c(5, 3), mar=c(5, 4, 4, 2) + 0.1)
mykap <- mykap2 <- NULL
for (i in 1:ncol(auc.call3.ccle)) {
	ff1 <- factor(x=auc.call3.ccle[ , i], levels=c("resistant", "intermediate", "sensitive"))
	ff2 <- factor(x=auc.call3.cgp[ , i], levels=c("resistant", "intermediate", "sensitive"))
  ## with intermediate
  levels(ff1) <- c("res", "inter", "sens")
  levels(ff2) <- c("res", "inter", "sens")
  nnn <- sum(complete.cases(ff1, ff2))
  tt <- table("CCLE"=ff1, "CGP"=ff2)
  tts <- vcd::assocstats(tt)
  ## check if all elements are on the diagonal
  ttt <- tt
  diag(ttt) <- 0
  if(sum(ttt) == 0) {
    rr <- rr.l <- rr.u <- 1
  } else {
    err <- try(rr <- epibasix::epiKappa(tt, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
		 rr <- rr.l <- rr.u <- NA
	 } else {
		rr.l <- rr$CIL
		rr.u <- rr$CIU
		rr <- rr$kappa
	}
  }
  mykap <- c(mykap, rr)
  ## without intermediate
  levels(ff1)[2] <- NA
  levels(ff2)[2] <- NA
  nnn2 <- sum(complete.cases(ff1, ff2))
  tt2 <- table("CCLE"=ff1, "CGP"=ff2)
  tts2 <- vcd::assocstats(tt2)
  ## check if all elements are on the diagonal
  ttt <- tt
  diag(ttt) <- 0
  if(sum(ttt) == 0) {
    rr2 <- rr2.l <- rr2.u <- 1
  } else {
    err <- try(rr2 <- epibasix::epiKappa(tt2, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
		 rr2 <- rr2.l <- rr2.u <- NA
	 } else {
		rr2.l <- rr2$CIL
		rr2.u <- rr2$CIU
		rr2 <- rr2$kappa
	}
  }
  mykap2 <- c(mykap2, rr2)
  plot.new()
  plotrix::addtable2plot(x=par("usr")[1] + 0.1, y=par("usr")[4] - 0.4, xjust=0, yjust=0, table=as.matrix(tt), bty="o", display.rownames=TRUE, display.colnames=TRUE, hlines=FALSE, vlines=FALSE, title="CCLE   vs   CGP", cex=1.25)
  title(main=sprintf("AUC sensitivity calling\n%s", gsub("drugid_", "", colnames(auc.call3.ccle)[i])), sub=sprintf("Kappa=%.2g, 95%%CI [%.2g,%.2g], p=%.1E\nKappa2=%.2g, 95%%CI [%.2g,%.2g], p=%.1E", rr, rr.l, rr.u, tts$chisq_tests[1,3], rr2, rr2.l, rr2.u, tts2$chisq_tests[1,3]))
  correlations[["auc.call"]][i, "drug.sensitivity"] <- rr
  ## statistics
  correlations.stats[["auc.call"]][i, ] <- c(rr, rr.l, rr.u, tts$chisq_tests[1,3], nnn)
}
names(mykap) <- names(mykap2) <- colnames(auc.call3.ccle)
dev.off()

########################
## gene drug associations from AUC extreme sensitivity calls, using only the common cell lines

message("Gene-drug association with tuned AUC sensitivity calls for common cell lines:")

myfn <- file.path(saveres, "cgp_ccle_auc_call_extreme_common_assoc.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC (CCLE)")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.call3.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.call3.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
		 auc <- factor(auc)
		levels(auc)[levels(auc) == "intermediate"] <- NA
      suppressWarnings(res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm))
      return(res)      
    }, data=data.ge.ccle, auc=auc.call3.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ge.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.ccle <- c(assoc.auc.ccle, list(mcres))
  }
  message("")
  names(assoc.auc.ccle) <- colnames(auc.call3.ccle)
  ## CGP
  message("Gene-drug association based on AUC (CGP)")
  assoc.auc.cgp <- NULL
  for(i in 1:ncol(auc.call3.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.call3.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      auc <- factor(auc)
		levels(auc)[levels(auc) == "intermediate"] <- NA
      suppressWarnings(res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm))
      return(res)      
    }, data=data.ge.cgp, auc=auc.call3.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ge.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.cgp <- c(assoc.auc.cgp, list(mcres))
  }
  message("")
  names(assoc.auc.cgp) <- colnames(auc.call3.cgp)
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.auc.ccle)) {
    tt <- cbind(assoc.auc.ccle[[i]], "EntrezID"=annot.ge[rownames(assoc.auc.ccle[[i]]), "EntrezID"], "symbol"=annot.ge[rownames(assoc.auc.ccle[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.ccle)
  # WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ccle_auc_call3_results_gene_drug_paper.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.auc.cgp)) {
    tt <- cbind(assoc.auc.cgp[[i]], "EntrezID"=annot.ge[rownames(assoc.auc.cgp[[i]]), "EntrezID"], "ymbol"=annot.ge[rownames(assoc.auc.cgp[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.cgp)
  # WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "cgp_auc_call3_results_gene_drug_paper.xls"), row.names=TRUE)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## sort gene drug assocation by quantile of p-values
auc.call3.assocs.conc <- NULL
for(i in 1:ncol(auc.call3.ccle)) {
  xx <- assoc.auc.cgp[[i]]
  yy <- assoc.auc.ccle[[i]]
  res <- t(sapply(quantile.cuts, function (qq, xx, yy, method=c("spearman", "cosine", "pearson")) {
    method <- match.arg(method)
    ## rank by p-values
    qqx <- quantile(x=xx[ , "pvalue"], probs=qq, na.rm=TRUE)
    qqy <- quantile(x=yy[ , "pvalue"], probs=qq, na.rm=TRUE)
    iix <- union(rownames(xx)[!is.na(xx[ , "pvalue"]) & xx[ , "pvalue"] <= qqx], rownames(yy)[!is.na(yy[ , "pvalue"]) & yy[ , "pvalue"] <= qqy])
    ## correlation
    switch(method,
      "spearman" = {
        nnn <- sum(complete.cases(xx[iix, "estimate"], yy[iix, "estimate"]))
        cc <- cor.test(x=xx[iix, "estimate"], y=yy[iix, "estimate"], method="spearman", use="pairwise.complete.obs", alternative="greater")
        cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
        res <- c(cc$estimate, cci[1], cci[2], nnn, cc$p.value)
      },
      "cosine" = {
        require(lsa)
        nnn <- sum(complete.cases(xx[iix, "estimate"], yy[iix, "estimate"]))
        cc <- lsa::cosine(x=xx[iix, "estimate"], y=yy[iix, "estimate"])
        ## TODO: compute confidence interval and significance using bootstrap
        cci <- c("lower"=NA, "upper"=NA, "p"=NA)
        res <- c(drop(cc), cci[1], cci[2], nnn, cci[3])
      },
      "pearson" = {
        nnn <- sum(complete.cases(xx[iix, "estimate"], yy[iix, "estimate"]))
        cc <- cor.test(x=xx[iix, "estimate"], y=yy[iix, "estimate"], method="pearson", use="pairwise.complete.obs", alternative="greater")
        res <- c(cc$estimate, cc$conf.int[1], cc$conf.int[2], nnn, cc$p.value)
      }
    )
    names(res) <- c("rho", "lower", "upper", "n", "p")   
    return (res)
  }, xx=xx, yy=yy, method=concordance.method))
  rownames(res) <- sprintf("PQUANTILE.%i", fdr.cuts * 100)
  auc.call3.assocs.conc <- c(auc.call3.assocs.conc, list(res))
}
names(auc.call3.assocs.conc) <- gsub("drugid_", "", colnames(auc.call3.ccle))
## barplot
pdf(file.path(saveres, "auc_call_extreme_common_assocs_conc_quantiles.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.call3.assocs.conc)) {
  xx <- auc.call3.assocs.conc[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.call3.assocs.conc[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.call3.assocs.conc[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.call3.assocs.conc[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.call3.assocs.conc[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.call3.assocs.conc)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  # text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Shared CGP and CCLE data\n\n\nP quantile cutoff", legend=sprintf("%s%%", gsub("PQUANTILE[.]", " ", rownames(auc.call3.assocs.conc[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()
# ## line plot
# plot(auc.call3.assocs.conc[[1]][ , "n"], auc.call3.assocs.conc[[1]][ , "rho"], type="l", lwd=1.5, ylim=c(0, 1), xlim=rev(range(auc.call3.assocs.conc[[1]][ , "n"])), col=drugs.color[names(auc.call3.assocs.conc)[1]], xlab="Number of gene-drug associations", ylab=concordance.method)
# for (i in 2:length(auc.call3.assocs.conc)) {
#   lines(auc.call3.assocs.conc[[i]][ , "n"], auc.call3.assocs.conc[[i]][ , "rho"], lwd=1.5, col=drugs.color[names(auc.call3.assocs.conc)[i]])
# }
# legend("topleft", legend=names(auc.call3.assocs.conc), bty="n", col=drugs.color[names(auc.call3.assocs.conc)], lwd=1.5, lty=1)

## sort gene drug assocation by fdr cutoffs
auc.call3.assocs.conc <- NULL
for(i in 1:ncol(auc.call3.ccle)) {
  xx <- assoc.auc.cgp[[i]]
  yy <- assoc.auc.ccle[[i]]
  res <- t(sapply(fdr.cuts, function (qq, xx, yy, method=c("spearman", "cosine", "pearson"), minsample=10) {
    method <- match.arg(method)
    ## rank by p-values
    iix <- union(rownames(xx)[!is.na(xx[ , "fdr"]) & xx[ , "fdr"] <= qq], rownames(yy)[!is.na(yy[ , "fdr"]) & yy[ , "fdr"] <= qq])
    nnn <- sum(complete.cases(xx[iix, "estimate"], yy[iix, "estimate"]))
    if (nnn >= minsample) {
      ## concordance
      switch(method,
        "spearman" = {
          cc <- cor.test(x=xx[iix, "estimate"], y=yy[iix, "estimate"], method="spearman", use="pairwise.complete.obs", alternative="greater")
          cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
          res <- c(cc$estimate, cci[1], cci[2], nnn, cc$p.value)
        },
        "cosine" = {
          require(lsa)
          cc <- lsa::cosine(x=xx[iix, "estimate"], y=yy[iix, "estimate"])
          ## TODO: compute confidence interval and significance using bootstrap
          cci <- c("lower"=NA, "upper"=NA, "p"=NA)
          res <- c(drop(cc), cci[1], cci[2], nnn, cci[3])
        },
        "pearson" = {
          cc <- cor.test(x=xx[iix, "estimate"], y=yy[iix, "estimate"], method="pearson", use="pairwise.complete.obs", alternative="greater")
          res <- c(cc$estimate, cc$conf.int[1], cc$conf.int[2], nnn, cc$p.value)
        }
      )
      names(res) <- c("rho", "lower", "upper", "n", "p")   
    } else {
      res <- c("rho"=NA, "lower"=NA, "upper"=NA, "n"=nnn, "p"=NA)
    }
    return (res)
  }, xx=xx, yy=yy, method=concordance.method, minsample=minsample))
  rownames(res) <- sprintf("FDR.%i", fdr.cuts * 100)
  auc.call3.assocs.conc <- c(auc.call3.assocs.conc, list(res))
}
names(auc.call3.assocs.conc) <- gsub("drugid_", "", colnames(auc.call3.ccle))
## barplot
pdf(file.path(saveres, "auc_call_extreme_common_assocs_conc_fdrs_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.call3.assocs.conc)) {
  xx <- auc.call3.assocs.conc[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.call3.assocs.conc[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.call3.assocs.conc[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.call3.assocs.conc[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.call3.assocs.conc[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.call3.assocs.conc)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Shared CGP and CCLE data\n\n\nFDR cutoff", legend=sprintf("%s%%", gsub("FDR[.]", " ", rownames(auc.call3.assocs.conc[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()

#################################################
#################################################
## all drugs and cell lines

load(file.path(saveres, "cdrug2_cgp_ccle_all.RData"))
tissue.ccle <- sampleinfo.ccle[ , "tissue.type"]
names(tissue.ccle) <- rownames(sampleinfo.ccle)
tissue.cgp <- sampleinfo.cgp[ , "tissue.type"]
names(tissue.cgp) <- rownames(sampleinfo.cgp)
gene.common <- intersect(colnames(data.ge.cgp), colnames(data.ge.ccle))
cell.common <- intersect(rownames(data.ge.cgp), rownames(data.ge.ccle))
data.ge.cgp <- data.ge.cgp[ , gene.common, drop=FALSE]
data.ge.ccle <- data.ge.ccle[ , gene.common, drop=FALSE]
drugsn <- rownames(drug.map)
annot.ge <- annot.ge[rownames(list.genes), , drop=FALSE]
data.ge.ccle <- data.ge.ccle[ , rownames(list.genes), drop=FALSE]
data.ge.cgp <- data.ge.cgp[ , rownames(list.genes), drop=FALSE]

########################
## boxplot of drug sensitivity (AUC)
########################

## drugs in common between CGP and CCLE
oo <- order(apply(drugpheno.cgp$AUC[cell.common, , drop=FALSE], 2, median, na.rm=TRUE), decreasing=FALSE)
commonix <- is.element(colnames(drugpheno.cgp$AUC[ , oo, drop=FALSE]), drug.map[ , "CGP"])
mycol <- rep("white", ncol(drugpheno.cgp$AUC))
mycol[commonix] <- "red"
pdf(file.path(saveres, "boxplot_auc_cgp_commoncells_drugs.pdf"), width=10, height=5)
par(las=3, mar=c(5, 4, 4, 2) + 0.1, xaxt="n", cex=0.8)
mp <- graphics::boxplot(drugpheno.cgp$AUC[cell.common, oo, drop=FALSE], outline=FALSE, ylab="AUC", main="Drug sensitivity (AUC)\nCGP", col=mycol, cex=0.5, ylim=c(0, 1))
axis(1, at=which(commonix), tick=TRUE, labels=T)
text(x=which(commonix) + 0.5, y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=rownames(drug.map), srt=45, xpd=NA, font=2, col="red")
# legend("topleft", legend=c("CGP", "CCLE"), col=c("lightgreen", "lightblue"), pch=15, bty="n")
dev.off()

## drugs in common between CCLE and CCLE
oo <- order(apply(drugpheno.ccle$AUC[cell.common, , drop=FALSE], 2, median, na.rm=TRUE), decreasing=FALSE)
commonix <- is.element(colnames(drugpheno.ccle$AUC[ , oo, drop=FALSE]), drug.map[ , "CCLE"])
mycol <- rep("white", ncol(drugpheno.ccle$AUC))
mycol[commonix] <- "red"
pdf(file.path(saveres, "boxplot_auc_ccle_commoncells_drugs.pdf"), width=10, height=5)
par(las=3, mar=c(5, 4, 4, 2) + 0.1, xaxt="n", cex=0.8)
mp <- graphics::boxplot(drugpheno.ccle$AUC[cell.common, oo, drop=FALSE], outline=FALSE, ylab="AUC", main="Drug sensitivity (AUC)\nCCLE", col=mycol, cex=0.5, ylim=c(0, 1))
axis(1, at=which(commonix), tick=TRUE, labels=T)
text(x=which(commonix) + 0.5, y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=rownames(drug.map), srt=45, xpd=NA, font=2, col="red")
# legend("topleft", legend=c("CCLE", "CCLE"), col=c("lightgreen", "lightblue"), pch=15, bty="n")
dev.off()

########################
## variance/MAD vs correlation
########################

drug.var.cgp <- apply(drugpheno.cgp$AUC[cell.common, drug.map[ , "CGP"]], 2, var, na.rm=TRUE)
drug.var.ccle <- apply(drugpheno.ccle$AUC[cell.common, drug.map[ , "CCLE"]], 2, var, na.rm=TRUE)
drug.mad.cgp <- apply(drugpheno.cgp$AUC[cell.common, drug.map[ , "CGP"]], 2, mad, na.rm=TRUE)
drug.mad.ccle <- apply(drugpheno.ccle$AUC[cell.common, drug.map[ , "CCLE"]], 2, mad, na.rm=TRUE)
names(drug.var.cgp) <- names(drug.var.ccle) <- names(drug.mad.cgp) <- names(drug.mad.ccle) <- drugsn
drug.cor <- correlations[["auc"]][ , "drug.sensitivity"]
drug.cor <- drug.cor[paste("drugid_", drugsn, sep="")]
drug.amcc <- mcc.auc.across[ , "mcc"]

pdf(file.path(saveres, "cgp_ccle_auc_var_vs_var.pdf"), width=5, height=5)
par(mar=c(5, 4, 1, 2) + 0.1, cex=0.7)
xxlim <- c(floor(min(drug.var.cgp, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.var.cgp, na.rm=TRUE) * 1200) / 1000)
yylim <- c(floor(min(drug.var.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.var.ccle, na.rm=TRUE) * 1200) / 1000)
llim <- c(min(xxlim[1], yylim[1]), max(xxlim[2], yylim[2]))
plot(x=drug.var.cgp, y=drug.var.ccle, xlim=llim, ylim=llim, pch=20, col=blues9[7], xlab="Variance of AUC in CGP", ylab="Variance of AUC in CCLE")
text(x=drug.var.cgp, y=drug.var.ccle, labels=drugsn, cex=0.7, font=1, srt=30, pos=4)
dev.off()

pdf(file.path(saveres, "cgp_ccle_auc_cor_vs_var.pdf"), width=10, height=5)
par(mfrow=c(1, 2), mar=c(5, 4, 1, 2) + 0.1, cex=0.7)
xxlim <- c(floor(min(drug.cor, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.cor, na.rm=TRUE) * 1200) / 1000)
## variance in cgp
yylim <- c(floor(min(drug.var.cgp, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.var.cgp, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.cor, y=drug.var.cgp, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="Correlation of AUC between CCLE and CGP", ylab="Variance of AUC in CGP")
text(x=drug.cor, y=drug.var.cgp, labels=drugsn, cex=0.7, font=1, srt=30, pos=4)
## variance in ccle
yylim <- c(floor(min(drug.var.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.var.ccle, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.cor, y=drug.var.ccle, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="Correlation of AUC between CCLE and CGP", ylab="Variance of AUC in CCLE")
text(x=drug.cor, y=drug.var.ccle, labels=drugsn, cex=0.7, font=1, srt=30, pos=4)
dev.off()

pdf(file.path(saveres, "cgp_ccle_auc_mad_vs_mad.pdf"), width=5, height=5)
par(mar=c(5, 4, 1, 2) + 0.1, cex=0.7)
xxlim <- c(floor(min(drug.mad.cgp, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.cgp, na.rm=TRUE) * 1200) / 1000)
yylim <- c(floor(min(drug.mad.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.ccle, na.rm=TRUE) * 1200) / 1000)
llim <- c(min(xxlim[1], yylim[1]), max(xxlim[2], yylim[2]))
plot(x=drug.mad.cgp, y=drug.mad.ccle, xlim=llim, ylim=llim, pch=20, col=blues9[7], xlab="MAD of AUC in CGP", ylab="MAD of AUC in CCLE")
abline(a=0, b=1, col="darkgrey", lwd=0.5)
abline(h=0.10, col="red", lty=2, lwd=0.5)
abline(v=0.10, col="red", lty=2, lwd=0.5)
text(x=drug.mad.cgp, y=drug.mad.ccle, labels=drugsn, cex=0.7, font=1, srt=45, pos=4)
dev.off()

pdf(file.path(saveres, "cgp_ccle_auc_cor_vs_mad.pdf"), width=10, height=5)
par(mfrow=c(1, 2), mar=c(5, 4, 1, 2) + 0.1, cex=0.7)
xxlim <- c(floor(min(drug.cor, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.cor, na.rm=TRUE) * 1200) / 1000)
## variability in cgp
cc <- cor.test(drug.mad.cgp, drug.cor, method=concordance.method, use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.cgp, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.cgp, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.cor, y=drug.mad.cgp, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="Correlation of AUC between CCLE and CGP", ylab="MAD of AUC in CGP")
text(x=drug.cor, y=drug.mad.cgp, labels=drugsn, cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=1, cex=1)
## variability in ccle
nnn <- sum(complete.cases(drug.mad.ccle, drug.cor))
cc <- cor.test(drug.mad.ccle, drug.cor, method=concordance.method, use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.ccle, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.cor, y=drug.mad.ccle, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="Correlation of AUC between CCLE and CGP", ylab="MAD of AUC in CCLE")
text(x=drug.cor, y=drug.mad.ccle, labels=drugsn, cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=1, cex=1)
dev.off()

pdf(file.path(saveres, "cgp_ccle_auc_amcc_vs_mad.pdf"), width=10, height=5)
par(mfrow=c(1, 2), mar=c(5, 4, 1, 2) + 0.1, cex=0.7)
xxlim <- c(floor(min(drug.amcc, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.amcc, na.rm=TRUE) * 1200) / 1000)
## variability in cgp
cc <- cor.test(drug.mad.cgp, drug.amcc, method=concordance.method, use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.cgp, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.cgp, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.amcc, y=drug.mad.cgp, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="AMCC of AUC between CCLE and CGP", ylab="MAD of AUC in CGP")
text(x=drug.amcc, y=drug.mad.cgp, labels=drugsn, cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=1, cex=1)
## variability in ccle
nnn <- sum(complete.cases(drug.mad.ccle, drug.amcc))
cc <- cor.test(drug.mad.ccle, drug.amcc, method=concordance.method, use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.ccle, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.amcc, y=drug.mad.ccle, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="AMCC of AUC between CCLE and CGP", ylab="MAD of AUC in CCLE")
text(x=drug.amcc, y=drug.mad.ccle, labels=drugsn, cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=1, cex=1)
dev.off()

########################
message("Correlation of tuned AUC sensitivity calls (15 most sensitive vs 55 most resistant for all cell lines):")

## new sensitivity call
auc.call3.cgp <- drugpheno.cgp$AUC[rownames(data.ge.cgp), drug.map[ , "CGP"], drop=FALSE]
auc.call3.ccle <- drugpheno.ccle$AUC[rownames(data.ge.ccle), drug.map[ , "CCLE"], drop=FALSE]

myfn <- function(x, n1=55, n2=15) {
	ss <- sort(x, na.last=NA)
	rr <- c(ss[n1], ss[length(ss) - n2 + 1])
	ss <- cut2(x=x, cuts=rr)
	if(length(levels(ss)) == 3) { levels(ss) <- c("resistant", "intermediate", "sensitive") }
	if(length(levels(ss)) == 2) { levels(ss) <- c("resistant", "sensitive") }
	if(length(levels(ss)) == 1) { levels(ss) <- c("resistant") }
	return (ss)
}
auc.call3.cgp <- apply(auc.call3.cgp, 2, myfn)
dimnames(auc.call3.cgp) <- list(rownames(data.ge.cgp), rownames(druginfo))
auc.call3.ccle <- apply(auc.call3.ccle, 2, myfn)
dimnames(auc.call3.ccle) <- list(rownames(data.ge.ccle), rownames(druginfo))

########################
## gene drug associations from AUC extreme sensitivity calls, using only the common cell lines

message("Gene-drug association with tuned AUC sensitivity calls for common cell lines:")

myfn <- file.path(saveres, "cgp_ccle_auc_call_extreme_all_assoc.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC (CCLE)")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.call3.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.call3.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
		 auc <- factor(auc)
		levels(auc)[levels(auc) == "intermediate"] <- NA
      suppressWarnings(res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm))
      return(res)      
    }, data=data.ge.ccle, auc=auc.call3.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ge.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.ccle <- c(assoc.auc.ccle, list(mcres))
  }
  message("")
  names(assoc.auc.ccle) <- colnames(auc.call3.ccle)
  ## CGP
  message("Gene-drug association based on AUC (CGP)")
  assoc.auc.cgp <- NULL
  for(i in 1:ncol(auc.call3.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.call3.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      auc <- factor(auc)
		levels(auc)[levels(auc) == "intermediate"] <- NA
      suppressWarnings(res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm))
      return(res)      
    }, data=data.ge.cgp, auc=auc.call3.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ge.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.cgp <- c(assoc.auc.cgp, list(mcres))
  }
  message("")
  names(assoc.auc.cgp) <- colnames(auc.call3.cgp)
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.auc.ccle)) {
    tt <- cbind(assoc.auc.ccle[[i]], "EntrezID"=annot.ge[rownames(assoc.auc.ccle[[i]]), "EntrezID"], "symbol"=annot.ge[rownames(assoc.auc.ccle[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.ccle)
  # WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ccle_auc_call3_results_gene_drug_paper.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.auc.cgp)) {
    tt <- cbind(assoc.auc.cgp[[i]], "EntrezID"=annot.ge[rownames(assoc.auc.cgp[[i]]), "EntrezID"], "ymbol"=annot.ge[rownames(assoc.auc.cgp[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.cgp)
  # WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "cgp_auc_call3_results_gene_drug_paper.xls"), row.names=TRUE)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## sort gene drug assocation by quantile of p-values
auc.call3.assocs.conc <- NULL
for(i in 1:ncol(auc.call3.ccle)) {
  xx <- assoc.auc.cgp[[i]]
  yy <- assoc.auc.ccle[[i]]
  res <- t(sapply(quantile.cuts, function (qq, xx, yy, method=c("spearman", "cosine", "pearson")) {
    method <- match.arg(method)
    ## rank by p-values
    qqx <- quantile(x=xx[ , "pvalue"], probs=qq, na.rm=TRUE)
    qqy <- quantile(x=yy[ , "pvalue"], probs=qq, na.rm=TRUE)
    iix <- union(rownames(xx)[!is.na(xx[ , "pvalue"]) & xx[ , "pvalue"] <= qqx], rownames(yy)[!is.na(yy[ , "pvalue"]) & yy[ , "pvalue"] <= qqy])
    ## correlation
    switch(method,
      "spearman" = {
        nnn <- sum(complete.cases(xx[iix, "estimate"], yy[iix, "estimate"]))
        cc <- cor.test(x=xx[iix, "estimate"], y=yy[iix, "estimate"], method="spearman", use="pairwise.complete.obs", alternative="greater")
        cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
        res <- c(cc$estimate, cci[1], cci[2], nnn, cc$p.value)
      },
      "cosine" = {
        require(lsa)
        nnn <- sum(complete.cases(xx[iix, "estimate"], yy[iix, "estimate"]))
        cc <- lsa::cosine(x=xx[iix, "estimate"], y=yy[iix, "estimate"])
        ## TODO: compute confidence interval and significance using bootstrap
        cci <- c("lower"=NA, "upper"=NA, "p"=NA)
        res <- c(drop(cc), cci[1], cci[2], nnn, cci[3])
      },
      "pearson" = {
        nnn <- sum(complete.cases(xx[iix, "estimate"], yy[iix, "estimate"]))
        cc <- cor.test(x=xx[iix, "estimate"], y=yy[iix, "estimate"], method="pearson", use="pairwise.complete.obs", alternative="greater")
        res <- c(cc$estimate, cc$conf.int[1], cc$conf.int[2], nnn, cc$p.value)
      }
    )
    names(res) <- c("rho", "lower", "upper", "n", "p")   
    return (res)
  }, xx=xx, yy=yy, method=concordance.method))
  rownames(res) <- sprintf("PQUANTILE.%i", fdr.cuts * 100)
  auc.call3.assocs.conc <- c(auc.call3.assocs.conc, list(res))
}
names(auc.call3.assocs.conc) <- gsub("drugid_", "", colnames(auc.call3.ccle))
## barplot
pdf(file.path(saveres, "auc_call_extreme_all_assocs_conc_quantiles.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.call3.assocs.conc)) {
  xx <- auc.call3.assocs.conc[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.call3.assocs.conc[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.call3.assocs.conc[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.call3.assocs.conc[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.call3.assocs.conc[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.call3.assocs.conc)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  # text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Full CGP and CCLE data\n\n\nP quantile cutoff", legend=sprintf("%s%%", gsub("PQUANTILE[.]", " ", rownames(auc.call3.assocs.conc[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()
# ## line plot
# plot(auc.call3.assocs.conc[[1]][ , "n"], auc.call3.assocs.conc[[1]][ , "rho"], type="l", lwd=1.5, ylim=c(0, 1), xlim=rev(range(auc.call3.assocs.conc[[1]][ , "n"])), col=drugs.color[names(auc.call3.assocs.conc)[1]], xlab="Number of gene-drug associations", ylab=concordance.method)
# for (i in 2:length(auc.call3.assocs.conc)) {
#   lines(auc.call3.assocs.conc[[i]][ , "n"], auc.call3.assocs.conc[[i]][ , "rho"], lwd=1.5, col=drugs.color[names(auc.call3.assocs.conc)[i]])
# }
# legend("topleft", legend=names(auc.call3.assocs.conc), bty="n", col=drugs.color[names(auc.call3.assocs.conc)], lwd=1.5, lty=1)

## sort gene drug assocation by fdr cutoffs
auc.call3.assocs.conc <- NULL
for(i in 1:ncol(auc.call3.ccle)) {
  xx <- assoc.auc.cgp[[i]]
  yy <- assoc.auc.ccle[[i]]
  res <- t(sapply(fdr.cuts, function (qq, xx, yy, method=c("spearman", "cosine", "pearson"), minsample=10) {
    method <- match.arg(method)
    ## rank by p-values
    iix <- union(rownames(xx)[!is.na(xx[ , "fdr"]) & xx[ , "fdr"] <= qq], rownames(yy)[!is.na(yy[ , "fdr"]) & yy[ , "fdr"] <= qq])
    nnn <- sum(complete.cases(xx[iix, "estimate"], yy[iix, "estimate"]))
    if (nnn >= minsample) {
      ## concordance
      switch(method,
        "spearman" = {
          cc <- cor.test(x=xx[iix, "estimate"], y=yy[iix, "estimate"], method="spearman", use="pairwise.complete.obs", alternative="greater")
          cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
          res <- c(cc$estimate, cci[1], cci[2], nnn, cc$p.value)
        },
        "cosine" = {
          require(lsa)
          cc <- lsa::cosine(x=xx[iix, "estimate"], y=yy[iix, "estimate"])
          ## TODO: compute confidence interval and significance using bootstrap
          cci <- c("lower"=NA, "upper"=NA, "p"=NA)
          res <- c(drop(cc), cci[1], cci[2], nnn, cci[3])
        },
        "pearson" = {
          cc <- cor.test(x=xx[iix, "estimate"], y=yy[iix, "estimate"], method="pearson", use="pairwise.complete.obs", alternative="greater")
          res <- c(cc$estimate, cc$conf.int[1], cc$conf.int[2], nnn, cc$p.value)
        }
      )
      names(res) <- c("rho", "lower", "upper", "n", "p")   
    } else {
      res <- c("rho"=NA, "lower"=NA, "upper"=NA, "n"=nnn, "p"=NA)
    }
    return (res)
  }, xx=xx, yy=yy, method=concordance.method, minsample=minsample))
  rownames(res) <- sprintf("FDR.%i", fdr.cuts * 100)
  auc.call3.assocs.conc <- c(auc.call3.assocs.conc, list(res))
}
names(auc.call3.assocs.conc) <- gsub("drugid_", "", colnames(auc.call3.ccle))
## barplot
pdf(file.path(saveres, "auc_call_extreme_all_assocs_conc_fdrs_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.call3.assocs.conc)) {
  xx <- auc.call3.assocs.conc[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.call3.assocs.conc[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.call3.assocs.conc[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.call3.assocs.conc[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.call3.assocs.conc[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.call3.assocs.conc)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Full CGP and CCLE data\n\n\nFDR cutoff", legend=sprintf("%s%%", gsub("FDR[.]", " ", rownames(auc.call3.assocs.conc[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()






