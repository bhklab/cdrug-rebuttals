########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## April 8, 2014
## March 1, 2015
########################

## TODO
## Fig 1: boxplot of correlation across and between cell lines
## Fig 2: gene-drug association using Geeleher's approach -> logistic regression with 15 most sensitive vs 55 most resistant
## Fig 3: new plot for published AUC with red points
## Explanation of AMCC in SI
## Fig 4: AMCC
## Fig 5: AMCC for CGP replicates camptothecin
## Fig 6: boxplot for AMCC for (full) gene expression and published AUC
## find biomarkers reported by Hunag in our suppl tables


require(amap) || stop("Library amap is not available!")
require(vcd) || stop("Library vcd is not available!")
require(gplots) || stop("Library gplots is not available!")
require(WriteXLS) || stop("Library WriteXLS is not available!")
require(xtable) || stop("Library xtable is not available!")
require(epibasix) || stop("Library gplots is not available!")
require(OptimalCutpoints) || stop("Library OptimalCutpoints is not available!")

## load data shared between CGP and CCLE
load(file.path(saveres, "cdrug2_cgp_ccle_common.RData"))
drugsn <- gsub("drugid_", "", rownames(druginfo))
drugs.color <- rainbow(length(drugsn), v=0.9)
names(drugs.color) <- drugsn
nature2013.common.cellines <- read.csv(file=file.path("code", "HaibeKains_Nature_2013_common_cellines.csv"))[ , 1]
sampleinfo <- sampleinfo.cgp
l1000.genes <- l1000.genes[intersect(rownames(annot.ge), rownames(l1000.genes)), , drop=FALSE]

########################
## tissue types
########################

tissue.cgp <- as.character(sampleinfo.cgp[ , "tissue.type"])
tissue.ccle <- as.character(sampleinfo.ccle[ , "tissue.type"])
tissue <- tissue.cgp
tissuen <- sort(unique(as.character(tissue)))
drugsn <- rownames(druginfo)

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
    xx <- -log10(drugpheno.cgp$IC50[, i] / 10^6)
    yy <- -log10(drugpheno.ccle$IC50[, i] / 10^6)
    ccix <- complete.cases(xx, yy)
    xx2 <- rank(-xx[ccix], ties.method="first")
    yy2 <- rank(-yy[ccix], ties.method="first")
    mm <- NULL
    for(j in 1:(min(max(xx2), max(yy2)) - 1)) {
      xx3 <- factor(ifelse (xx2 <= j, "sensitive", "resistant"))
      yy3 <- factor(ifelse (yy2 <= j, "sensitive", "resistant"))
      tt <- table("CGP"=xx3, "CCLE"=yy3)
      ## res <- mcc(x=xx3, y=yy3, nperm=1000, nthread=nbcore)
      res <- mcc(x=xx3, y=yy3, nperm=0, nthread=1)
      mm <- rbind(mm, res)
    }
    ## remove extreme indices
    rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 1):nrow(mm))
    mccix <- max(which(mm[-rmix, "mcc"] == max(mm[-rmix, "mcc"], na.rm=TRUE))) + (min.cat - 1)
    ## compute significance only for the AMCC
    xx3 <- factor(ifelse (xx2 <= mccix, "sensitive", "resistant"))
    yy3 <- factor(ifelse (yy2 <= mccix, "sensitive", "resistant"))
    mm[mccix, "p"] <- mcc(x=xx3, y=yy3, nperm=10^3, nthread=nbcore)["p"]
    ## bonferronni correction
    # mm[mccix, "p"] <- mm[mccix, "p"] * length(xx3)
    if (!is.na(mm[mccix, "p"]) && mm[mccix, "p"] > 1) { mm[mccix, "p"] <- 1 }
    mcc.ic50 <- rbind(mcc.ic50, c(mm[mccix, ], "n1"=mccix, "n2"=nrow(mm) - mccix, "n"=nrow(mm)))
    par(mfrow=c(1, 2))
    ## scatterplot with sperman correlation
    myScatterPlot(xx, yy, main=drugn, xlab="IC50 (CGP)", ylab="IC50 (CCLE)")
    rs <- cor.test(x=xx, y=yy, method="spearman", use="pairwise.complete.obs")
	 oo <- sort(xx[ccix], decreasing=TRUE)
	 abline(v=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
	 oo <- sort(yy[ccix], decreasing=TRUE)
	 abline(h=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
    legend("topleft", legend=c(sprintf("Rs=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# cell lines=%i", sum(ccix))), col="white", pch=0, bty="n")
    ## mcc plot
    plot(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], ylim=c(min(mm[-rmix, "mcc"]), 1), xlab="# sensitive cell lines", ylab="AMCC", pch=20, col="lightgrey", main=drugn)
    lines(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], col="darkgrey")
    legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "mcc"], mm[mccix, "p"]), sprintf("# sensitive cell lines=%i", mccix)), col="white", pch=0, bty="n")
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
    xx2 <- rank(-xx[ccix], ties.method="first")
    yy2 <- rank(-yy[ccix], ties.method="first")
    mm <- NULL
    for(j in 1:(min(max(xx2), max(yy2)) - 1)) {
      xx3 <- factor(ifelse (xx2 <= j, "sensitive", "resistant"))
      yy3 <- factor(ifelse (yy2 <= j, "sensitive", "resistant"))
      # tt <- table("CGP"=xx3, "CCLE"=yy3)
      ## res <- mcc(x=xx3, y=yy3, nperm=1000, nthread=nbcore)
      res <- mcc(x=xx3, y=yy3, nperm=0, nthread=1)
      mm <- rbind(mm, res)
    }
    ## remove extreme indices
    rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 1):nrow(mm))
    mccix <- max(which(mm[-rmix, "mcc"] == max(mm[-rmix, "mcc"], na.rm=TRUE))) + (min.cat - 1)
    ## compute significance only for the AMCC
    xx3 <- factor(ifelse (xx2 <= mccix, "sensitive", "resistant"))
    yy3 <- factor(ifelse (yy2 <= mccix, "sensitive", "resistant"))
    mm[mccix, "p"] <- mcc(x=xx3, y=yy3, nperm=10^3, nthread=nbcore)["p"]
    ## bonferronni correction
    # mm[mccix, "p"] <- mm[mccix, "p"] * length(xx3)
    if (!is.na(mm[mccix, "p"]) && mm[mccix, "p"] > 1) { mm[mccix, "p"] <- 1 }
    mcc.auc <- rbind(mcc.auc, c(mm[mccix, ], "n1"=mccix, "n2"=nrow(mm) - mccix, "n"=nrow(mm)))
    par(mfrow=c(1, 2))
    ## scatterplot with sperman correlation
    myScatterPlot(xx, yy, main=drugn, xlab="AUC (CGP)", ylab="AUC (CCLE)")
    rs <- cor.test(x=xx, y=yy, method="spearman", use="pairwise.complete.obs")
	 oo <- sort(xx[ccix], decreasing=TRUE)
	 abline(v=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
	 oo <- sort(yy[ccix], decreasing=TRUE)
	 abline(h=(oo[mccix] + oo[mccix + 1]) / 2, col="red", lty=2)
    legend("topleft", legend=c(sprintf("Rs=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# cell lines=%i", sum(ccix))), col="white", pch=0, bty="n")
    ## mcc plot
    plot(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], ylim=c(min(mm[-rmix, "mcc"]), 1), xlab="# sensitive cell lines", ylab="AMCC", pch=20, col="lightgrey", main=drugn)
    lines(x=1:nrow(mm[-rmix, ]), y=mm[-rmix, "mcc"], col="darkgrey")
    legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "mcc"], mm[mccix, "p"]), sprintf("# sensitive cell lines=%i", mccix)), col="white", pch=0, bty="n")
  }
  rownames(mcc.auc) <- rownames(druginfo)
  dev.off()
  save(list=c("mcc.auc"), compress=TRUE, file=myfn)
} else { load(myfn) }
mcc.auc.across <- mcc.auc

## AMCC for gene expression across cell lines
myfn <- file.path(saveres, "ge_cgp_ccle_amcc_across.RData")
if (!file.exists(myfn)) {
  splitix <- parallel::splitIndices(nx=nrow(l1000.genes), ncl=nthread)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(x, xx, yy, min.cat) {
    res <- t(sapply(x, function(x, xx, yy, min.cat) {
      xx <- xx[ , x]
      yy <- yy[ , x]
      ccix <- complete.cases(xx, yy)
      xx2 <- rank(-xx[ccix], ties.method="first")
      yy2 <- rank(-yy[ccix], ties.method="first")
      mm <- NULL
      for(j in 1:(min(max(xx2), max(yy2)) - 1)) {
        xx3 <- factor(ifelse (xx2 <= j, "sensitive", "resistant"))
        yy3 <- factor(ifelse (yy2 <= j, "sensitive", "resistant"))
        # tt <- table("CGP"=xx3, "CCLE"=yy3)
        ## res <- mcc(x=xx3, y=yy3, nperm=1000, nthread=nbcore)
        res <- mcc(x=xx3, y=yy3, nperm=0, nthread=1)
        mm <- rbind(mm, res)
      }
      ## remove extreme indices
      rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 1):nrow(mm))
      mccix <- max(which(mm[-rmix, "mcc"] == max(mm[-rmix, "mcc"], na.rm=TRUE))) + (min.cat - 1)
      ## compute significance only for the AMCC
      # xx3 <- factor(ifelse (xx2 <= mccix, "sensitive", "resistant"))
      # yy3 <- factor(ifelse (yy2 <= mccix, "sensitive", "resistant"))
      # mm[mccix, "p"] <- mcc(x=xx3, y=yy3, nperm=10^3, nthread=nbcore)["p"]
      ## bonferronni correction
      # mm[mccix, "p"] <- mm[mccix, "p"] * length(xx3)
      # if (!is.na(mm[mccix, "p"]) && mm[mccix, "p"] > 1) { mm[mccix, "p"] <- 1 }
      return (mm[mccix, ])
    }, xx=xx, yy=yy, min.cat=min.cat))
    return (res)
  }, xx=data.ge.cgp[ , rownames(l1000.genes), drop=FALSE], yy=data.ge.ccle[ , rownames(l1000.genes), drop=FALSE], min.cat=min.cat)
  mcc.ge <- do.call(rbind, mcres)
  # rownames(mcc.ge) <- rownames(annot.ge)
  rownames(mcc.ge) <- rownames(l1000.genes)
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
    ccix <- complete.cases(xx, yy)
    if (sum(ccix) >= (2 * min.cat)) {
      xx2 <- rank(-xx[ccix], ties.method="first")
      yy2 <- rank(-yy[ccix], ties.method="first")
      mm <- NULL
      for(j in 1:(min(max(xx2), max(yy2)) - 1)) {
        xx3 <- factor(ifelse (xx2 <= j, "sensitive", "resistant"))
        yy3 <- factor(ifelse (yy2 <= j, "sensitive", "resistant"))
        # tt <- table("CGP"=xx3, "CCLE"=yy3)
        ## res <- mcc(x=xx3, y=yy3, nperm=1000, nthread=nbcore)
        res <- mcc(x=xx3, y=yy3, nperm=0, nthread=1)
        mm <- rbind(mm, res)
      }
      ## remove extreme indices
      rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
      mccix <- max(which(mm[-rmix, "mcc"] == max(mm[-rmix, "mcc"], na.rm=TRUE))) + (min.cat - 1)
      ## compute significance only for the AMCC
      xx3 <- factor(ifelse (xx2 <= mccix, "sensitive", "resistant"))
      yy3 <- factor(ifelse (yy2 <= mccix, "sensitive", "resistant"))
      mm[mccix, "p"] <- mcc(x=xx3, y=yy3, nperm=10^3, nthread=nbcore)["p"]
      ## bonferronni correction
      # mm[mccix, "p"] <- mm[mccix, "p"] * length(xx3)
      if (!is.na(mm[mccix, "p"]) && mm[mccix, "p"] > 1) { mm[mccix, "p"] <- 1 }
      mcc.ic50 <- rbind(mcc.ic50, c(mm[mccix, ], "n1"=mccix, "n2"=nrow(mm) - mccix, "n"=nrow(mm)))
      # par(mfrow=c(1, 2))
      # ## scatterplot with sperman correlation
      # myScatterPlot(xx, yy, main=celln, xlab="IC50 (CGP)", ylab="IC50 (CCLE)")
      # rs <- cor.test(x=xx, y=yy, method="spearman", use="pairwise.complete.obs")
      # legend("topright", legend=c(sprintf("Rs=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# drugs=%i", sum(ccix))), col="white", pch=0, bty="n")
      # ## mcc plot
      # plot(x=(1:nrow(mm))[-rmix], y=mm[-rmix, "mcc"], ylim=c(min(mm[-rmix, "mcc"]), 1), xlab="# drugs", ylab="MCC", pch=20, col="lightgrey", main=celln)
      # lines(x=(1:nrow(mm))[-rmix], y=mm[-rmix, "mcc"], col="grey")
      # legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "mcc"], mm[mccix, "p"]), sprintf("# effective drugs=%i", mccix)), col="white", pch=0, bty="n")
    } else {
      mcc.ic50 <- rbind(mcc.ic50, c("mcc"=NA, "p"=NA, "n1"=0, "n2"=0, "n"=sum(ccix)))
    }
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
    ccix <- complete.cases(xx, yy)
    if (sum(ccix) >= (2 * min.cat)) {
      xx2 <- rank(-xx[ccix], ties.method="first")
      yy2 <- rank(-yy[ccix], ties.method="first")
      mm <- NULL
      for(j in 1:(min(max(xx2), max(yy2)) - 1)) {
        xx3 <- factor(ifelse (xx2 <= j, "sensitive", "resistant"))
        yy3 <- factor(ifelse (yy2 <= j, "sensitive", "resistant"))
        # tt <- table("CGP"=xx3, "CCLE"=yy3)
        ## res <- mcc(x=xx3, y=yy3, nperm=1000, nthread=nbcore)
        res <- mcc(x=xx3, y=yy3, nperm=0, nthread=1)
        mm <- rbind(mm, res)
      }
      ## remove extreme indices
      rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
      mccix <- max(which(mm[-rmix, "mcc"] == max(mm[-rmix, "mcc"], na.rm=TRUE))) + (min.cat - 1)
      ## compute significance only for the AMCC
      xx3 <- factor(ifelse (xx2 <= mccix, "sensitive", "resistant"))
      yy3 <- factor(ifelse (yy2 <= mccix, "sensitive", "resistant"))
      mm[mccix, "p"] <- mcc(x=xx3, y=yy3, nperm=10^3, nthread=nbcore)["p"]
      ## bonferronni correction
      # mm[mccix, "p"] <- mm[mccix, "p"] * length(xx3)
      if (!is.na(mm[mccix, "p"]) && mm[mccix, "p"] > 1) { mm[mccix, "p"] <- 1 }
      mcc.auc <- rbind(mcc.auc, c(mm[mccix, ], "n1"=mccix, "n2"=nrow(mm) - mccix, "n"=nrow(mm)))
      # par(mfrow=c(1, 2))
      # ## scatterplot with sperman correlation
      # myScatterPlot(xx, yy, main=celln, xlab="AUC (CGP)", ylab="AUC (CCLE)")
      # rs <- cor.test(x=xx, y=yy, method="spearman", use="pairwise.complete.obs")
      # legend("topright", legend=c(sprintf("Rs=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# drugs=%i", sum(ccix))), col="white", pch=0, bty="n")
      # ## mcc plot
      # plot(x=(1:nrow(mm))[-rmix], y=mm[-rmix, "mcc"], ylim=c(min(mm[-rmix, "mcc"]), 1), xlab="# drugs", ylab="MCC", pch=20, col="lightgrey", main=celln)
      # lines(x=(1:nrow(mm))[-rmix], y=mm[-rmix, "mcc"], col="grey")
      # legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "mcc"], mm[mccix, "p"]), sprintf("# effective drugs=%i", mccix)), col="white", pch=0, bty="n")
    } else {
      mcc.auc <- rbind(mcc.auc, c("mcc"=NA, "p"=NA, "n1"=0, "n2"=0, "n"=sum(ccix)))
    }
  }
  rownames(mcc.auc) <- rownames(sampleinfo)
  save(list=c("mcc.auc"), compress=TRUE, file=myfn)
} else { load(myfn) }
mcc.auc.between <- mcc.auc

## AMCC for gene expression between cell lines
myfn <- file.path(saveres, "ge_cgp_ccle_amcc_between.RData")
if (!file.exists(myfn)) {
  mcc.ge <- NULL
  for(i in 1:nrow(sampleinfo)) {
    celln <- rownames(sampleinfo)[i]
    # message(sprintf("compute AMCC for %s", celln))
    xx <- data.ge.cgp[i, rownames(l1000.genes)]
    yy <- data.ge.ccle[i, rownames(l1000.genes)]
    ccix <- complete.cases(xx, yy)
    if (sum(ccix) >= (2 * min.cat)) {
      xx2 <- rank(-xx[ccix], ties.method="first")
      yy2 <- rank(-yy[ccix], ties.method="first")
      mm <- NULL
      for(j in 1:(min(max(xx2), max(yy2)) - 1)) {
        xx3 <- factor(ifelse (xx2 <= j, "sensitive", "resistant"))
        yy3 <- factor(ifelse (yy2 <= j, "sensitive", "resistant"))
        # tt <- table("CGP"=xx3, "CCLE"=yy3)
        ## res <- mcc(x=xx3, y=yy3, nperm=1000, nthread=nbcore)
        res <- mcc(x=xx3, y=yy3, nperm=0, nthread=1)
        mm <- rbind(mm, res)
      }
      ## remove extreme indices
      rmix <- c(1:(min.cat - 1), (nrow(mm) - min.cat + 2):nrow(mm))
      mccix <- max(which(mm[-rmix, "mcc"] == max(mm[-rmix, "mcc"], na.rm=TRUE))) + (min.cat - 1)
      ## compute significance only for the AMCC
      xx3 <- factor(ifelse (xx2 <= mccix, "sensitive", "resistant"))
      yy3 <- factor(ifelse (yy2 <= mccix, "sensitive", "resistant"))
      mm[mccix, "p"] <- mcc(x=xx3, y=yy3, nperm=10^3, nthread=nbcore)["p"]
      ## bonferronni correction
      # mm[mccix, "p"] <- mm[mccix, "p"] * length(xx3)
      if (!is.na(mm[mccix, "p"]) && mm[mccix, "p"] > 1) { mm[mccix, "p"] <- 1 }
      mcc.ge <- rbind(mcc.ge, c(mm[mccix, ], "n1"=mccix, "n2"=nrow(mm) - mccix, "n"=nrow(mm)))
      # par(mfrow=c(1, 2))
      # ## scatterplot with sperman correlation
      # myScatterPlot(xx, yy, main=celln, xlab="AUC (CGP)", ylab="AUC (CCLE)")
      # rs <- cor.test(x=xx, y=yy, method="spearman", use="pairwise.complete.obs")
      # legend("topright", legend=c(sprintf("Rs=%.2g, p=%.1E", rs$estimate, rs$p.value), sprintf("# drugs=%i", sum(ccix))), col="white", pch=0, bty="n")
      # ## mcc plot
      # plot(x=(1:nrow(mm))[-rmix], y=mm[-rmix, "mcc"], ylim=c(min(mm[-rmix, "mcc"]), 1), xlab="# drugs", ylab="MCC", pch=20, col="lightgrey", main=celln)
      # lines(x=(1:nrow(mm))[-rmix], y=mm[-rmix, "mcc"], col="grey")
      # legend("topright", legend=c(sprintf("AMCC=%.2g, p=%.1E", mm[mccix, "mcc"], mm[mccix, "p"]), sprintf("# effective drugs=%i", mccix)), col="white", pch=0, bty="n")
    } else {
      mcc.ge <- rbind(mcc.ge, c("mcc"=NA, "p"=NA, "n1"=0, "n2"=0, "n"=sum(ccix)))
    }
  }
  rownames(mcc.ge) <- rownames(sampleinfo)
  save(list=c("mcc.ge"), compress=TRUE, file=myfn)
} else { load(myfn) }
mcc.ge.between <- mcc.ge


## boxplot of MCC
pdf(file.path(saveres, "boxplot_amcc_across.pdf"), width=6, height=6)
ll <- list("GE"=mcc.ge.across[ , "mcc"], "AUC"=mcc.auc.across[ , "mcc"], "IC50"=mcc.ic50.across[ , "mcc"])
kt <- kruskal.test(x=ll)
wt1 <- wilcox.test(x=ll$GE, y=ll$AUC)
wt2 <- wilcox.test(x=ll$GE, y=ll$IC50)
boxplot(ll, ylab="AMCC", main="Concordance across cell lines", pch=20, col="lightgrey", sub=sprintf("GE vs. AUC=%.1E\nGE vs. IC50=%.1E", wt1$p.value, wt2$p.value), border=c("black", "red", "red"))
dev.off()

pdf(file.path(saveres, "boxplot_amcc_between.pdf"), width=6, height=6)
ll <- list("GE"=mcc.ge.between[ , "mcc"], "AUC"=mcc.auc.between[ , "mcc"], "IC50"=mcc.ic50.between[ , "mcc"])
kt <- kruskal.test(x=ll)
wt1 <- wilcox.test(x=ll$GE, y=ll$AUC)
wt2 <- wilcox.test(x=ll$GE, y=ll$IC50)
boxplot(ll, ylab="AMCC", main="Concordance between cell lines", pch=20, col="lightgrey", sub=sprintf("GE vs. AUC=%.1E\nGE vs. IC50=%.1E", wt1$p.value, wt2$p.value), border=c("black", "red", "red"))
dev.off()

## create tables summarizing all the correlations
tt <- matrix(NA, nrow=length(drugsn), ncol=2, dimnames=list(drugsn, c("drug.sensitivity", "gene.drug")))
correlations <- list("ic50"=tt, "ic50.call"=tt, "auc"=tt, "auc.call"=tt)
## correlation statistics
tt <- matrix(NA, nrow=length(drugsn), ncol=5, dimnames=list(drugsn, c("rho", "lower", "upper", "p", "n")))
correlations.stats <- list("ic50"=tt, "ic50.call"=tt, "auc"=tt, "auc.call"=tt)

## consistency between IC50 with spearman correlation
pdf(file.path(saveres, "cgp_ccle_scatterplot_ic50_spearman_pres.pdf"), height=9, width=16)
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
    cc <- cor.test(xx, yy, method="spearman", use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  par(mar=c(4, 4, 3, 1) + 0.1)
  mycol <- rep(blues9[7], length(xx))
  names(mycol) <- names(xx)
  mycol[!is.element(names(mycol), nature2013.common.cellines)] <- "red3"
  myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 10, "-log10 IC50 (CGP)", ""), ylab=ifelse((i %% 5) == 1, "-log10 IC50 (CCLE)", ""), main=drugn, xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75, col=mycol)
  abline(a=0, b=1, col="black")
  legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2, cex=1)
  correlations[["ic50"]][i, "drug.sensitivity"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  correlations.stats[["ic50"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
}
dev.off()
## scatter plot with barplot
pdf(file.path(saveres, "cgp_ccle_scatterplot_ic50_spearman_paper.pdf"), height=14, width=14)
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
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.stats[["ic50"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["ic50"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["ic50"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["ic50"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(correlations.stats[["auc"]])
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 0.71)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
# plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

## consistency between AUC with spearman correlation
pdf(file.path(saveres, "cgp_ccle_scatterplot_auc_spearman_pres.pdf"), height=9, width=16)
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
    cc <- cor.test(xx, yy, method="spearman", use="complete.obs", alternative="greater")
  } else {
    cc <- list("estimate"=NA, "p.value"=NA)
  }
  par(mar=c(4, 4, 3, 1) + 0.1)
  mycol <- rep(blues9[7], length(xx))
  names(mycol) <- names(xx)
  mycol[!is.element(names(mycol), nature2013.common.cellines)] <- "red3"
  myScatterPlot(x=xx, y=yy, xlab=ifelse(i > 10, "AUC (CGP)", ""), ylab=ifelse((i %% 5) == 1, "AUC (CCLE)", ""), main=drugn, xlim=xxlim, ylim=yylim, pch=16, method="transparent", transparency=0.75, col=mycol)
  abline(a=0, b=1, col="black")
  legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2, cex=1)
  correlations[["auc"]][i, "drug.sensitivity"] <- cc$estimate
  ## correlation statistics
  cci <- spearmanCI(x=cc$estimate, n=nnn, alpha=0.05)
  correlations.stats[["auc"]][i, ] <- c(cc$estimate, cci[1], cci[2], cc$p.value, nnn)
}
dev.off()

## scatter plot with barplot
pdf(file.path(saveres, "cgp_ccle_scatterplot_auc_spearman_paper.pdf"), height=14, width=14)
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
par(mar=c(6, 4, 3, 1) + 0.1, xaxt="n", las=1)
xx <- correlations.stats[["auc"]][ , "rho"]
xx[!is.na(xx) & xx < 0] <- 0
ll <- correlations.stats[["auc"]][ , "lower"]
ll[!is.na(ll) & ll < 0] <- 0
uu <- correlations.stats[["auc"]][ , "upper"]
uu[!is.na(uu) & uu < 0] <- 0
uu[!is.na(uu) & uu > 1] <- 1
pp <- correlations.stats[["auc"]][ , "p"]
names(xx) <- names(ll) <- names(uu) <- names(pp) <- gsub("drugid_", "", rownames(correlations.stats[["auc"]]))
# yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
yylim <- c(0, 1)
mp <- barplot(height=xx, space=0.3, col=rainbow(length(xx), v=0.9), ylab="Rs", ylim=yylim)
axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
# plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
dev.off()

########################
## correlation across cell lines
########################

myfn <- file.path(saveres, "cgp_ccle_concordance_across_cellines.RData")
if (!file.exists(myfn)) {
  ## correlation for gene expression across cell lines
  # gg <- rownames(annot.ge)
  gg <- rownames(l1000.genes)
  ge.cor <- sapply(gg, function (x, d1, d2) {
    return (cor(d1[ , x], d2[ , x], method="spearman", use="pairwise.complete.obs"))
  }, d1=data.ge.cgp[ , rownames(l1000.genes), drop=FALSE], d2=data.ge.ccle[ , rownames(l1000.genes), drop=FALSE])
  ## correlation for ic50 across cell lines
  dd <- rownames(druginfo)
  ic50.cor <- sapply(dd, function (x, d1, d2) {
    return (cor(d1[ , x], d2[ , x], method="spearman", use="pairwise.complete.obs"))
  }, d1=drugpheno.cgp$IC50, d2=drugpheno.ccle$IC50)
  ## correlation for auc across cell lines
  dd <- rownames(druginfo)
  auc.cor <- sapply(dd, function (x, d1, d2) {
    return (cor(d1[ , x], d2[ , x], method="spearman", use="pairwise.complete.obs"))
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

pdf(file.path(saveres, "cgp_ccle_cor_across_cellines_boxplot.pdf"))
## test significance of the difference between genomic and drug sensitivity data
w1 <- wilcox.test(x=ge.cor, y=auc.cor, conf.int=TRUE)
w2 <- wilcox.test(x=ge.cor, y=ic50.cor, conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E", w1$p.value, w2$p.value)
boxplot(list("GE"=ge.cor, "AUC"=auc.cor, "IC50"=ic50.cor), main="Concordance across cell lines", ylab=expression(R[s]), sub=ss, ylim=yylim)
dev.off()

pdf(file.path(saveres, "cgp_ccle_kappa_across_cellines_boxplot.pdf"))
## test significance of the difference between genomic and drug sensitivity data
w1 <- wilcox.test(x=mut.kappa, y=ic50.call.kappa, conf.int=TRUE)
w2 <- wilcox.test(x=mut.kappa, y=auc.call.kappa, conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("Mutation vs. IC50 calls = %.1E\nMutation vs. AUC calls = %.1E", w1$p.value, w2$p.value)
boxplot(list("Mutations"=mut.kappa, "IC50 calls"=ic50.call.kappa, "AUC calls"=auc.call.kappa), main="Concordance across cell lines", ylab=expression(kappa), sub=ss, ylim=yylim)
dev.off()

########################
## correlation between cell lines
########################

myfn <- file.path(saveres, "cgp_ccle_concordance_between_cellines.RData")
if (!file.exists(myfn)) {
  cellid <- rownames(data.ge.cgp)
  ## correlation for gene expression across cell lines
  ge.cor <- sapply(cellid, function (x, d1, d2) {
    return (cor(d1[x, ], d2[x, ], method="spearman", use="pairwise.complete.obs"))
  }, d1=data.ge.cgp[ , rownames(l1000.genes), drop=FALSE], d2=data.ge.ccle[ , rownames(l1000.genes), drop=FALSE])
  ## correlation for ic50 across cell lines
  ic50.cor <- sapply(cellid, function (x, d1, d2) {
    return (cor(d1[x, ], d2[x, ], method="spearman", use="pairwise.complete.obs"))
  }, d1=drugpheno.cgp$IC50, d2=drugpheno.ccle$IC50)
  ## correlation for auc across cell lines
  auc.cor <- sapply(cellid, function (x, d1, d2) {
    return (cor(d1[x, ], d2[x, ], method="spearman", use="pairwise.complete.obs"))
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

pdf(file.path(saveres, "cgp_ccle_cor_between_cellines_boxplot.pdf"))
## test significance of the difference between genomic and drug sensitivity data
w1 <- wilcox.test(x=ge.cor, y=auc.cor, conf.int=TRUE)
w2 <- wilcox.test(x=ge.cor, y=ic50.cor, conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("GE vs. AUC = %.1E\nGE vs. IC50 = %.1E", w1$p.value, w2$p.value)
boxplot(list("GE"=ge.cor, "AUC"=auc.cor, "IC50"=ic50.cor), main="Concordance between cell lines", ylab=expression(R[s]), sub=ss, ylim=yylim)
dev.off()

pdf(file.path(saveres, "cgp_ccle_kappa_between_cellines_boxplot.pdf"))
## test significance of the difference between genomic and drug sensitivity data
w1 <- wilcox.test(x=mut.kappa, y=ic50.call.kappa, conf.int=TRUE)
w2 <- wilcox.test(x=mut.kappa, y=auc.call.kappa, conf.int=TRUE)
yylim <- c(-1, 1)
ss <- sprintf("Mutation vs. IC50 calls = %.1E\nMutation vs. AUC calls = %.1E", w1$p.value, w2$p.value)
boxplot(list("Mutations"=mut.kappa, "IC50 calls"=ic50.call.kappa, "AUC calls"=auc.call.kappa), main="Concordance between cell lines", ylab=expression(kappa), sub=ss, ylim=yylim)
dev.off()



########################
## all drugs and cell lines

load(file.path(saveres, "cdrug2_cgp_ccle_all.RData"))
gene.common <- intersect(colnames(data.ge.cgp), colnames(data.ge.ccle))
cell.common <- intersect(rownames(data.ge.cgp), rownames(data.ge.ccle))
data.ge.cgp <- data.ge.cgp[ , gene.common, drop=FALSE]
data.ge.ccle <- data.ge.ccle[ , gene.common, drop=FALSE]
drugsn <- rownames(drug.map)

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
# pointLabel(x=drug.cor, y=drug.var.ccle, labels=drugsn, cex=0.7, font=1, srt=30, pos=4, allowSmallOverlap=TRUE, method="SANN", offset=0)
# dfr <- data.frame("cor"=drug.cor, "var"=drug.var.ccle)
# dfr$t <- c(paste("A",1:10,sep=""),paste("B",1:10,sep=""))
# direct.label(xyplot(var ~ cor, data=dfr, group=drugsn, col=blues9[7]))
# wordcloud::
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
cc <- cor.test(drug.mad.cgp, drug.cor, method="spearman", use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.cgp, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.cgp, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.cor, y=drug.mad.cgp, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="Correlation of AUC between CCLE and CGP", ylab="MAD of AUC in CGP")
text(x=drug.cor, y=drug.mad.cgp, labels=drugsn, cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=2, cex=1)
## variability in ccle
nnn <- sum(complete.cases(drug.mad.ccle, drug.cor))
cc <- cor.test(drug.mad.ccle, drug.cor, method="spearman", use="complete.obs", alternative="two.sided")
yylim <- c(floor(min(drug.mad.ccle, na.rm=TRUE) * 1000) / 1000, ceiling(max(drug.mad.ccle, na.rm=TRUE) * 1200) / 1000)
plot(x=drug.cor, y=drug.mad.ccle, xlim=xxlim, ylim=yylim, pch=20, col=blues9[7], xlab="Correlation of AUC between CCLE and CGP", ylab="MAD of AUC in CCLE")
text(x=drug.cor, y=drug.mad.ccle, labels=drugsn, cex=0.7, font=1, srt=30, pos=4)
legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E", cc$estimate, cc$p.value), text.font=2, cex=1)

# pointLabel(x=drug.cor, y=drug.mad.ccle, labels=drugsn, cex=0.7, font=1, srt=30, pos=4, allowSmallOverlap=TRUE, method="SANN", offset=0)
# dfr <- data.frame("cor"=drug.cor, "var"=drug.mad.ccle)
# dfr$t <- c(paste("A",1:10,sep=""),paste("B",1:10,sep=""))
# direct.label(xyplot(var ~ cor, data=dfr, group=drugsn, col=blues9[7]))
# wordcloud::
dev.off()

########################
## gene-drug associations using Spearman correlation
########################

guided.p <- file.path(saveres, "guided_pvalue")
if(!file.exists(guided.p)) { dir.create(guided.p, showWarnings=FALSE, recursive=TRUE) }
guided.var <- file.path(saveres, "guided_variance")
if(!file.exists(guided.var)) { dir.create(guided.var, showWarnings=FALSE, recursive=TRUE) }

  
## filtering by variance
ge.var.cgp <- rank(apply(data.ge.cgp, 2, var, na.rm=TRUE), na.last=NA)
ge.var.ccle <- rank(apply(data.ge.ccle, 2, var, na.rm=TRUE), na.last=NA)

## filtering by p-values
## AUC
myfn <- file.path(saveres, "cgp_ccle_auc_gene_assocs_spearman.RData")
if (!file.exists(myfn)) {
  splitix <- parallel::splitIndices(nx=nrow(drug.map), ncl=nbcore)
  splitix <- splitix[sapply(splitix, length) > 0]
  ## CGP
  ge.assoc.cgp <- parallel::mclapply(splitix, function(x, drug.map, data.ge, drugpheno) {
    dix <- drug.map[x]
    ge.assoc <- t(apply(data.ge, 2, function (x, y) {
      rr <- cor.test(x, y, method="spearman", use="pairwise.complete.obs", exact=FALSE)
      rr <- c(rr$estimate, rr$p.value)
      return (rr)
    }, y=drugpheno[ , dix]))
    colnames(ge.assoc) <- c("rho", "p")
    return (ge.assoc)
  }, drug.map=drug.map[ , "CGP"], data.ge=data.ge.cgp, drugpheno=drugpheno.cgp$AUC)
  names(ge.assoc.cgp) <- rownames(drug.map)
  ## CCLE
  ge.assoc.ccle <- parallel::mclapply(splitix, function(x, drug.map, data.ge, drugpheno) {
    dix <- drug.map[x]
    ge.assoc <- t(apply(data.ge, 2, function (x, y) {
      rr <- cor.test(x, y, method="spearman", use="pairwise.complete.obs", exact=FALSE)
      rr <- c(rr$estimate, rr$p.value)
      return (rr)
    }, y=drugpheno[ , dix]))
    colnames(ge.assoc) <- c("rho", "p")
    return (ge.assoc)
  }, drug.map=drug.map[ , "CCLE"], data.ge=data.ge.ccle, drugpheno=drugpheno.ccle$AUC)
  names(ge.assoc.ccle) <- rownames(drug.map)
  save(list=c("ge.assoc.cgp", "ge.assoc.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## plot for AUC
load(file.path(saveres, "cgp_ccle_auc_gene_assocs_spearman.RData"))
## guided by p-value
for (i in 1:nrow(drug.map)) {
  ## no filtering
  pdf(file.path(guided.p, sprintf("cgp_ccle_qq_plot_guided_pvalue_auc_%s.pdf", rownames(drug.map)[i])), width=10, height=5)
  par(mfrow=c(1, 2))
  ## CCLE guided by CGP
  prop.conc.coef <- NULL
  tt <- table(sign(ge.assoc.ccle[[i]][ , "rho"]), sign(ge.assoc.cgp[[i]][ , "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.ccle[[i]]
  p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
  plot(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="blue", main=sprintf("Q-Q plot\n%s [AUC] in CCLE", rownames(drug.map)[i]), cex=0.8)
  ## CGP p<0.05
  myx <- !is.na(ge.assoc.cgp[[i]][ , "p"]) & ge.assoc.cgp[[i]][ , "p"] < 0.05
  tt <- table(sign(ge.assoc.ccle[[i]][myx, "rho"]), sign(ge.assoc.cgp[[i]][myx, "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.ccle[[i]][myx, , drop=FALSE]
  if (nrow(tt) > 0) {
    p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
    points(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), pch=20, col="red")
  }
  ## CGP p<0.001
  myx <- !is.na(ge.assoc.cgp[[i]][ , "p"]) & ge.assoc.cgp[[i]][ , "p"] < 0.001
  tt <- table(sign(ge.assoc.ccle[[i]][myx, "rho"]), sign(ge.assoc.cgp[[i]][myx, "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.ccle[[i]][myx, , drop=FALSE]
  if (nrow(tt) > 0) {
    p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
    points(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), pch=20, col="green")
  }
  abline(a=0, b=1, col="black")
  # legend("topleft", title="Filter", legend=c("CGP-guided (p < 0.001)", "CGP-guided (p < 0.05)", "None"), col=c("green", "red", "blue"), pch=c(20, 20, 20), bty="n", cex=0.8)
  # legend("bottomright", title="%% of associations\n with concordant direction (sign)", legend=c(sprintf("%i%%", round(prop.conc.coef[1] * 100)), sprintf("%i%%", round(prop.conc.coef[2] * 100)), sprintf("%i%%", round(prop.conc.coef[3] * 100))), col=c("green", "red", "blue"), pch=c(20, 20, 20), bty="n", cex=0.8)
  legend("topleft", title=sprintf("Filter / %% concordant sign(Rs)"), legend=c(sprintf("CGP-guided (p < 0.001) / %i%%", round(prop.conc.coef[3] * 100)), sprintf("CGP-guided (p < 0.05) / %i%%", round(prop.conc.coef[2] * 100)), sprintf("None / %i%%", round(prop.conc.coef[1] * 100))), col=c("green", "red", "blue"), pch=c(20, 20, 20), bty="n", cex=0.8)
  ## CGP guided by CCLE
  prop.conc.coef <- NULL
  tt <- table(sign(ge.assoc.ccle[[i]][ , "rho"]), sign(ge.assoc.cgp[[i]][ , "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.cgp[[i]]
  p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
  plot(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="blue", main=sprintf("Q-Q plot\n%s [AUC] in CGP", rownames(drug.map)[i]), cex=0.8)
  ## CCLE p<0.05
  myx <- !is.na(ge.assoc.ccle[[i]][ , "p"]) & ge.assoc.ccle[[i]][ , "p"] < 0.05
  tt <- table(sign(ge.assoc.ccle[[i]][myx, "rho"]), sign(ge.assoc.cgp[[i]][myx, "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.cgp[[i]][myx, , drop=FALSE]
  if (nrow(tt) > 0) {
    p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
    points(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), pch=20, col="red")
  }
  ## CCLE p<0.001
  myx <- !is.na(ge.assoc.ccle[[i]][ , "p"]) & ge.assoc.ccle[[i]][ , "p"] < 0.001
  tt <- table(sign(ge.assoc.ccle[[i]][myx, "rho"]), sign(ge.assoc.cgp[[i]][myx, "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.cgp[[i]][myx, , drop=FALSE]
  if (nrow(tt) > 0) {
    p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
    points(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), pch=20, col="green")
  }
  abline(a=0, b=1, col="black")
   legend("topleft", title=sprintf("Filter / %% concordant sign(Rs)"), legend=c(sprintf("CCLE-guided (p < 0.001) / %i%%", round(prop.conc.coef[3] * 100)), sprintf("CCLE-guided (p < 0.05) / %i%%", round(prop.conc.coef[2] * 100)), sprintf("None / %i%%", round(prop.conc.coef[1] * 100))), col=c("green", "red", "blue"), pch=c(20, 20, 20), bty="n", cex=0.8)
  dev.off()
}
## guided by variance
for (i in 1:nrow(drug.map)) {
  ## no filtering
  pdf(file.path(guided.var, sprintf("cgp_ccle_qq_plot_guided_variance_auc_%s.pdf", rownames(drug.map)[i])), width=10, height=5)
  par(mfrow=c(1, 2))
  ## CCLE guided by CGP
  prop.conc.coef <- NULL
  tt <- table(sign(ge.assoc.ccle[[i]][ , "rho"]), sign(ge.assoc.cgp[[i]][ , "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.ccle[[i]]
  p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
  plot(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="blue", main=sprintf("Q-Q plot\n%s [AUC] in CCLE", rownames(drug.map)[i]), cex=0.8)
  ## CGP 10% most variant
  myx <- names(ge.var.cgp)[order(ge.var.cgp, decreasing=TRUE)[1:round(length(ge.var.cgp) / 10)]]
  tt <- table(sign(ge.assoc.ccle[[i]][myx, "rho"]), sign(ge.assoc.cgp[[i]][myx, "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.ccle[[i]][myx, , drop=FALSE]
  if (nrow(tt) > 0) {
    p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
    points(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), pch=20, col="red")
  }
  ## CGP p<0.001
  myx <- names(ge.var.cgp)[order(ge.var.cgp, decreasing=TRUE)[1:round(length(ge.var.cgp) / 100)]]
  tt <- table(sign(ge.assoc.ccle[[i]][myx, "rho"]), sign(ge.assoc.cgp[[i]][myx, "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.ccle[[i]][myx, , drop=FALSE]
  if (nrow(tt) > 0) {
    p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
    points(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), pch=20, col="green")
  }
  abline(a=0, b=1, col="black")
  # legend("topleft", title="Filter", legend=c("CGP-guided (p < 0.001)", "CGP-guided (p < 0.05)", "None"), col=c("green", "red", "blue"), pch=c(20, 20, 20), bty="n", cex=0.8)
  # legend("bottomright", title="%% of associations\n with concordant direction (sign)", legend=c(sprintf("%i%%", round(prop.conc.coef[1] * 100)), sprintf("%i%%", round(prop.conc.coef[2] * 100)), sprintf("%i%%", round(prop.conc.coef[3] * 100))), col=c("green", "red", "blue"), pch=c(20, 20, 20), bty="n", cex=0.8)
  legend("topleft", title=sprintf("Filter / %% concordant sign(Rs)"), legend=c(sprintf("CGP-guided (1%% most variant) / %i%%", round(prop.conc.coef[3] * 100)), sprintf("CGP-guided (10%% most variant) / %i%%", round(prop.conc.coef[2] * 100)), sprintf("None / %i%%", round(prop.conc.coef[1] * 100))), col=c("green", "red", "blue"), pch=c(20, 20, 20), bty="n", cex=0.8)
  ## CGP guided by CCLE
  prop.conc.coef <- NULL
  tt <- table(sign(ge.assoc.ccle[[i]][ , "rho"]), sign(ge.assoc.cgp[[i]][ , "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.cgp[[i]]
  p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
  plot(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), xlab="Expected -log10(p)", ylab="Observed -log10(p)", pch=20, col="blue", main=sprintf("Q-Q plot\n%s [AUC] in CGP", rownames(drug.map)[i]), cex=0.8)
  ## CCLE p<0.05
  myx <- names(ge.var.ccle)[order(ge.var.ccle, decreasing=TRUE)[1:round(length(ge.var.ccle) / 10)]]
  tt <- table(sign(ge.assoc.ccle[[i]][myx, "rho"]), sign(ge.assoc.cgp[[i]][myx, "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.cgp[[i]][myx, , drop=FALSE]
  if (nrow(tt) > 0) {
    p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
    points(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), pch=20, col="red")
  }
  ## CCLE p<0.001
  myx <- names(ge.var.ccle)[order(ge.var.ccle, decreasing=TRUE)[1:round(length(ge.var.ccle) / 100)]]
  tt <- table(sign(ge.assoc.ccle[[i]][myx, "rho"]), sign(ge.assoc.cgp[[i]][myx, "rho"]))
  prop.conc.coef <- c(prop.conc.coef, sum(diag(tt)) / sum(tt))
  tt <- ge.assoc.cgp[[i]][myx, , drop=FALSE]
  if (nrow(tt) > 0) {
    p.expected <- seq(0, 1, length.out=nrow(tt)+2)[-c(1, nrow(tt)+2)]
    points(x=-log10(p.expected), y=-log10(sort(tt[ , "p"], decreasing=FALSE)), pch=20, col="green")
  }
  abline(a=0, b=1, col="black")
   legend("topleft", title=sprintf("Filter / %% concordant sign(Rs)"), legend=c(sprintf("CCLE-guided (p < 0.001) / %i%%", round(prop.conc.coef[3] * 100)), sprintf("CCLE-guided (p < 0.05) / %i%%", round(prop.conc.coef[2] * 100)), sprintf("None / %i%%", round(prop.conc.coef[1] * 100))), col=c("green", "red", "blue"), pch=c(20, 20, 20), bty="n", cex=0.8)
  dev.off()
}


########################
## additional analyses
########################

## gene-drug associations from AUC values

## all data
load(file.path(saveres, "cdrug2_cgp_ccle_all.RData"))
tissue.cgp <- as.character(sampleinfo.cgp[ , "tissue.type"])
tissue.ccle <- as.character(sampleinfo.ccle[ , "tissue.type"])
## keep only the common drugs
## CGP
for (i in 1:length(drugpheno.cgp)) {
  drugpheno.cgp[[i]] <- drugpheno.cgp[[i]][ , drug.map[ , "CGP"], drop=FALSE]
  colnames(drugpheno.cgp[[i]]) <- drug.map[ , "CCLE"]
}
auc.cgp <- drugpheno.cgp$AUC
## CCLE
for (i in 1:length(drugpheno.ccle)) {
  drugpheno.ccle[[i]] <- drugpheno.ccle[[i]][ , drug.map[ , "CCLE"], drop=FALSE]
}
auc.ccle <- drugpheno.ccle$AUC

message("Gene-drug association with tuned AUC sensitvity calls, shared data:")

myfn <- file.path(saveres, "cgp_ccle_auc_assoc_all.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC (CCLE)")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      levels(auc)[levels(auc) == "intermediate"] <- NA
      suppressWarnings(res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm))
      return(res)      
    }, data=data.ge.ccle, auc=auc.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ge.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.ccle <- c(assoc.auc.ccle, list(mcres))
  }
  message("")
  names(assoc.auc.ccle) <- colnames(auc.ccle)
  ## CGP
  message("Gene-drug association based on AUC (CGP)")
  assoc.auc.cgp <- NULL
  for(i in 1:ncol(auc.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      levels(auc)[levels(auc) == "intermediate"] <- NA
      suppressWarnings(res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm))
      return(res)      
    }, data=data.ge.cgp, auc=auc.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ge.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.cgp <- c(assoc.auc.cgp, list(mcres))
  }
  message("")
  names(assoc.auc.cgp) <- colnames(auc.cgp)
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.auc.ccle)) {
    tt <- cbind(assoc.auc.ccle[[i]], "EntrezID"=annot.ge.ccle[rownames(assoc.auc.ccle[[i]]), "EntrezID"], "symbol"=annot.ge.ccle[rownames(assoc.auc.ccle[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.ccle)
 #  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ccle_auc_results_gene_drug_paper.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.auc.cgp)) {
    tt <- cbind(assoc.auc.cgp[[i]], "EntrezID"=annot.ge.cgp[rownames(assoc.auc.cgp[[i]]), "EntrezID"], "ymbol"=annot.ge.cgp[rownames(assoc.auc.cgp[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.cgp)
  # WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "cgp_auc_results_gene_drug_paper.xls"), row.names=TRUE)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

gix <- intersect(colnames(data.ge.cgp), colnames(data.ge.ccle))

## sort gene drug assocation by quantile of p-values
auc.assoc.conc.all <- NULL
for(i in 1:ncol(auc.ccle)) {
  xx <- assoc.auc.cgp[[i]][gix, , drop=FALSE]
  yy <- assoc.auc.ccle[[i]][gix, , drop=FALSE]
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
  auc.assoc.conc.all <- c(auc.assoc.conc.all, list(res))
}
names(auc.assoc.conc.all) <- gsub("drugid_", "", colnames(auc.ccle))
## barplot
pdf(file.path(saveres, "auc_assoc_conc_quantiles_all.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.assoc.conc.all)) {
  xx <- auc.assoc.conc.all[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.assoc.conc.all[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.assoc.conc.all[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.assoc.conc.all[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.assoc.conc.all[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.assoc.conc.all)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  # text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Full CGP and CCLE data\n\n\nP quantile cutoff", legend=sprintf("%s%%", gsub("PQUANTILE[.]", " ", rownames(auc.assoc.conc.all[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()

## sort gene drug assocation by fdr cutoffs
auc.assoc.conc.all <- NULL
for(i in 1:ncol(auc.ccle)) {
  xx <- assoc.auc.cgp[[i]][gix, , drop=FALSE]
  yy <- assoc.auc.ccle[[i]][gix, , drop=FALSE]
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
  auc.assoc.conc.all <- c(auc.assoc.conc.all, list(res))
}
names(auc.assoc.conc.all) <- gsub("drugid_", "", colnames(auc.ccle))
## barplot
pdf(file.path(saveres, "auc_assoc_conc_fdrs_all_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.assoc.conc.all)) {
  xx <- auc.assoc.conc.all[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.assoc.conc.all[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.assoc.conc.all[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.assoc.conc.all[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.assoc.conc.all[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.assoc.conc.all)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Full CGP and CCLE data\n\n\nFDR cutoff", legend=sprintf("%s%%", gsub("FDR[.]", " ", rownames(auc.assoc.conc.all[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()



## shared data
load(file.path(saveres, "cdrug2_cgp_ccle_common.RData"))
drugsn <- gsub("drugid_", "", rownames(druginfo))
tissue.cgp <- as.character(sampleinfo.cgp[ , "tissue.type"])
tissue.ccle <- as.character(sampleinfo.ccle[ , "tissue.type"])
tissue <- tissue.cgp
tissuen <- sort(unique(as.character(tissue)))

## gene-drug associations using linear models
message("Gene-drug association with AUC values:")

auc.cgp <- drugpheno.cgp$AUC
auc.ccle <- drugpheno.ccle$AUC

myfn <- file.path(saveres, "cgp_ccle_auc_assoc.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC (CCLE)")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      levels(auc)[levels(auc) == "intermediate"] <- NA
      suppressWarnings(res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm))
      return(res)      
    }, data=data.ge.ccle, auc=auc.ccle[ ,i], tissue=tissue.ccle)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ge.ccle), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.ccle <- c(assoc.auc.ccle, list(mcres))
  }
  message("")
  names(assoc.auc.ccle) <- colnames(auc.ccle)
  ## CGP
  message("Gene-drug association based on AUC (CGP)")
  assoc.auc.cgp <- NULL
  for(i in 1:ncol(auc.cgp)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.cgp)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.cgp), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
      levels(auc)[levels(auc) == "intermediate"] <- NA
      suppressWarnings(res <- apply(X=data[ , x, drop=FALSE], MARGIN=2, FUN=gene.drug.assocs, y=auc, z=tissue, method=genedrugm))
      return(res)      
    }, data=data.ge.cgp, auc=auc.cgp[ ,i], tissue=tissue.cgp)
    mcres <- t(do.call(cbind, mcres))
    mcres <- mcres[colnames(data.ge.cgp), , drop=FALSE]
    mcres <- cbind(mcres, "fdr"=p.adjust(mcres[ ,"pvalue"], method="fdr"))
    assoc.auc.cgp <- c(assoc.auc.cgp, list(mcres))
  }
  message("")
  names(assoc.auc.cgp) <- colnames(auc.cgp)
  ## save all associations
  ## CCLE
  rr <- NULL
  for(i in 1:length(assoc.auc.ccle)) {
    tt <- cbind(assoc.auc.ccle[[i]], "EntrezID"=annot.ge[rownames(assoc.auc.ccle[[i]]), "EntrezID"], "symbol"=annot.ge[rownames(assoc.auc.ccle[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.ccle)
  # WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ccle_auc_results_gene_drug_paper.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.auc.cgp)) {
    tt <- cbind(assoc.auc.cgp[[i]], "EntrezID"=annot.ge[rownames(assoc.auc.cgp[[i]]), "EntrezID"], "ymbol"=annot.ge[rownames(assoc.auc.cgp[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.cgp)
  # WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "cgp_auc_results_gene_drug_paper.xls"), row.names=TRUE)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

## sort gene drug assocation by quantile of p-values
auc.assocs.conc <- NULL
for(i in 1:ncol(auc.ccle)) {
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
  auc.assocs.conc <- c(auc.assocs.conc, list(res))
}
names(auc.assocs.conc) <- gsub("drugid_", "", colnames(auc.ccle))
## barplot
pdf(file.path(saveres, "auc_assocs_conc_quantiles.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.assocs.conc)) {
  xx <- auc.assocs.conc[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.assocs.conc[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.assocs.conc[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.assocs.conc[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.assocs.conc[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.assocs.conc)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  # text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Shared CGP and CCLE data\n\n\nP quantile cutoff", legend=sprintf("%s%%", gsub("PQUANTILE[.]", " ", rownames(auc.assocs.conc[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()

## sort gene drug assocation by fdr cutoffs
auc.assocs.conc <- NULL
for(i in 1:ncol(auc.ccle)) {
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
  auc.assocs.conc <- c(auc.assocs.conc, list(res))
}
names(auc.assocs.conc) <- gsub("drugid_", "", colnames(auc.ccle))
## barplot
pdf(file.path(saveres, "auc_assocs_conc_fdrs_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.assocs.conc)) {
  xx <- auc.assocs.conc[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.assocs.conc[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.assocs.conc[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.assocs.conc[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.assocs.conc[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.assocs.conc)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Shared CGP and CCLE data\n\n\nFDR cutoff", legend=sprintf("%s%%", gsub("FDR[.]", " ", rownames(auc.assocs.conc[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()

## sensitivity calls

## load full data
load(file.path(saveres, "cdrug2_cgp_ccle_all.RData"))
tissue.cgp <- as.character(sampleinfo.cgp[ , "tissue.type"])
tissue.ccle <- as.character(sampleinfo.ccle[ , "tissue.type"])

## tune sensitivity calls to get only extremely sensitive cell lines for targeted drugs
# median absolute deviation
## CGP
myx <- !is.na(druginfo.cgp[, "Drug.class.II"]) & druginfo.cgp[, "Drug.class.II"] == "Cytotoxic"
cytix.cgp <- factor(myx)
levels(cytix.cgp) <- c("Targeted", "Cytotoxic")
mad.drug.cgp <- apply(drugpheno.cgp$AUC, 2, mad, na.rm=TRUE)
## CCLE
myx <- !is.na(druginfo.ccle[, "Class"]) & druginfo.ccle[, "Class"] == "Cytotoxic"
cytix.ccle <- factor(myx)
levels(cytix.ccle) <- c("Targeted", "Cytotoxic")
mad.drug.ccle <- apply(drugpheno.ccle$AUC, 2, mad, na.rm=TRUE)

yylim <- c(0, max(c(mad.drug.cgp, mad.drug.ccle), na.rm=TRUE))
## CGP: sensitivity variability for cytotoxic drugs
pdf(file.path(saveres, "cgp_mad_cytotoxic_drugs.pdf"))
w1 <- wilcox.test(mad.drug.cgp ~ cytix.cgp, conf.int=TRUE)
ss <- sprintf("Wilcoxon rank sum test p = %.1E", w1$p.value)
boxplot(mad.drug.cgp ~ cytix.cgp, ylab="MAD of AUC", main="Distribution of drug sensitivity (AUC) in CGP", ylim=yylim, sub=ss)
abline(h=0.10, col="red", lty=2)
dev.off()
## CCLE: sensitivity variability for cytotoxic drugs
pdf(file.path(saveres, "ccle_mad_cytotoxic_drugs.pdf"))
w1 <- wilcox.test(mad.drug.ccle ~ cytix.ccle, conf.int=TRUE)
ss <- sprintf("Wilcoxon rank sum test p = %.1E", w1$p.value)
boxplot(mad.drug.ccle ~ cytix.ccle, ylab="MAD of AUC", main="Distribution of drug sensitivity (AUC) in CCLE", ylim=yylim, sub=ss)
abline(h=0.10, col="red", lty=2)
dev.off()
## cutoff for MAD of cytotoxic drug, quantile 10%
# mad.cytotoxic <- quantile(mad.drug.cgp[!is.na(cytix.ccle) & cytix.ccle == "Cytotoxic"], probs=0.10, na.rm=TRUE)
mad.cytotoxic <- 0.10

## departure from the pareto distribution
myfn <- file.path(saveres, "cgp_ccle_auc_gpd.RData")
if (!file.exists(myfn)) {
  ## CGP
  myx <- !is.na(druginfo.cgp[, "Drug.class.II"]) & druginfo.cgp[, "Drug.class.II"] == "Cytotoxic"
  cytix.cgp <- factor(myx)
  levels(cytix.cgp) <- c("Targeted", "Cytotoxic")
  splitix <- parallel::splitIndices(nx=ncol(drugpheno.cgp$AUC), ncl=nbcore)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(iix, drugpheno, minsample) {
    res <- apply(drugpheno[ , iix, drop=FALSE], 2, function (x, minsample) {
      x <- x[!is.na(x)]
      rr <- NA
      if (length(x) > minsample) {
        # rr <- gPdtest::gpd.test(x=x, J=1000)$boot.test$p.value
        rr <- pareto.test(x=x, B=1000)$D
      }
      return (rr)
    }, minsample=minsample)
    return (res)
  }, drugpheno=drugpheno.cgp$AUC, minsample=minsample)
  gpd.drug.cgp <- unlist(mcres)
  ## CCLE
  myx <- !is.na(druginfo.ccle[, "Class"]) & druginfo.ccle[, "Class"] == "Cytotoxic"
  cytix.ccle <- factor(myx)
  levels(cytix.ccle) <- c("Targeted", "Cytotoxic")
  splitix <- parallel::splitIndices(nx=ncol(drugpheno.ccle$AUC), ncl=nbcore)
  splitix <- splitix[sapply(splitix, length) > 0]
  mcres <- parallel::mclapply(splitix, function(iix, drugpheno, minsample) {
    res <- apply(drugpheno[ , iix, drop=FALSE], 2, function (x, minsample) {
      x <- x[!is.na(x)]
      rr <- NA
      if (length(x) > minsample) {
        # rr <- gPdtest::gpd.test(x=x, J=1000)$boot.test$p.value
        rr <- pareto.test(x=x, B=1000)$D
      }
      return (rr)
    }, minsample=minsample)
    return (res)
  }, drugpheno=drugpheno.ccle$AUC, minsample=minsample)
  gpd.drug.ccle <- unlist(mcres)
  save(list=c("gpd.drug.ccle", "gpd.drug.cgp"), compress=TRUE, file=myfn)
} else { load(myfn) }

## boxplots
xx <- c(gpd.drug.cgp, gpd.drug.ccle)
xx <- xx[is.finite(xx)]
yylim <- c(0, max(xx, na.rm=TRUE))
## CGP: sensitivity variability for cytotoxic drugs
pdf(file.path(saveres, "cgp_gpd_cytotoxic_drugs.pdf"))
w1 <- wilcox.test(gpd.drug.cgp ~ cytix.cgp, conf.int=TRUE)
ss <- sprintf("Wilcoxon rank sum test p = %.1E", w1$p.value)
boxplot(gpd.drug.cgp ~ cytix.cgp, ylab="KS stat", main="Distribution of drug sensitivity (AUC) in CGP", ylim=yylim, sub=ss)
abline(h=0.10, col="red", lty=2)
dev.off()
## CCLE: sensitivity variability for cytotoxic drugs
pdf(file.path(saveres, "ccle_gpd_cytotoxic_drugs.pdf"))
w1 <- wilcox.test(gpd.drug.ccle ~ cytix.ccle, conf.int=TRUE)
ss <- sprintf("Wilcoxon rank sum test p = %.1E", w1$p.value)
boxplot(gpd.drug.ccle ~ cytix.ccle, ylab="KS stat", main="Distribution of drug sensitivity (AUC) in CCLE", ylim=yylim, sub=ss)
abline(h=0.10, col="red", lty=2)
dev.off()
## CGP+CCLE: sensitivity variability for cytotoxic drugs
pdf(file.path(saveres, "cgp_ccle_gpd_cytotoxic_drugs_paper.pdf"))
dd <- data.frame("gpd"=c(gpd.drug.cgp, gpd.drug.ccle), "cytix"=c(as.character(cytix.cgp), as.character(cytix.ccle)))
dd$cytix <- factor(dd$cytix, levels=c("Targeted", "Cytotoxic"))
dd <- dd[complete.cases(dd), , drop=FALSE]
w1 <- wilcox.test(gpd ~ cytix, data=dd, conf.int=TRUE)
ss <- sprintf("Wilcoxon rank sum test p = %.1E", w1$p.value)
boxplot(gpd ~ cytix, data=dd, ylab="Kolmogorov-Smirnov statistic", main="Distribution of KS statistics (AUC) in CGP+CCLE\nDeviation from the Pareto distribution for cytotoxic drugs", ylim=yylim, sub=ss)
gpd.cutoff <- OptimalCutpoints::optimal.cutpoints(X="gpd", status="cytix", data=dd, methods="Youden", tag.healthy="Targeted")
gpd.cutoff <- gpd.cutoff$Youden$Global$optimal.cutoff$cutoff
abline(h=gpd.cutoff, col="red", lty=2)
dev.off()

## auc sensitivity calling using the new waterfall approach
myfn <- file.path(saveres, "cgp_ccle_auc_call3.RData")
if (!file.exists(myfn)) {
  ## CGP
  auc.call3.cgp <- data.frame(matrix(NA, ncol=ncol(drugpheno.cgp$AUC), nrow=nrow(drugpheno.cgp$AUC), dimnames=dimnames(drugpheno.cgp$AUC)))
  ## AUC sensitivity calling 3 levels
  pdf(file.path(saveres, "cgp_auc_sensitivity_calling3_drugs_tuned.pdf"), width=10, height=10)
  drugn.cgp <- druginfo.cgp[ , "drug.name"]
  drugn.cgp[is.na(druginfo.cgp[ , "drug.name"])] <- rownames(druginfo.cgp)[is.na(druginfo.cgp[ , "drug.name"])]
  for(i in 1:ncol(drugpheno.cgp$AUC)) {
    auc.call3.cgp [ , i] <- callingWaterfallAUC(x=drugpheno.cgp$AUC[ , i], drug.class="Unspecified", drug.ks=gpd.cutoff, intermediate.fold=c("Targeted"=4, "Cytotoxic"=1.2), name=sprintf("%s (CGP)", drugn.cgp[i]), plot=TRUE)
  }
  dev.off()
  dimnames(auc.call3.cgp) <- dimnames(drugpheno.cgp$AUC)
  ## CCLE
  auc.call3.ccle <- data.frame(matrix(NA, ncol=ncol(drugpheno.ccle$AUC), nrow=nrow(drugpheno.ccle$AUC), dimnames=dimnames(drugpheno.ccle$AUC)))
  ## AUC sensitivity calling 3 levels
  pdf(file.path(saveres, "ccle_auc_sensitivity_calling3_drugs_tuned.pdf"), width=10, height=10)
  drugn.ccle <- druginfo.ccle[ , "Compound..code.or.generic.name."]
  drugn.ccle[is.na(druginfo.ccle[ , "Compound..code.or.generic.name."])] <- rownames(druginfo.ccle)[is.na(druginfo.ccle[ , "Compound..code.or.generic.name."])]
  for(i in 1:ncol(drugpheno.ccle$AUC)) {
    auc.call3.ccle [ , i] <- callingWaterfallAUC(x=drugpheno.ccle$AUC[ , i], drug.class="Unspecified", drug.ks=gpd.cutoff, intermediate.fold=c("Targeted"=4, "Cytotoxic"=1.2), name=sprintf("%s (CCLE)", drugn.ccle[i]), plot=TRUE)
  }
  dev.off()
  dimnames(auc.call3.ccle) <- dimnames(drugpheno.ccle$AUC)
  save(list=c("auc.call3.cgp", "auc.call3.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }


########################
## gene-drug associations from AUC extreme sensitivity calls

## all data
load(file.path(saveres, "cdrug2_cgp_ccle_all.RData"))
tissue.cgp <- as.character(sampleinfo.cgp[ , "tissue.type"])
tissue.ccle <- as.character(sampleinfo.ccle[ , "tissue.type"])
## keep only the common drugs
## CGP
for (i in 1:length(drugpheno.cgp)) {
  drugpheno.cgp[[i]] <- drugpheno.cgp[[i]][ , drug.map[ , "CGP"], drop=FALSE]
  colnames(drugpheno.cgp[[i]]) <- drug.map[ , "CCLE"]
}
## CCLE
for (i in 1:length(drugpheno.ccle)) {
  drugpheno.ccle[[i]] <- drugpheno.ccle[[i]][ , drug.map[ , "CCLE"], drop=FALSE]
}
## auc sensitivity call
auc.call3.cgp <- auc.call3.cgp[ , drug.map[ , "CGP"]]
auc.call3.ccle <- auc.call3.ccle[ , drug.map[ , "CCLE"]]
colnames(auc.call3.cgp) <- drug.map[ , "CCLE"]

message("Gene-drug association with tuned AUC sensitvity calls, all data:")

myfn <- file.path(saveres, "cgp_ccle_auc_call3_assoc_all.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC (CCLE)")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.call3.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.call3.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
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
    tt <- cbind(assoc.auc.ccle[[i]], "EntrezID"=annot.ge.ccle[rownames(assoc.auc.ccle[[i]]), "EntrezID"], "symbol"=annot.ge.ccle[rownames(assoc.auc.ccle[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.ccle)
 #  WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "ccle_auc_call3_results_gene_drug_paper.xls"), row.names=TRUE)
  ## CGP
  rr <- NULL
  for(i in 1:length(assoc.auc.cgp)) {
    tt <- cbind(assoc.auc.cgp[[i]], "EntrezID"=annot.ge.cgp[rownames(assoc.auc.cgp[[i]]), "EntrezID"], "ymbol"=annot.ge.cgp[rownames(assoc.auc.cgp[[i]]), "symbol"])
    rr <- c(rr, list(data.frame(tt)))
  }
  names(rr) <- names(assoc.auc.cgp)
  # WriteXLS::WriteXLS("rr", ExcelFileName=file.path(saveres, "cgp_auc_call3_results_gene_drug_paper.xls"), row.names=TRUE)
  save(list=c("assoc.auc.cgp", "assoc.auc.ccle"), compress=TRUE, file=myfn)
} else { load(myfn) }

gix <- intersect(colnames(data.ge.cgp), colnames(data.ge.ccle))

## sort gene drug assocation by quantile of p-values
auc.call3.assoc.conc.all <- NULL
for(i in 1:ncol(auc.call3.ccle)) {
  xx <- assoc.auc.cgp[[i]][gix, , drop=FALSE]
  yy <- assoc.auc.ccle[[i]][gix, , drop=FALSE]
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
  auc.call3.assoc.conc.all <- c(auc.call3.assoc.conc.all, list(res))
}
names(auc.call3.assoc.conc.all) <- gsub("drugid_", "", colnames(auc.call3.ccle))
## barplot
pdf(file.path(saveres, "auc_call3_assoc_conc_quantiles_all.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.call3.assoc.conc.all)) {
  xx <- auc.call3.assoc.conc.all[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.call3.assoc.conc.all[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.call3.assoc.conc.all[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.call3.assoc.conc.all[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.call3.assoc.conc.all[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.call3.assoc.conc.all)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  # text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Full CGP and CCLE data\n\n\nP quantile cutoff", legend=sprintf("%s%%", gsub("PQUANTILE[.]", " ", rownames(auc.call3.assoc.conc.all[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()

## sort gene drug assocation by fdr cutoffs
auc.call3.assoc.conc.all <- NULL
for(i in 1:ncol(auc.call3.ccle)) {
  xx <- assoc.auc.cgp[[i]][gix, , drop=FALSE]
  yy <- assoc.auc.ccle[[i]][gix, , drop=FALSE]
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
  auc.call3.assoc.conc.all <- c(auc.call3.assoc.conc.all, list(res))
}
names(auc.call3.assoc.conc.all) <- gsub("drugid_", "", colnames(auc.call3.ccle))
## barplot
pdf(file.path(saveres, "auc_call3_assoc_conc_fdrs_all_paper.pdf"), height=14, width=14)
par(mfrow=c(4, 4), mar=c(3, 4, 3, 1) + 0.1, xaxt="n", las=1)
for (i in 1:length(auc.call3.assoc.conc.all)) {
  xx <- auc.call3.assoc.conc.all[[i]][ , "rho"]
  xx[!is.na(xx) & xx < 0] <- 0
  xx[!is.na(xx) & xx > 1] <- 1
  ll <- auc.call3.assoc.conc.all[[i]][ , "lower"]
  ll[!is.na(ll) & ll < 0] <- 0
  ll[!is.na(ll) & ll > 1] <- 1
  uu <- auc.call3.assoc.conc.all[[i]][ , "upper"]
  uu[!is.na(uu) & uu < 0] <- 0
  uu[!is.na(uu) & uu > 1] <- 1
  pp <- auc.call3.assoc.conc.all[[i]][ , "p"]
  names(xx) <- names(ll) <- names(uu) <- names(pp) <- rownames(auc.call3.assoc.conc.all[[i]])
  # yylim <- round(range(c(ll, xx, uu), na.rm=TRUE) * 10) / 10
  yylim <- c(0, 1)
  mp <- barplot(height=xx, space=0.3, col=heat.colors(length(xx), alpha=0.9), ylab=concordance.method, ylim=yylim, main=names(auc.call3.assoc.conc.all)[i])
  axis(1, at=seq(1, length(xx), by=1), labels=FALSE)
  # text(x=mp + (max(mp) * 0.0515), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=toupper(names(xx)), srt=45, xpd=NA, cex=0.8, font=2)
  # plotrix::plotCI(x=mp, y=xx, li=ll, ui=uu, err="y", pch=".", add=TRUE)
  text(x=mp + (max(mp) * 0.0515), y=xx, pos=2, labels=ifelse(pp < 0.05, "*", " "), xpd=NA, font=2, cex=1.5)
}
plot.new()
legend("center", bty="n", title="Full CGP and CCLE data\n\n\nFDR cutoff", legend=sprintf("%s%%", gsub("FDR[.]", " ", rownames(auc.call3.assoc.conc.all[[1]]))), col=heat.colors(length(xx), alpha=0.9), pch=15, cex=1.5, pt.cex=3, text.font=1)
dev.off()


########################
message("Correlation of tuned AUC sensitivity calls:")

## shared data
load(file.path(saveres, "cdrug2_cgp_ccle_common.RData"))
drugsn <- gsub("drugid_", "", rownames(druginfo))
tissue.cgp <- as.character(sampleinfo.cgp[ , "tissue.type"])
tissue.ccle <- as.character(sampleinfo.ccle[ , "tissue.type"])
tissue <- tissue.cgp
tissuen <- sort(unique(as.character(tissue)))

## new sensitivity call
auc.call3.cgp <- auc.call3.cgp[rownames(data.ge.cgp), ]
auc.call3.ccle <- auc.call3.ccle[rownames(data.ge.ccle), ]

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
    if(class(err) == "try-error") { rr <- rr.l <- rr.u <- NA }
    rr.l <- rr$CIL
    rr.u <- rr$CIU
    rr <- rr$kappa
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
    if(class(err) == "try-error") { rr2 <- rr2.l <- rr2.u <- NA }
    rr2.l <- rr2$CIL
    rr2.u <- rr2$CIU
    rr2 <- rr2$kappa
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

message("Gene-drug association with tuned AUC sensitvity calls, shared data:")

myfn <- file.path(saveres, "cgp_ccle_auc_call3_assoc.RData")
if(!file.exists(myfn)) {
  ## CCLE
  message("Gene-drug association based on AUC (CCLE)")
  assoc.auc.ccle <- NULL
  for(i in 1:ncol(auc.call3.ccle)) {
    message("Computation for drug ", gsub("drugid_", "", colnames(auc.call3.ccle)[i]))
    splitix <- parallel::splitIndices(nx=ncol(data.ge.ccle), ncl=nbcore)
    mcres <- parallel::mclapply(splitix, function(x, data, auc, tissue) {
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
pdf(file.path(saveres, "auc_call3_assocs_conc_quantiles.pdf"), height=14, width=14)
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
pdf(file.path(saveres, "auc_call3_assocs_conc_fdrs_paper.pdf"), height=14, width=14)
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


## end


