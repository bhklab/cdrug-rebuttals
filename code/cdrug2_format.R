########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## April 8, 2014
########################

# rm(list=ls())

require(vcd) || stop("Library vcd is not available")
require(epibasix) || stop("Library epibasix is not available")
require(R.utils) || stop("Library R.utils is not available")
require(amap) || stop("Library amap is not available")
require(gplots) || stop("Library gplots is not available")
require(VennDiagram) || stop("Library VennDiagram is not available")
require(PharmacoGx) || stop("Library PharmacoGx is not available")

myfn <- file.path(saveres, "cdrug2_cgp_ccle_all.RData")
if (!file.exists(myfn)) {
  ## load processed data
  ## CGP
  load(file.path(saveres, "cgp_data.RData"))
  ## CCLE
  load(file.path(saveres, "ccle_data.RData"))

  ## drug ids in common
  ## manual mapping CCLE vs CGP
  drug.map <- read.csv(file.path("data", "matching_drug_CCLE_CGP"), row.names=1)
  drug.map[ , "CCLE"] <- paste("drugid", rownames(drug.map), sep="_")
  drug.map <- drug.map[complete.cases(drug.map[ , c("CCLE", "CGP")]), c("CCLE", "CGP"), drop=FALSE]
  
  ## cell line annotations
  ## CGP cell line collection
  # celline.cgp <- read.csv(file=file.path(saveres, "cell_line_collection_cgp.csv"))
  # rownames(celline.cgp) <- as.character(celline.cgp[ , "cellid"])
  # celline.cgp[ , "cellid"] <- as.character(celline.cgp[ , "cellid"])
  if(any(!is.element(rownames(data.ge.cgp), rownames(celline.cgp)))) { stop("some cell lines in CGP are not part of CGP cell line collection") }
  ## CCLE cell line collection
  # celline.ccle <- read.csv(file=file.path(saveres, "cell_line_collection_ccle.csv"))
  # rownames(celline.ccle) <- as.character(celline.ccle[ , "cellid"])
  # celline.ccle[ , "cellid"] <- as.character(celline.ccle[ , "cellid"])
  if(any(!is.element(rownames(data.ge.ccle), rownames(celline.ccle)))) { stop("some cell lines in CCLE are not part of CCLE cell line collection") }

  ## read manual matching for cell lines in CCLE and CGP
  match.cellines <- read.csv(file=file.path("code", "cell_line_annotation_all.csv"), row.names=1)

  ## CCLE
  ## update ccle cell line collection
  nn <- intersect(as.character(match.cellines[ , "CCLE"]), as.character(celline.ccle[ , "cellid"]))
  iix0 <- which(is.element(as.character(match.cellines[ , "CCLE.cellid"]), nn))
  iix <- match(as.character(match.cellines[iix0, "CCLE.cellid"]), as.character(celline.ccle[ , "cellid"]))
  celline.ccle[iix, "cellid"] <- rownames(match.cellines)[iix0]
  celline.ccle <- data.frame(celline.ccle, "tissueid"=celline.ccle[ , "Site.Primary"])
  # celline.ccle[iix, "tissueid"] <- as.character(match.cellines[iix0, "unique.tissueid"])
  rownames(celline.ccle) <- as.character(celline.ccle[ , "cellid"])
  ## update ccle gene expression
  nn <- intersect(as.character(match.cellines[ , "CCLE.cellid"]), rownames(data.ge.ccle))
  iix0 <- which(is.element(as.character(match.cellines[ , "CCLE.cellid"]), nn))
  iix <- match(as.character(match.cellines[iix0, "CCLE.cellid"]), rownames(data.ge.ccle))
  rownames(data.ge.ccle)[iix] <- rownames(match.cellines)[iix0]
  ## update ccle sample information
  nn <- intersect(as.character(match.cellines[ , "CCLE.cellid"]), rownames(sampleinfo.ccle))
  iix0 <- which(is.element(as.character(match.cellines[ , "CCLE.cellid"]), nn))
  iix <- match(as.character(match.cellines[iix0, "CCLE.cellid"]), rownames(sampleinfo.ccle))
  rownames(sampleinfo.ccle)[iix] <- rownames(match.cellines)[iix0]
  sampleinfo.ccle[ , "cellid"] <- rownames(sampleinfo.ccle)
  ## update ccle drug phenotypes
  for (i in 1:length(drugpheno.ccle)) {
    nn <- intersect(as.character(match.cellines[ , "CCLE.cellid"]), rownames(drugpheno.ccle[[i]]))
    iix0 <- which(is.element(as.character(match.cellines[ , "CCLE.cellid"]), nn))
    iix <- match(as.character(match.cellines[iix0, "CCLE.cellid"]), rownames(drugpheno.ccle[[i]]))
    rownames(drugpheno.ccle[[i]])[iix] <- rownames(match.cellines)[iix0]
  }
  ## update ccle drug concentrations
  nn <- intersect(as.character(match.cellines[ , "CCLE.cellid"]), drugconc.ccle[ , "cellid"])
  iix0 <- which(is.element(as.character(match.cellines[ , "CCLE.cellid"]), nn))
  for(i in 1:length(iix0)) {
    myx <- !is.na(drugconc.ccle[ , "cellid"]) & drugconc.ccle[ , "cellid"] == match.cellines[iix0[i], "CCLE.cellid"]
    drugconc.ccle[myx, "cellid"] <- rownames(match.cellines)[iix0[i]]
  }
  rownames(drugconc.ccle) <- paste(as.character(drugconc.ccle[ , "drugid"]), as.character(drugconc.ccle[ , "cellid"]), sep="...")
  ## update ccle mutations
  nn <- intersect(as.character(match.cellines[ , "CCLE.cellid"]), rownames(mutation.ccle))
  iix0 <- which(is.element(as.character(match.cellines[ , "CCLE.cellid"]), nn))
  iix <- match(as.character(match.cellines[iix0, "CCLE.cellid"]), rownames(mutation.ccle))
  rownames(mutation.ccle)[iix] <- rownames(match.cellines)[iix0]
  
  ## CGP
  ## update cgp cell line collection
  nn <- intersect(as.character(match.cellines[ , "CGP"]), as.character(celline.cgp[ , "cellid"]))
  iix0 <- which(is.element(as.character(match.cellines[ , "CGP"]), nn))
  iix <- match(as.character(match.cellines[iix0, "CGP"]), as.character(celline.cgp[ , "cellid"]))
  celline.cgp[iix, "cellid"] <- rownames(match.cellines)[iix0]
  celline.cgp <- data.frame(celline.cgp, "tissueid"=celline.cgp[ , "Primary.site"])
  # celline.cgp[iix, "tissueid"] <- as.character(match.cellines[iix0, "unique.tissueid"])
  ## remove possible duplicates
  celline.cgp <- celline.cgp[!duplicated(as.character(celline.cgp[ , "cellid"])), , drop=FALSE]
  rownames(celline.cgp) <- as.character(celline.cgp[ , "cellid"])
  ## update cgp gene expression
  nn <- intersect(as.character(match.cellines[ , "CGP"]), rownames(data.ge.cgp))
  iix0 <- which(is.element(as.character(match.cellines[ , "CGP"]), nn))
  iix <- match(as.character(match.cellines[iix0, "CGP"]), rownames(data.ge.cgp))
  rownames(data.ge.cgp)[iix] <- rownames(match.cellines)[iix0]
  ## update cgp sample information
  nn <- intersect(as.character(match.cellines[ , "CGP"]), rownames(sampleinfo.cgp))
  iix0 <- which(is.element(as.character(match.cellines[ , "CGP"]), nn))
  iix <- match(as.character(match.cellines[iix0, "CGP"]), rownames(sampleinfo.cgp))
  rownames(sampleinfo.cgp)[iix] <- rownames(match.cellines)[iix0]
  sampleinfo.cgp[ , "cellid"] <- rownames(sampleinfo.cgp)
  ## update cgp drug phenotypes
  for (i in 1:length(drugpheno.cgp)) {
    nn <- intersect(as.character(match.cellines[ , "CGP"]), rownames(drugpheno.cgp[[i]]))
    iix0 <- which(is.element(as.character(match.cellines[ , "CGP"]), nn))
    iix <- match(as.character(match.cellines[iix0, "CGP"]), rownames(drugpheno.cgp[[i]]))
    rownames(drugpheno.cgp[[i]])[iix] <- rownames(match.cellines)[iix0]
  }
  ## update cgp drug concentrations
  nn <- intersect(as.character(match.cellines[ , "CGP"]), drugconc.cgp[ , "cellid"])
  iix0 <- which(is.element(as.character(match.cellines[ , "CGP"]), nn))
  for(i in 1:length(iix0)) {
    myx <- !is.na(drugconc.cgp[ , "cellid"]) & drugconc.cgp[ , "cellid"] == match.cellines[iix0[i], "CGP"]
    drugconc.cgp[myx, "cellid"] <- rownames(match.cellines)[iix0[i]]
  }
  rownames(drugconc.cgp) <- paste(as.character(drugconc.cgp[ , "drugid"]), as.character(drugconc.cgp[ , "cellid"]), sep="...")
  ## update cgp mutations
  nn <- intersect(as.character(match.cellines[ , "CGP"]), rownames(mutation.cgp))
  iix0 <- which(is.element(as.character(match.cellines[ , "CGP"]), nn))
  iix <- match(as.character(match.cellines[iix0, "CGP"]), rownames(mutation.cgp))
  rownames(mutation.cgp)[iix] <- rownames(match.cellines)[iix0]
  
  ## intersection between CGP and CCLE
  ll <- list("CGP"=rownames(sampleinfo.cgp), "CCLE"=rownames(sampleinfo.ccle))
  rr <- sapply(overlap(ll), length)
  pdf(file.path(file.path(saveres, "intersection_cellines_cgp_ccle_paper.pdf")), width=5, height=5)
  ww <- VennDiagram::draw.pairwise.venn(area1=rr[1], area2=rr[2], cross.area=rr[3], category=names(ll), fill=c("blue3", "yellow3"), lty="blank", cex=rep(2, length(rr)), cat.cex=rep(2, length(ll)), cat.col=c("blue3", "yellow3"), scaled=FALSE, aplha=rep(0.5, length(ll)), cat.dist=rep(0.1, length(ll)), ind=FALSE, mar=rep(0.1, length(ll)))
  grid::grid.draw(ww)
  dev.off()

  ## merge cell line annotations from CCLE and CGP
  celln <- sort(unique(c(rownames(sampleinfo.ccle), rownames(sampleinfo.cgp))))
  celline.collection <- data.frame(matrix(NA, nrow=length(celln), ncol=7, dimnames=list(celln, c("cellid", "CGP.cell.line", "CCLE.cell.line", "CGP.tissue.type", "CCLE.tissue.type", "CGP.link", "CCLE.link"))))
  celline.collection[ , "cellid"] <- celln
  celline.collection[rownames(data.ge.cgp), c("CGP.cell.line", "CGP.tissue.type", "CGP.link")] <- celline.cgp[rownames(data.ge.cgp), c("Sample.name", "tissueid", "link")]
  celline.collection[rownames(data.ge.ccle), c("CCLE.cell.line", "CCLE.tissue.type", "CCLE.link")] <- celline.ccle[rownames(data.ge.ccle), c("Cell.line.primary.name", "tissueid", "link")]
  tissue.cgp <- celline.cgp[rownames(sampleinfo.cgp), "tissueid"]
  tissue.ccle <- celline.ccle[rownames(sampleinfo.ccle), "tissueid"]
  ## save cell line collection
  tt <- celline.collection[ , "CGP.tissue.type"]
  tt[is.na(tt)] <- celline.collection[is.na(tt), "CCLE.tissue.type"]
  celline.collection <- data.frame(celline.collection, "tissue.type"=tt)
  write.csv(celline.collection, file=file.path(saveres, "cell_line_collection_all.csv"), row.names=FALSE)
  ## update sample information with tissue types
  sampleinfo.cgp <- data.frame(sampleinfo.cgp, "tissue.type"=tissue.cgp)
  sampleinfo.ccle <- data.frame(sampleinfo.ccle, "tissue.type"=tissue.ccle)

  ## mutation data
  message("\tFormat (missense) mutation data")
  write.csv(mutation.cgp, file=file.path(saveres, "mutation_cgp_common.csv"))
  write.csv(mutation.ccle, file=file.path(saveres, "mutation_ccle_common.csv"))

  ## drug sensitivity measures for CGP
  message("\tFormat drug sensitivity measures")

  ## IC50 in micro molar
  message("\tIC50 (CGP)")

  drugn.cgp <- druginfo.cgp[ , "drug.name"]
  drugn.cgp[is.na(druginfo.cgp[ , "drug.name"])] <- gsub("drugid_", "", rownames(druginfo.cgp)[is.na(druginfo.cgp[ , "drug.name"])])
  ## boxplot
  xx <- -log10(drugpheno.cgp$IC50 / 10^6)
  oo <- order(apply(xx, 2, median, na.rm=TRUE), apply(xx, 2, IQR, na.rm=TRUE), decreasing=FALSE)
  xx <- xx[, oo, drop=FALSE]
  colnames(xx) <- drugn.cgp
  ## all drugs
  pdf(file.path(saveres, "boxplot_ic50_cgp.pdf"), width=20, height=10)
  par(las=3, mar=c(10, 4, 4, 2) + 0.1, cex=0.8)
  graphics::boxplot(xx, outline=FALSE, ylab="-log10(IC50)", main="Drug sensitivity\nCGP")
  dev.off()
  ## drugs in common between CGP and CCLE
  commoncell <- intersect(rownames(data.ge.cgp), rownames(data.ge.ccle))
  commonix <- is.element(colnames(xx), drug.map[ , "CGP"])
  mycol <- rep("white", ncol(xx))
  mycol[commonix] <- "red"
  pdf(file.path(saveres, "boxplot_ic50_cgp_commondrugs.pdf"), width=17, height=7)
  par(las=3, mar=c(5, 4, 4, 2) + 0.1, xaxt="n")
  mp <- graphics::boxplot(xx[commoncell, , drop=FALSE], outline=FALSE, ylab="-log10(IC50)", main="Drugs sensitivity (IC50)\nCGP", col=mycol, cex=0.5)
  axis(1, at=which(commonix), tick=TRUE, labels=T)
  text(x=which(commonix) + 0.5, y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=rownames(drug.map), srt=45, xpd=NA, font=2, col="red", cex=0.75)
  # legend("topleft", legend=c("CGP", "CCLE"), col=c("lightgreen", "lightblue"), pch=15, bty="n")
  dev.off()

  ## filtering based on concentration range
  ## ic50 larger or equal to the maximum tested drug concentration are filtered out
  ic50.filt.cgp <- matrix(FALSE, nrow=nrow(drugpheno.cgp$IC50), ncol=ncol(drugpheno.cgp$IC50), dimnames=dimnames(drugpheno.cgp$IC50))
  drugc <- drugconc.cgp[is.element(drugconc.cgp[ , "cellid"], rownames(drugpheno.cgp$IC50)) & is.element(drugconc.cgp[ , "drugid"], colnames(drugpheno.cgp$IC50)), , drop=FALSE]
  maxdose <- drugc[ , "max.Dose.uM"]
  ffilt <- matrix(NA, nrow=nrow(drugpheno.cgp$IC50), ncol=ncol(drugpheno.cgp$IC50), dimnames=dimnames(drugpheno.cgp$IC50))
  ffilt[as.matrix(drugc[ , c("cellid", "drugid")])] <- maxdose
  ic50.filt.cgp[!is.na(drugpheno.cgp$IC50) & drugpheno.cgp$IC50 < ffilt] <- TRUE
  dimnames(ic50.filt.cgp) <- dimnames(drugpheno.cgp$IC50)
  ## ic50.filt is set to TRUE if the IC50 measurement passes the filter

  ## activity area
  message("\tAUC (CGP)")

  ## boxplot for all drugs
  pdf(file.path(saveres, "boxplot_auc_cgp.pdf"), width=20, height=10)
  par(las=3, mar=c(10, 4, 4, 2) + 0.1, cex=0.8)
  oo <- order(apply(drugpheno.cgp$AUC, 2, median, na.rm=TRUE), decreasing=FALSE)
  graphics::boxplot(drugpheno.cgp$AUC[ , oo, drop=FALSE], outline=FALSE, ylab="AUC", main="Drug sensitivity (AUC)\nCGP")
  dev.off()
  ## drugs in common between CGP and CCLE
  commoncell <- intersect(rownames(drugpheno.cgp$AUC), rownames(drugpheno.ccle$AUC))
  commonix <- is.element(colnames(drugpheno.cgp$AUC[ , oo, drop=FALSE]), drug.map[ , "CGP"])
  mycol <- rep("white", ncol(drugpheno.cgp$AUC))
  mycol[commonix] <- "red"
  pdf(file.path(saveres, "boxplot_auc_cgp_commondrugs.pdf"), width=17, height=7)
  par(las=3, mar=c(5, 4, 4, 2) + 0.1, xaxt="n")
  mp <- graphics::boxplot(drugpheno.cgp$AUC[commoncell, oo, drop=FALSE], outline=FALSE, ylab="AUC", main="Drugs sensitivity (AUC)\nCGP", col=mycol, cex=0.5)
  axis(1, at=which(commonix), tick=TRUE, labels=T)
  text(x=which(commonix) + 0.5, y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=rownames(drug.map), srt=45, xpd=NA, font=2, col="red", cex=0.75)
  # legend("topleft", legend=c("CGP", "CCLE"), col=c("lightgreen", "lightblue"), pch=15, bty="n")
  dev.off()


  ## IC50 in micro molar
  message("\tIC50 (CCLE)")

  drugn.ccle <- druginfo.ccle[ , "Compound..code.or.generic.name."]
  drugn.ccle[is.na(druginfo.ccle[ , "Compound..code.or.generic.name."])] <- gsub("drugid_", "", rownames(druginfo.cgp)[is.na(druginfo.cgp[ , "drug.name"])])
  ## boxplot
  xx <- -log10(drugpheno.ccle$IC50 / 10^6)
  oo <- order(apply(xx, 2, median, na.rm=TRUE), apply(xx, 2, IQR, na.rm=TRUE), decreasing=FALSE)
  xx <- xx[, oo, drop=FALSE]
  colnames(xx) <- drugn.ccle
  ## all drugs
  pdf(file.path(saveres, "boxplot_ic50_ccle.pdf"), width=20, height=10)
  par(las=3, mar=c(10, 4, 4, 2) + 0.1, cex=0.8)
  graphics::boxplot(xx, outline=FALSE, ylab="-log10(IC50)", main="Drug sensitivity\nCCLE")
  dev.off()
  ## drugs in common between CCLE and CCLE
  commoncell <- intersect(rownames(drugpheno.ccle$IC50), rownames(drugpheno.cgp$IC50))
  commonix <- is.element(colnames(xx), drug.map[ , "CCLE"])
  mycol <- rep("white", ncol(xx))
  mycol[commonix] <- "red"
  pdf(file.path(saveres, "boxplot_ic50_ccle_commondrugs.pdf"), width=17, height=7)
  par(las=3, mar=c(5, 4, 4, 2) + 0.1, xaxt="n")
  mp <- graphics::boxplot(xx[commoncell, , drop=FALSE], outline=FALSE, ylab="-log10(IC50)", main="Drugs sensitivity (IC50)\nCCLE", col=mycol, cex=0.5)
  axis(1, at=which(commonix), tick=TRUE, labels=T)
  text(x=which(commonix) + 0.5, y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=rownames(drug.map), srt=45, xpd=NA, font=2, col="red", cex=0.75)
  # legend("topleft", legend=c("CCLE", "CCLE"), col=c("lightgreen", "lightblue"), pch=15, bty="n")
  dev.off()

  ## filtering based on concentration range
  ## ic50 larger or equal to the maximum tested drug concentration are filtered out
  ic50.filt.ccle <- matrix(FALSE, nrow=nrow(drugpheno.ccle$IC50), ncol=ncol(drugpheno.ccle$IC50), dimnames=dimnames(drugpheno.ccle$IC50))
  drugc <- drugconc.ccle[is.element(drugconc.ccle[ , "cellid"], rownames(drugpheno.ccle$IC50)) & is.element(drugconc.ccle[ , "drugid"], colnames(drugpheno.ccle$IC50)), , drop=FALSE]
  maxdose <- apply(drugc[ , grep("^Dose", colnames(drugconc.ccle))], 1, function (x) { 
    rr <- NA
    if(!all(is.na(x))) { rr <- max(x, na.rm=TRUE) }
    return (rr)
  })
  ffilt <- matrix(NA, nrow=nrow(drugpheno.ccle$IC50), ncol=ncol(drugpheno.ccle$IC50), dimnames=dimnames(drugpheno.ccle$IC50))
  ffilt[as.matrix(drugc[ , c("cellid", "drugid")])] <- maxdose
  ic50.filt.ccle[!is.na(drugpheno.ccle$IC50) & drugpheno.ccle$IC50 < ffilt] <- TRUE
  dimnames(ic50.filt.ccle) <- dimnames(drugpheno.ccle$IC50)
  ## ic50.filt is set to TRUE if the IC50 measurement passes the filter

  ## activity area
  message("\tAUC (CCLE)")

  ## boxplot for all drugs
  pdf(file.path(saveres, "boxplot_auc_ccle.pdf"), width=20, height=10)
  par(las=3, mar=c(10, 4, 4, 2) + 0.1, cex=0.8)
  oo <- order(apply(drugpheno.ccle$AUC, 2, median, na.rm=TRUE), decreasing=FALSE)
  graphics::boxplot(drugpheno.ccle$AUC[ , oo, drop=FALSE], outline=FALSE, ylab="AUC", main="Drug sensitivity (AUC)\nCCLE")
  dev.off()
  ## drugs in common between CCLE and CCLE
  commoncell <- intersect(rownames(drugpheno.ccle$AUC), rownames(drugpheno.ccle$AUC))
  commonix <- is.element(colnames(drugpheno.ccle$AUC[ , oo, drop=FALSE]), drug.map[ , "CCLE"])
  mycol <- rep("white", ncol(drugpheno.ccle$AUC))
  mycol[commonix] <- "red"
  pdf(file.path(saveres, "boxplot_auc_ccle_commondrugs.pdf"), width=17, height=7)
  par(las=3, mar=c(5, 4, 4, 2) + 0.1, xaxt="n")
  mp <- graphics::boxplot(drugpheno.ccle$AUC[commoncell, oo, drop=FALSE], outline=FALSE, ylab="AUC", main="Drugs sensitivity (AUC)\nCCLE", col=mycol, cex=0.5)
  axis(1, at=which(commonix), tick=TRUE, labels=T)
  text(x=which(commonix) + 0.5, y=par("usr")[3] - (par("usr")[4] * 0.02), pos=2, labels=rownames(drug.map), srt=45, xpd=NA, font=2, col="red", cex=0.75)
  # legend("topleft", legend=c("CCLE", "CCLE"), col=c("lightgreen", "lightblue"), pch=15, bty="n")
  dev.off()

  ## save data
  save(list=c("data.ge.ccle", "data.ge.cgp", "mutation.cgp", "mutation.ccle", "annot.ge.cgp", "annot.ge.ccle", "druginfo.ccle", "druginfo.cgp", "drugpheno.ccle", "drugpheno.cgp", "sampleinfo.cgp", "sampleinfo.ccle", "celline.collection", "drug.map"), compress=TRUE, file=file.path(saveres, "cdrug2_cgp_ccle_full_all.RData"))
  
  ## gene centric data
  ## CGP
  myx <- which(annot.ge.cgp[ ,"best"])
  myx <- myx[!duplicated(annot.ge.cgp[myx,"EntrezID"])]
  data.ge.cgp <- data.ge.cgp[ ,myx,drop=FALSE]
  annot.ge.cgp <- annot.ge.cgp[myx, ,drop=FALSE]
  colnames(data.ge.cgp) <- rownames(annot.ge.cgp) <- paste("geneid", annot.ge.cgp[ ,"EntrezID"], sep="_")
  ## CCLE
  myx <- which(annot.ge.ccle[ ,"best"])
  myx <- myx[!duplicated(annot.ge.ccle[myx,"EntrezID"])]
  data.ge.ccle <- data.ge.ccle[ ,myx,drop=FALSE]
  annot.ge.ccle <- annot.ge.ccle[myx, ,drop=FALSE]
  colnames(data.ge.ccle) <- rownames(annot.ge.ccle) <- paste("geneid", annot.ge.ccle[ ,"EntrezID"], sep="_")
  ## save data
  save(list=c("data.ge.ccle", "data.ge.cgp", "mutation.cgp", "mutation.ccle", "annot.ge.cgp", "annot.ge.ccle", "druginfo.ccle", "druginfo.cgp", "drugpheno.ccle", "drugpheno.cgp", "drugconc.cgp", "drugconc.ccle", "sampleinfo.cgp", "sampleinfo.ccle", "celline.collection", "drug.map"), compress=TRUE, file=myfn)
  
} else { load(myfn) }

########################
## comparison between CCLE and CGP
## cell line id
cellid.common <- fold(intersect, rownames(data.ge.cgp), rownames(data.ge.ccle), rownames(mutation.cgp), rownames(mutation.ccle))
ge.common <- intersect(colnames(data.ge.ccle), colnames(data.ge.cgp))
mut.common <- intersect(colnames(mutation.ccle), colnames(mutation.cgp))
## drug ids; use CCLE drug names for consistency
iix <- drug.map[complete.cases(drug.map[ , c("CGP", "CCLE")]), "CGP"]
names(iix) <- drug.map[complete.cases(drug.map[ , c("CGP", "CCLE")]), "CCLE"]
## drug pheno
drugpheno.common <- intersect(names(drugpheno.cgp), names(drugpheno.ccle))
## CGP
## druginfo
druginfo.cgp <- druginfo.cgp[iix, , drop=FALSE]
rownames(druginfo.cgp) <- names(iix)
## drugpheno
for (i in 1:length(drugpheno.cgp)) {
  drugpheno.cgp[[i]] <- drugpheno.cgp[[i]][cellid.common, iix, drop=FALSE]
  colnames(drugpheno.cgp[[i]]) <- names(iix)
}
drugpheno.cgp <- drugpheno.cgp[drugpheno.common]
## drug concentrations
drugconc.cgp <- drugconc.cgp[is.element(drugconc.cgp[ ,"drugid"], iix), , drop=FALSE]
drugconc.cgp[ , "drugid"] <- names(iix)[sapply(1:nrow(drugconc.cgp), function(x, y, z) { return(which(y[x] == z)) }, y=drugconc.cgp[ , "drugid"], z=iix)]
rownames(drugconc.cgp) <- paste(drugconc.cgp[ , "cellid"], drugconc.cgp[ , "drugid"], sep=".")
## sample information
sampleinfo.cgp <- sampleinfo.cgp[cellid.common, , drop=FALSE]
## gene expression
data.ge.cgp <- data.ge.cgp[cellid.common, ge.common, drop=FALSE]
annot.ge.cgp <- annot.ge.cgp[ge.common, , drop=FALSE]
## mutation
mutation.cgp <- mutation.cgp[cellid.common, mut.common, drop=FALSE]
## CCLE
## druginfo
druginfo.ccle <- druginfo.ccle[names(iix), , drop=FALSE]
## drugpheno
for (i in 1:length(drugpheno.ccle)) {
  drugpheno.ccle[[i]] <- drugpheno.ccle[[i]][cellid.common, names(iix), drop=FALSE]
}
drugpheno.ccle <- drugpheno.ccle[drugpheno.common]
## drug concentrations
drugconc.ccle <- drugconc.ccle[is.element(drugconc.ccle[ ,"drugid"], names(iix)), , drop=FALSE]
## sample information
sampleinfo.ccle <- sampleinfo.ccle[cellid.common, , drop=FALSE]
## gene expression
data.ge.ccle <- data.ge.ccle[cellid.common, ge.common, drop=FALSE]
annot.ge.ccle <- annot.ge.ccle[ge.common, , drop=FALSE]
## mutation
mutation.ccle <- mutation.ccle[cellid.common, mut.common, drop=FALSE]

## annotations for gene expression
annot.ge <- annot.ge.ccle
## annotations for cell lines
celline.collection <- celline.collection[cellid.common, , drop=FALSE]
## drug information
tt1 <- druginfo.cgp
colnames(tt1) <- paste("CGP", colnames(tt1), sep=".")
tt2 <- druginfo.ccle
colnames(tt2) <- paste("CCLE", colnames(tt2), sep=".")
druginfo <- cbind(tt1, tt2)

save(list=c("data.ge.ccle", "data.ge.cgp", "mutation.cgp", "mutation.ccle", "annot.ge", "druginfo.ccle", "druginfo.cgp", "drugpheno.ccle", "drugpheno.cgp", "sampleinfo.cgp", "sampleinfo.ccle", "celline.collection", "druginfo"), compress=TRUE, file=file.path(saveres, "cdrug2_cgp_ccle_common.RData"))



## end


