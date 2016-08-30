### Scatterplot with transparency
myScatterPlot <- function(x, y, method=c("plain", "transparent", "smooth"), transparency=0.10, smooth.pch=".", pch=16, minp=50, col=blues9[7], smooth.col=c("white", blues9), ...) {
  require(grDevices) || stop("Library grDevices is not available!")
  method <- match.arg(method)
  if (length(col) != length(x)) {
    col <- rep(col, length.out=length(x))
  }
  ccix <- complete.cases(x, y)
  x <- x[ccix]
  y <- y[ccix]
  col <- col[ccix]
  
  if (sum(ccix) < minp) {
    ## too few points, no transparency, no smoothing
    if (sum(ccix) > 0) { rr <- plot(x=x, y=y, col=col, pch=pch, ...) } else { rr <- plot(x=x, y=y, col=col, pch=pch, ...) }
  } else {
    ## enough data points
    switch(method,
           "plain"={
             rr <- plot(x=x, y=y, col=col, pch=pch, ...)
           },
           "transparent"={
             myrgb <- sapply(col, grDevices::col2rgb, alpha=FALSE) / 255
             myrgb <- apply(myrgb, 2, function (x, transparency) {
               return (rgb(red=x[1], green=x[2], blue=x[3], alpha=transparency, maxColorValue=1))
             }, transparency=transparency)
             rr <- plot(x=x, y=y, pch=pch, col=myrgb, ...)
           },
           "smooth"={
             rr <- smoothScatter(x=x, y=y, col="lightgray", colramp=colorRampPalette(smooth.col), pch=smooth.pch, ...)
           }
    )
  }
  
  invisible(rr)
}

# ### check known biomarkers
# knownBiomarkersCheck <-
#   function(ccle.sig.rna, gdsc.sig.rna, ccle.sig.mutation, gdsc.sig.mutation, gdsc.sig.fusion, ccle.sig.cnv, gdsc.sig.cnv, method, cell)
#   {
#
#     known.biomarkers <- read.csv(file="known_biomarkers.csv", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, na.strings=c("", " ", "NA"))
#     known.biomarkers <- cbind(known.biomarkers, "GDSC effect size"=NA, "GDSC pvalue"=NA, "GDSC FDR"=NA, "CCLE effect size"=NA, "CCLE pvalue"=NA, "CCLE FDR"=NA, "Reproducible"=NA)
#     known.biomarkers <- known.biomarkers[which(!is.na(known.biomarkers[ ,"type"])),]
#
#     for(i in 1:nrow(known.biomarkers)) {
#       if(!is.na(known.biomarkers[i ,"type"])) {
#         if(known.biomarkers[i ,"type"] == "expression") {
#           feature <- rownames(featureInfo(CCLE, "rna"))[which(featureInfo(CCLE, "rna")$Symbol == known.biomarkers[i ,"gene"])]
#           known.biomarkers[i ,c("CCLE effect size", "CCLE pvalue", "CCLE FDR")] <- ccle.sig.rna[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
#           known.biomarkers[i ,c("GDSC effect size", "GDSC pvalue", "GDSC FDR")] <- gdsc.sig.rna[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
#           known.biomarkers[i, "Reproducible"] <- ifelse(known.biomarkers[i ,"CCLE pvalue"] < 0.05 &
#                                                           known.biomarkers[i, "GDSC pvalue"] < 0.05 &
#                                                           sign(known.biomarkers[i ,"CCLE effect size"]) == sign(known.biomarkers[i ,"GDSC effect size"]), "YES", "NO")
#         }else if(known.biomarkers[i ,"type"] == "mutation") {
#           feature <- known.biomarkers[i ,"gene"]
#           known.biomarkers[i ,c("CCLE effect size", "CCLE pvalue", "CCLE FDR")] <- ccle.sig.mutation[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
#           known.biomarkers[i ,c("GDSC effect size", "GDSC pvalue", "GDSC FDR")] <- gdsc.sig.mutation[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
#           known.biomarkers[i, "Reproducible"] <- ifelse(known.biomarkers[i ,"CCLE pvalue"] < 0.05 &
#                                                           known.biomarkers[i, "GDSC pvalue"] < 0.05 &
#                                                           sign(known.biomarkers[i ,"CCLE effect size"]) == sign(known.biomarkers[i ,"GDSC effect size"]), "YES", "NO")
#         }else if(known.biomarkers[i ,"type"] == "fusion") {
#           feature <- known.biomarkers[i ,"gene"]
#           known.biomarkers[i ,c("GDSC effect size", "GDSC pvalue", "GDSC FDR")] <- gdsc.sig.fusion[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
#           #known.biomarkers[i, "Reproducible"] <- "NO"
#         }else if(known.biomarkers[i ,"type"] == "amplification") {
#           feature <- rownames(featureInfo(CCLE, "cnv"))[which(featureInfo(CCLE, "cnv")$Symbol == known.biomarkers[i ,"gene"])]
#           known.biomarkers[i ,c("CCLE effect size", "CCLE pvalue", "CCLE FDR")] <- ccle.sig.cnv[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
#           known.biomarkers[i ,c("GDSC effect size", "GDSC pvalue", "GDSC FDR")] <- gdsc.sig.cnv[feature, known.biomarkers[i ,"drug"], c("estimate","pvalue", "fdr")]
#           known.biomarkers[i, "Reproducible"] <- ifelse(known.biomarkers[i ,"CCLE pvalue"] < 0.05 &
#                                                           known.biomarkers[i, "GDSC pvalue"] < 0.05 &
#                                                           sign(known.biomarkers[i ,"CCLE effect size"]) == sign(known.biomarkers[i ,"GDSC effect size"]), "YES", "NO")
#         }
#       }
#     }
#     colnames(known.biomarkers)[1:3] <- capitalize(colnames(known.biomarkers)[1:3])
#     xtable::print.xtable(xtable::xtable(known.biomarkers[, c(1:3, 7:8, 10:11, 13)], digits=c(0, 0, 0, 0, 2, -1, 2, -1, 0)), include.rownames=FALSE, floating=FALSE, table.placement="!h", file=sprintf("known_biomarkers_%s_%s.tex", method, cell), append=FALSE)
#
#
#   }
#
# knownBiomarkersCheck(ccle.sig.rna=ccle.sig.rna,
#                      gdsc.sig.rna=gdsc.sig.rna2,
#                      ccle.sig.mutation=ccle.sig.mutation,
#                      gdsc.sig.mutation=gdsc.sig.mutation,
#                      gdsc.sig.fusion=gdsc.sig.fusion,
#                      ccle.sig.cnv=ccle.sig.cnv,
#                      gdsc.sig.cnv=gdsc.sig.cnv,
#                      cell="all", method="continuous")
# knownBiomarkersCheck(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2.bin, cell="common", method="continuous")
#
#
# ###Put all biomarkers in excel file
# load("signatures_data.RData")
#
# drugBasedBiomarkers(ccle.sig.rna=ccle.sig.rna, gdsc.sig.rna=gdsc.sig.rna2, cell="all", method="continuous", drugs=drugs, features=features, cut.off=0.05)
# drugBasedBiomarkers(ccle.sig.rna=common.ccle.sig.rna, gdsc.sig.rna=common.gdsc.sig.rna2, cell="common", method="continuous", drugs=drugs, features=features, cut.off=0.05)
# drugBasedBiomarkers(ccle.sig.rna=ccle.sig.mutation, gdsc.sig.rna=gdsc.sig.mutation, cell="mutation", method="continuous", drugs=drugs, features=features.mutation, cut.off=0.05)
# drugBasedBiomarkers(ccle.sig.rna=ccle.sig.cnv, gdsc.sig.rna=gdsc.sig.cnv, cell="cnv", method="continuous", drugs=drugs, features=cnv.fetures, cut.off=0.05)
#
# drugBasedBiomarkers <-
#   function(ccle.sig.rna, gdsc.sig.rna, cell, method, drugs, features, cut.off) {
#     require(WriteXLS)
#     all.biomarkers <- list()
#     for(drug in drugs) {
#       ccle.biomarkers <- ccle.sig.rna[features, drug, ]
#       colnames(ccle.biomarkers) <- paste0("CCLE_", colnames(ccle.biomarkers))
#
#       gdsc.biomarkers <- gdsc.sig.rna[features, drug, ]
#       colnames(gdsc.biomarkers) <- paste0("GDSC_", colnames(gdsc.biomarkers))
#
#       biomarkers <- cbind("Symbol"=NA, gdsc.biomarkers, ccle.biomarkers, "Specificity"="Non significant")
#       biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) < cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) < cut.off), "Specificity"] = "Both"
#       biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) < cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) >= cut.off), "Specificity"] = "CCLE"
#       biomarkers[which(as.numeric(biomarkers[,"CCLE_fdr"]) >= cut.off & as.numeric(biomarkers[,"GDSC_fdr"]) < cut.off), "Specificity"] = "GDSC"
#       #biomarkers[,"Symbol"] <- featureInfo(CCLE, "rna")[rownames(biomarkers), "Symbol"]
#       biomarkers[,"Symbol"] <- rownames(biomarkers)
#       all.biomarkers[[drug]] <- as.data.frame(biomarkers, stringsAsFactors=FALSE)
#     }
#
#     WriteXLS::WriteXLS("all.biomarkers", ExcelFileName=sprintf("all_biomarkers_%s_%s.xlsx", method, cell), row.names=TRUE)
#
#   }
