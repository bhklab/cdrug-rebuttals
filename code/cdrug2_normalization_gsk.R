########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## September 15, 2013
########################

# rm(list=ls())

require(gdata) || stop("Library gdata is not available!")
require(R.utils) || stop("Library R.utils is not available!")
require(PharmacoGx) || stop("Library PharmacoGx is not available")

path.data <- file.path("data", "GSK")
path.ge <- file.path(path.data, "ge")
path.drug <- file.path(path.data, "drug")
path.cell <- file.path(path.data, "celline")
path.mut <- file.path(path.data, "mutation")
## create directories
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.ge)) { dir.create(path.ge, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.drug)) { dir.create(path.drug, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.cell)) { dir.create(path.cell, showWarnings=FALSE, recursive=TRUE) }
if(!file.exists(path.mut)) { dir.create(path.mut, showWarnings=FALSE, recursive=TRUE) }
  

########################
## download data
########################
ftpdir <- "ftp://caftpd.nci.nih.gov/pub/caARRAY/transcript_profiling"
myfn <- file.path(path.ge, "celfile_timestamp.RData")
if(!file.exists(myfn)) {
  message("Download genomic data")
  
  dir.create(file.path(path.ge, "dwl"), showWarnings=FALSE)
  
  ## download and compress CEL files
  celfile.timestamp <- celfn <- NULL
  i <- 1
  while(i <= 8) {
    ## assuming there are only 9 zip archives (need to check if the update version has more)
   dwl.status <- download.file(url=sprintf("%s/cel/cel_0%i.zip", ftpdir, i), destfile=file.path(path.ge, "dwl", sprintf("cel_0%i.zip", i)), quiet=TRUE)
   if(dwl.status != 0) {
     message("\t-> download failed, let's try again ...")
     file.remove(file.path(path.ge, "dwl", sprintf("cel_0%i.zip", i)))
     i <- i - 1
    } else {
       ## unzip archive
       fff <- unzip(zipfile=file.path(path.ge, "dwl", sprintf("cel_0%i.zip", i)), list=TRUE)
       celfile.timestamp <- c(celfile.timestamp, as.character(fff[ ,"Date"]))
       celfn <- c(celfn, as.character(fff[ ,"Name"]))
       res <- unzip(zipfile=file.path(path.ge, "dwl", sprintf("cel_0%i.zip", i)), exdir=path.ge)
       ## compress each CEL file individually using gzip
       library(R.utils)
       sapply(file.path(path.ge, as.character(fff[ ,"Name"])), R.utils::gzip, overwrite=TRUE)
       i <- i + 1
     }
  }
  celfile.timestamp <- t(sapply(strsplit(celfile.timestamp, split=" "), function(x) { return(x) }))
  dimnames(celfile.timestamp) <- list(celfn, c("file.day", "file.hour"))
   
  # unlink(file.path(path.ge, "dwl"), recursive=TRUE)
  write.csv(celfile.timestamp, file=file.path(path.ge, "celfile_timestamp.csv"))
  save(list=c("celfile.timestamp"), compress=TRUE, file=myfn)
}

## download sample information
message("Download sample information")
myfn <- file.path(path.ge, "GSK_RNA.sdrf")
if (!file.exists(myfn)) {
  dir.create(file.path(path.ge, "dwl"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url=sprintf("%s/GSK_RNA.sdrf", ftpdir), destfile=file.path(path.ge, "dwl", "GSK_RNA.sdrf"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.ge, "dwl", "GSK_RNA.sdrf"), to=myfn)
}
  
## download drug sensitivity
message("Download drug sensitivity measurements")
myfn <- file.path(path.drug, "gsk_drug_sensitivity.xls")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "dwl"), showWarnings=FALSE)
  dwl.status <- download.file(url="http://cancerres.aacrjournals.org/content/suppl/2010/04/19/0008-5472.CAN-09-3788.DC1/stab_2.xls", destfile=file.path(path.drug, "dwl", "stab_2.xls"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.drug, "dwl", "stab_2.xls"), to=myfn)
}

########################
## normalize and format data
########################
myfn <- file.path(saveres, "gskcellines_frma.RData")
if(!file.exists(myfn)) {

  require(affy) || stop("Library affy is not available!")
  require(Hmisc) || stop("Library Hmisc is not available!")
  require(genefu) || stop("Library genefu is not available!")
  require(frma) || stop("Library frma is not available!")
  require(hgu133plus2frmavecs) || stop("Library hgu133plus2frmavecs is not available!")
  data(hgu133plus2frmavecs)
  require(hgu133plus2cdf) || stop("Library hgu133plus2cdf is not available!")
  data(hgu133plus2cdf)
  
  ## load timestamp
  load(file.path(path.ge, "celfile_timestamp.RData"))

  ## CEL file names
  celfn <- list.celfiles(path.ge, full.names=TRUE)
  celfns <- list.celfiles(path.ge, full.names=FALSE)
  ## experiments' names
  names(celfn) <- names(celfns) <- gsub(".CEL.gz", "", celfns)
  ## chip type and date
  chipt <- sapply(celfn, celfileChip)
  chipd <- t(sapply(celfn, celfileDateHour))
  ## reorder CEL files by hybridization time or timestamp
  myx <- NULL
  if(any(!complete.cases(chipd))) {
    ## all hybridization dates are not available
    load(file.path(path.ge, "celfile_timestamp.RData"))
    if(!all(is.element(celfns, paste(rownames(celfile.timestamp), "gz", sep=".")))) { stop("Timestamp is not available for all CEL files!") }
      celfile.timestamp <- celfile.timestamp[match(celfns, paste(rownames(celfile.timestamp), "gz", sep=".")), , drop=FALSE] 
      myx <- order(celfile.timestamp[ ,"file.day"], celfile.timestamp[ ,"file.hour"], decreasing=FALSE)
  } else {
      myx <- order(chipd[ ,"day"], chipd[ ,"hour"], decreasing=FALSE)
  }
  celfn <- celfn[myx]
  celfns <- celfns[myx]
  chipt <- chipt[myx]
  chipd <- chipd[myx, , drop=FALSE]
  celfile.timestamp <- celfile.timestamp[myx, , drop=FALSE]

  ## read info about drugs and experiments

  ## info about each experiment
  message("Read sample information")
  sampleinfo <- read.csv(file.path(path.ge, "GSK_RNA.sdrf"), sep="\t", comment.char="#")
  sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
  ## curate cell line names
  rownames(sampleinfo) <- gsub(" - Replicate ", "_rep", sampleinfo[ , "Source.Name"])
  sampleinfo <- data.frame("samplename"=rownames(sampleinfo), "cellid"=sampleinfo[ , "Characteristics.Cell.Line.Name."], "filename"=sprintf("%s.gz", sampleinfo[ , "Array.Data.File"]), "tissue.type"=tolower(gsub("_$", "", gsub(badchars, "_", sampleinfo[ , "Characteristics.OrganismPart."]))), sampleinfo)
  
  ## read drug phenotypes
  drugpheno <- read.xls(xls=file.path(path.drug, "gsk_drug_sensitivity.xls"), sheet=1)
  drugpheno[drugpheno == ""] <- NA
  cn <- gsub(badchars, ".", drugpheno[6, ])
  cn <- gsub("CL.ID", "CL_ID", cn)
  drugpheno <- drugpheno[-(1:6), , drop=FALSE]
  drugpheno <- drugpheno[!is.na(drugpheno[ , 2]), , drop=FALSE]
  dimnames(drugpheno) <- list(drugpheno[ , 2], cn)
  celline.info <- drugpheno[ , c("CL_ID", "Cell.Line", "Site", "Dx"), drop=FALSE]
  mutation <- drugpheno[ , c("HRAS.", "KRAS.", "NRAS.", "BRAF.", "PIK3CA.", "PTEN."), drop=FALSE]
  drugpheno <- data.matrix(drugpheno[ , grep("IC50", cn), drop=FALSE])
  ## IC50 in nano molar

  ## keep only the CEL files present in sampleinfo
  myx <- intersect(sampleinfo[ , "filename"], celfns)
  celfn <- celfn[match(myx, celfns)]
  celfns <- celfns[match(myx, celfns)]
  sampelinfo <- sampleinfo[match(myx, sampleinfo[ , "filename"]), , drop=FALSE]

  ## normalization
  message("Normalize gene expression data")
  # rr <- just.rma(filenames=celfn)
  ## frma normalization using parallel
  splitix <- parallel::splitIndices(nx=length(celfn), ncl=nbcore)
  splitix <- splitix[sapply(splitix, length) > 0]
  res <- parallel::mclapply(splitix, function(x, celfn) {
    ## fRMA
    tt <- celfn[x]
    names(tt) <- NULL
    abatch <- affy::read.affybatch(filenames=tt)
    rr <- frma::frma(object=abatch, summarize="robust_weighted_average", verbose=TRUE, input.vecs=hgu133plus2frmavecs)
    rr <- exprs(rr)
  }, celfn=celfn)
  datat <- t(do.call(cbind, res))
  ## match the experiment labels
  rownames(datat) <- rownames(sampleinfo)[match(rownames(datat), as.character(sampleinfo[ ,"filename"]))]

  ## build annotation matrix
  message("Build annotation matrix")
  myfn2 <- file.path(saveres, "gsk_ge_annot.RData")
  if(!file.exists(myfn2)) {
    ## select the best probe for a single gene
    library(jetset)
    js <- jetset::jscores(chip="hgu133plus2", probeset=colnames(genexprs))
    js <- js[colnames(genexprs), , drop=FALSE]
    ## identify the best probeset for each entrez gene id
    geneid1 <- as.character(js[ ,"EntrezID"])
    names(geneid1) <- rownames(js)
    geneid2 <- sort(unique(geneid1))
    names(geneid2) <- paste("geneid", geneid2, sep=".")
    gix1 <- !is.na(geneid1)
    gix2 <- !is.na(geneid2)
    geneid.common <- intersect(geneid1[gix1], geneid2[gix2])
    ## probes corresponding to common gene ids
    gg <- names(geneid1)[is.element(geneid1, geneid.common)]
    gid <- geneid1[is.element(geneid1, geneid.common)]
    ## duplicated gene ids
    gid.dupl <- unique(gid[duplicated(gid)])
    gg.dupl <- names(geneid1)[is.element(geneid1, gid.dupl)]
    ## unique gene ids
    gid.uniq <- gid[!is.element(gid, gid.dupl)]
    gg.uniq <- names(geneid1)[is.element(geneid1, gid.uniq)]
    ## which are the best probe for each gene
    js <- data.frame(js, "best"=FALSE)
    js[gg.uniq, "best"] <- TRUE
    ## data for duplicated gene ids
    if(length(gid.dupl) > 0) {
    	library(jetset)
    	## use jetset oevrall score to select the best probeset
    	myscore <- js[gg.dupl,"overall"]
    	myscore <- cbind("probe"=gg.dupl, "gid"=geneid1[gg.dupl], "score"=myscore)
    	myscore <- myscore[order(as.numeric(myscore[ , "score"]), decreasing=TRUE, na.last=TRUE), , drop=FALSE]
    	myscore <- myscore[!duplicated(myscore[ , "gid"]), , drop=FALSE]
    	js[myscore[ ,"probe"], "best"] <- TRUE
    }
    annot <- data.frame("probe"=rownames(js), "EntrezGene.ID"=js[ ,"EntrezID"], js)
    annot <- annot[colnames(genexprs), , drop=FALSE]
    save(list=c("annot"), compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  
  ########################
  ## merge data and put NA when the data are anot avaialble
  ########################
  message("Merge data")
  
  ## cell lines and drugs
  cellnall <- sort(unique(c(as.character(sampleinfo[ ,"cellid"]), as.character(drugpheno[ ,"cellid"]), rownames(genexprs), rownames(mutation))))
  drugnall <- sort(unique(c(rownames(druginfo), as.character(drugconc[ , "drugid"]), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
  
  ## update sampleinfo
  dd <- data.frame(matrix(NA, nrow=length(cellnall), ncol=ncol(sampleinfo), dimnames=list(cellnall, colnames(sampleinfo))))
  newlev <- sapply(sampleinfo, levels)
  newlev$cellid <- cellnall
  dd <- genefu::setcolclass.df(df=dd, colclass=sapply(sampleinfo, class), factor.levels=newlev)
  dd[rownames(sampleinfo), colnames(sampleinfo)] <- sampleinfo
  dd[ , "cellid"] <- rownames(dd)
  sampleinfo<- dd
  
  ## update gene expression
  dd <- matrix(NA, ncol=ncol(genexprs), nrow=length(cellnall), dimnames=list(cellnall, colnames(genexprs)))
  dd[rownames(genexprs),colnames(genexprs)] <- genexprs
  genexprs <- dd
  
  ## update mutation
  dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
  dd[rownames(mutation),colnames(mutation)] <- mutation
  mutation <- dd
  ## update drugpheno
  ## IC50 in microM
  iix <- grep("_IC50$", colnames(drugpheno))
  iixn <- gsub("_IC50$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.ic50 <- dd
  ## EC50 in microM
  iix <- grep("_EC50$", colnames(drugpheno))
  iixn <- gsub("_EC50$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.ec50 <- dd
  ## Amax in microM
  iix <- grep("_Amax$", colnames(drugpheno))
  iixn <- gsub("_Amax$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.amax <- dd
  ## AUC (higher value represents sensitivity)
  iix <- grep("_ActivityArea$", colnames(drugpheno))
  iixn <- gsub("_ActivityArea$", "", colnames(drugpheno)[iix])
  auct <- drugpheno[ , iix]
  colnames(auct) <- iixn
  ## division by the number of concentrations tested
  myx <- sapply(strsplit(colnames(drugpheno), "_"), function(x) { return(x[length(x)] == "Doses") })
  ndoses <- drugpheno[ ,myx,drop=FALSE]
  ndoses <- apply(ndoses, 2, function(x) {
    return(sapply(x, function(x) {
      return(sapply(strsplit(as.character(x), ","), function(x) { if(is.na(x[1])) { return(NA) } else { return(length(x)) } }))
    }))
  })
  colnames(ndoses) <- gsub("_Doses", "", colnames(ndoses))
  ndoses <- ndoses[rownames(auct), colnames(auct)]
  auct <- auct / ndoses
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(auct), colnames(auct)] <- data.matrix(auct)
  dd[dd < 0] <- 0
  dd[dd > 1] <- 1
  drugpheno.auc <- dd
  
  ## compute drug sensitivity calls using waterfall plot
  ## IC50
  ic50.call2 <- ic50.call3 <- data.frame(matrix(NA, ncol=ncol(drugpheno.ic50), nrow=nrow(drugpheno.ic50), dimnames=dimnames(drugpheno.ic50)))
  ## IC50 sensitivity calling 2 levels
  pdf(file.path(saveres, "ccle_ic50_sensitivity_calling2_drugs.pdf"), width=10, height=10)
  for(i in 1:ncol(drugpheno.ic50)) {
    ic50.call2[ , i] <- PharmacoGx::callingWaterfall(x=drugpheno.ic50[ , i] / 10^6, type="IC50", intermediate.fold=0, cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CCLE)", colnames(drugpheno.ic50)[i]))
  }
  dev.off()
  ## IC50 sensitivity calling 3 levels
  pdf(file.path(saveres, "ccle_ic50_sensitivity_calling3_drugs.pdf"), width=10, height=10)
  for(i in 1:ncol(drugpheno.ic50)) {
    ic50.call3[ , i] <- PharmacoGx::callingWaterfall(x=drugpheno.ic50[ , i] / 10^6, type="IC50", intermediate.fold=4, cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CCLE)", colnames(drugpheno.ic50)[i]))
  }
  dev.off()
  dimnames(ic50.call2) <- dimnames(ic50.call3) <- dimnames(drugpheno.ic50)
  ## AUC
  auc.call2 <- auc.call3 <- data.frame(matrix(NA, ncol=ncol(drugpheno.auc), nrow=nrow(drugpheno.auc), dimnames=dimnames(drugpheno.auc)))
  ## AUC sensitivity calling 2 levels
  pdf(file.path(saveres, "ccle_auc_sensitivity_calling2_drugs.pdf"), width=10, height=10)
  for(i in 1:ncol(drugpheno.auc)) {
    auc.call2[ , i] <- PharmacoGx::callingWaterfall(x=drugpheno.auc[ , i] / 10^6, type="AUC", intermediate.fold=0, cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CCLE)", colnames(drugpheno.auc)[i]))
  }
  dev.off()
  ## AUC sensitivity calling 3 levels
  pdf(file.path(saveres, "ccle_auc_sensitivity_calling3_drugs.pdf"), width=10, height=10)
  for(i in 1:ncol(drugpheno.auc)) {
    auc.call3[ , i] <- PharmacoGx::callingWaterfall(x=drugpheno.auc[ , i] / 10^6, type="AUC", intermediate.fold=1.2, cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CCLE)", colnames(drugpheno.auc)[i]))
  }
  dev.off()
  dimnames(auc.call2) <- dimnames(auc.call3) <- dimnames(drugpheno.auc)

  ## save drug phenotypes
  drugpheno.orig <- drugpheno
  drugpheno <- list("AUC"=drugpheno.auc, "IC50"=drugpheno.ic50, "EC50"=drugpheno.ec50, "AMAX"=drugpheno.amax, "IC50.CALL2"=ic50.call2, "IC50.CALL3"=ic50.call3, "AUC.CALL2"=auc.call2, "AUC.CALL3"=auc.call3)
  
  ## update drugconc
  drugconcnall <- as.vector(outer(drugnall, cellnall, paste, sep="..."))
  dd <- data.frame(matrix(NA, nrow=length(drugconcnall), ncol=ncol(drugconc), dimnames=list(drugconcnall, colnames(drugconc))))
  newlev <- sapply(drugconc, levels)
  newlev$drugid <- drugnall
  newlev$cellid <- cellnall
  dd <- genefu::setcolclass.df(df=dd, colclass=sapply(drugconc, class), factor.levels=newlev)
  dd[rownames(drugconc), colnames(drugconc)] <- drugconc
  dd[ , "cellid"] <- sapply(strsplit(x=rownames(dd), split="[.][.][.]"), function (x) { return (x[2]) })
  dd[ , "drugid"] <- sapply(strsplit(x=rownames(dd), split="[.][.][.]"), function (x) { return (x[1]) })
  drugconc <- dd
  
  ## update druginfo
  dd <- data.frame(matrix(NA, nrow=length(drugnall), ncol=ncol(druginfo), dimnames=list(drugnall, colnames(druginfo))))
  newlev <- sapply(druginfo, levels)
  newlev$drugid <- sapply(strsplit(drugnall, split="_"), function(x) { return(x[2]) })
  dd <- genefu::setcolclass.df(df=dd, colclass=sapply(druginfo, class), factor.levels=newlev)
  dd[match(rownames(druginfo), drugnall), colnames(druginfo)] <- druginfo
  dd[ , "drugid"] <- newlev$drugid
  druginfo <- dd
  
  ## rename objects
  data.ge.ccle <- genexprs
  annot.ge.ccle <- annot
  sampleinfo.ccle <- sampleinfo
  drugpheno.ccle <- drugpheno
  druginfo.ccle <- druginfo
  drugconc.ccle <- drugconc
  mutation.ccle <- mutation

  ## make sure that cellid are not factors
  celline.ccle[, "cellid"] <- as.character(celline.ccle[, "cellid"])
  sampleinfo.ccle[, "cellid"] <- as.character(sampleinfo.ccle[, "cellid"])
  drugconc.ccle[, "cellid"] <- as.character(drugconc.ccle[, "cellid"])

  message("Save data")
  write.csv(celline.ccle, file=file.path(saveres, "cell_line_collection_ccle.csv"), row.names=FALSE)
  write.csv(annot.ge.ccle, file=file.path(saveres, "annot_ge_ccle.csv"), row.names=FALSE)
  write.csv(t(data.ge.ccle), file=file.path(saveres, "data_ge_ccle.csv"))
  write.csv(drugpheno.orig, file=file.path(saveres, "drugpheno_ccle.csv"), row.names=FALSE)
  write.csv(drugpheno.ccle$IC50, file=file.path(saveres, "drugpheno_ic50_ccle.csv"), row.names=TRUE)
  write.csv(drugpheno.ccle$EC50, file=file.path(saveres, "drugpheno_ec50_ccle.csv"), row.names=TRUE)
  write.csv(drugpheno.ccle$AMAX, file=file.path(saveres, "drugpheno_amax_ccle.csv"), row.names=TRUE)
  write.csv(drugpheno.ccle$AUC, file=file.path(saveres, "drugpheno_auc_ccle.csv"), row.names=TRUE)
  write.csv(druginfo.ccle, file=file.path(saveres, "druginfo_ccle.csv"), row.names=FALSE)
  write.csv(drugconc.ccle, file=file.path(saveres, "drugconc_ccle.csv"), row.names=FALSE)
  write.csv(sampleinfo.ccle, file=file.path(saveres, "sampleinfo_ccle.csv"), row.names=FALSE)
  write.csv(mutation.ccle, file=file.path(saveres, "mutation_ccle.csv"), row.names=TRUE)
  save(list=c("data.ge.ccle", "annot.ge.ccle", "sampleinfo.ccle", "drugpheno.ccle", "druginfo.ccle", "drugconc.ccle", "mutation.ccle", "celline.ccle"), compress=TRUE, file=myfn)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  sampleinfo.ge.gsk <- sampleinfo
  data.ge.gsk <- datat
  annot.ge.gsk <- annot
  drugpheno.gsk <- drugpheno
  mutation.gsk <- mutation
  
  save(list=c("data.ge.gsk", "annot.ge.gsk", "sampleinfo.ge.gsk", "drugpheno.gsk", "mutation.gsk"), compress=TRUE, file=file.path(saveres, "gskcellines_complete_frma.RData"))
  
  ## average replicates
  message("Averaging replicates")
  ucell <- as.character(unique(sampleinfo.ge.gsk[ , "cellid"]))
  splitix <- parallel::splitIndices(nx=length(ucell), ncl=nbcore)
  res <- parallel::mclapply(splitix, function(x, ucell, data) {
    dd <- t(sapply(ucell[x], function(x, y) {
      iix <- grep(x, rownames(y))
      dd <- apply(y[iix, , drop=FALSE], 2, mean, na.rm=FALSE)
      return(dd)
    }, y=data))
    # dd <- matrix(NA, ncol=ncol(data.ge.gsk), nrow=length(ucell), dimnames=list(ucell, colnames(data.ge.gsk)))
    # for(i in 1:length(ucell)) {
    #   iix <- grep(ucell[i], rownames(data.gsk))
    #   dd[ucell[i], ] <- apply(data.ge.gsk[iix, , drop=FALSE], 2, mean, na.rm=FALSE)
    # }
    return(dd)
  }, ucell=ucell, data=data.ge.gsk)
  data.ge.gsk <- do.call(rbind, res)
  sampleinfo.ge.gsk <- sampleinfo.ge.gsk[!duplicated(sampleinfo.ge.gsk[ , "cellid"]), , drop=FALSE]
  rownames(sampleinfo.ge.gsk) <- as.character(sampleinfo.ge.gsk[ , "cellid"])
  sampleinfo.ge.gsk <- sampleinfo.ge.gsk[rownames(data.ge.gsk), , drop=FALSE]
  
  ## match cell line names
  ## no hits with partial matching
  # myx <- !is.element(rownames(sampleinfo.ge.gsk), gsub(badchars, "", rownames(drugpheno.gsk)))
  # myx2 <- !is.element(gsub(badchars, "", rownames(drugpheno.gsk)), rownames(sampleinfo.ge.gsk))
  # tt <- lapply(rownames(sampleinfo.ge.gsk)[myx], agrep, x=rownames(drugpheno.gsk)[myx2], value=TRUE)
  # tt <- mapply(function(x, y) { return(sprintf("%s -> %s", x, paste(y, collapse=":::"))) }, x=as.list(rownames(sampleinfo.ge.gsk)[myx]), y=tt)
  # write.csv(tt, "temp.csv")
  drugpheno.ge.gsk <- matrix(NA, ncol=ncol(drugpheno.gsk), nrow=nrow(sampleinfo.ge.gsk), dimnames=list(rownames(sampleinfo.ge.gsk), colnames(drugpheno.gsk)))
  matchix <- match(rownames(drugpheno.ge.gsk), gsub(badchars, "", rownames(drugpheno.gsk)))
  matchix.na <- is.na(matchix)
  drugpheno.ge.gsk[!matchix.na, ] <- drugpheno.gsk[matchix[!matchix.na], , drop=FALSE]
  colnames(drugpheno.ge.gsk) <- paste("drugid", toupper(sapply(strsplit(gsub(badchars, "", colnames(drugpheno.ge.gsk)), split="[(]gIC50"), function(x) { return(x[[1]]) })), sep="_")
  ## transofrm into microM
  drugpheno.ge.gsk <- drugpheno.ge.gsk / 1000
  ## drug concentrations
  iix <- which(!is.na(drugpheno.ge.gsk))
  posix <- pos2coord(pos=iix, dim.mat=dim(drugpheno.ge.gsk))
  drugconc.ge.gsk <- t(apply(posix, 1, function(x, y, z) {
    cc <- c(0.0003, 0.0032, 0.01, 0.032, 0.1, 0.317, 1, 3.16, 10)
    names(cc) <- sprintf("Dose%i.uM", 1:length(cc))
    return(c("cellid"=y[x[1]], "drugid"=z[x[2]], "nbr.conc.tested"=length(cc), cc))
  }, y=rownames(drugpheno.ge.gsk), z=colnames(drugpheno.ge.gsk)))
  rownames(drugconc.ge.gsk) <- paste(drugconc.ge.gsk[ , "cellid"], drugconc.ge.gsk[ , "drugid"], sep="...")
  
  message("Save GSK data")
  write.csv(annot.ge.gsk, file=file.path(saveres, "annot_ge_gsk.csv"), row.names=FALSE)
  write.csv(t(data.ge.gsk), file=file.path(saveres, "data_ge_gsk.csv"))
  write.csv(sampleinfo.ge.gsk, file=file.path(saveres, "sampleinfo_ge_gsk.csv"), row.names=FALSE)
  write.csv(drugpheno.ge.gsk, file=file.path(saveres, "drugpheno_ge_gsk.csv"))
  write.csv(drugconc.ge.gsk, file=file.path(saveres, "drugconc_ge_gsk.csv"))
  save(list=c("data.ge.gsk", "annot.ge.gsk", "sampleinfo.ge.gsk", "drugconc.ge.gsk", "drugpheno.ge.gsk", "mutation.gsk"), compress=TRUE, file=myfn)
}





## end




