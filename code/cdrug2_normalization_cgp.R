########################
## Benjamin Haibe-Kains
## Code under License Artistic-2.0
## April 8, 2014
########################

# rm(list=ls())

path.data <- file.path("data", "CGP")
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

require(gdata) || stop("Library gdata is not available!")
require(R.utils) || stop("Library R.utils is not available!")
require(PharmacoGx) || stop("Library PharmacoGx is not available")
require(XML) || stop("Library XML is not available")


########################
## download data
########################
ftpdir <- "ftp://ftp.ebi.ac.uk//pub/databases/microarray/data/experiment/MTAB/E-MTAB-783/"
myfn <- file.path(path.ge, "celfile_timestamp.RData")
if(!file.exists(myfn)) {
  message("Download genomic data")
  
  require(R.utils) || stop("Library R.utils is not available!")
  
  dir.create(file.path(path.ge, "dwl"), showWarnings=FALSE, recursive=TRUE)
  
  ## download and compress CEL files
  celfile.timestamp <- celfn <- NULL
  i <- 1
  while(i <= 9) {
    ## assuming there are only 9 zip archives (need to check if the update version has more)
   dwl.status <- download.file(url=sprintf("%s/E-MTAB-783.raw.%i.zip", ftpdir, i), destfile=file.path(path.ge, "dwl", sprintf("E-MTAB-783.raw.%i.zip", i)), quiet=TRUE)
   if(dwl.status != 0) {
     message("\t-> download failed, let's try again ...")
     file.remove(file.path(path.ge, "dwl", sprintf("E-MTAB-783.raw.%i.zip", i)))
     i <- i - 1
    } else {
       ## unzip archive
       fff <- unzip(zipfile=file.path(path.ge, "dwl", sprintf("E-MTAB-783.raw.%i.zip", i)), list=TRUE)
       celfile.timestamp <- c(celfile.timestamp, as.character(fff[ ,"Date"]))
       celfn <- c(celfn, as.character(fff[ ,"Name"]))
       res <- unzip(zipfile=file.path(path.ge, "dwl", sprintf("E-MTAB-783.raw.%i.zip", i)), exdir=path.ge)
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
} else { load(myfn) }

## download sample information
message("Download sample information")
myfn <- file.path(path.ge, "cgp_ge_sampleinfo.txt")
if (!file.exists(myfn)) {
  dir.create(file.path(path.ge, "dwl"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url=sprintf("%s/E-MTAB-783.sdrf.txt", ftpdir), destfile=file.path(path.ge, "dwl", "E-MTAB-783.sdrf.txt"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.ge, "dwl", "E-MTAB-783.sdrf.txt"), to=myfn)
}
  
## download drug sensitivity
message("Download drug sensitivity measurements")
myfn <- file.path(path.drug, "cgp_drug_sensitivity.csv")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "dwl"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_manova_input_w5.csv", destfile=file.path(path.drug, "dwl", "gdsc_manova_input_w5.csv"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.drug, "dwl", "gdsc_manova_input_w5.csv"), to=myfn)
}

## download drug concentration
message("Download screening drug concentrations")
myfn <- file.path(path.drug, "cgp_drug_concentration.csv")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "dwl"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_compounds_conc_w5.csv", destfile=file.path(path.drug, "dwl", "gdsc_compounds_conc_w5.csv"), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
  file.copy(from=file.path(path.drug, "dwl", "gdsc_compounds_conc_w5.csv"), to=myfn)
}

## download cell line annotations and COSMIC IDs
## annotations from COSMIC cell line project
myfn <- file.path(path.cell, "tmp", "cosmic_annotations.RData")
if(!file.exists(myfn)) {
  message("Download COSMIC annotations for cell lines")
  myfn2 <- file.path(path.cell, "tmp", "cosmic_cell_line_collection.txt")
  if(!file.exists(myfn2)) {
    dir.create(file.path(path.cell, "tmp"), showWarnings=FALSE, recursive=TRUE)
    dwl.status <- getCosmic(em="bhk.labgroup@gmail.com", passw="pharmacogenomics", directory=file.path(path.cell, "tmp"))
    # dwl.status <- download.file(url=sprintf("http://cancer.sanger.ac.uk/files/cosmic/current_release/CosmicCompleteExport.tsv.gz"), destfile=file.path(path.cell, "tmp", sprintf("CosmicCompleteExport.tsv.gz")), quiet=TRUE)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline") }
    ## untar
    res <- R.utils::gunzip(filename=file.path(path.cell, "tmp", sprintf("CosmicCompleteExport.tsv.gz")), overwrite=TRUE)
    file.copy(from=file.path(path.cell, "tmp", "CosmicCompleteExport.tsv"), to=myfn2)
  }
  message("Process COSMIC annotations")
  cosmic.celline <- read.csv(file=file.path(path.cell, "tmp", "cosmic_cell_line_collection.txt"), sep="\t")
  # cosmic.celline <- cosmic.celline[- c(grep("row selected", cosmic.celline[ ,1]), grep("rows selected", cosmic.celline[ ,1])), , drop=FALSE]
  cosmic.celline <- cosmic.celline[complete.cases(cosmic.celline[ , c("Sample.name", "Sample.source")]) & cosmic.celline[ , "Sample.source"] == "cell-line", , drop=FALSE]
  cosmic.celline[cosmic.celline == "NS" | cosmic.celline == "" | cosmic.celline == " " | cosmic.celline == "  "] <- NA
  ## merge the gene targets
  dupln <- sort(unique(cosmic.celline[ , "Sample.name"][duplicated(cosmic.celline[ , "Sample.name"])]))
  tt <- cosmic.celline
  ## select unique cell lines
  iix.rm <- NULL
  for(i in 1:length(dupln)) {
    duplix <- cosmic.celline[ ,"Sample.name"] == dupln[i]
    iix <- sort((which(duplix)), decreasing=FALSE)[1]
    iix.rm <- c(iix.rm, setdiff(which(duplix), iix))
    ## get the most frequent tissue type
    tissuet <- table(cosmic.celline[duplix, "Primary.site"])
    if (length(tissuet) == 0) {
      tt[iix, "Primary.site"] <- NA
    } else {
      tt[iix, "Primary.site"] <- names(sort(tissuet, decreasing=TRUE))[1]
    }
    # tt[iix, "Gene.name"] <- paste(cosmic.celline[duplix, "Gene.name"], collapse="///")
    # tt[iix, "UniProt.ID"] <- paste(cosmic.celline[duplix, "UniProt.ID"], collapse="///")
    # tt[iix, "Zygosity"] <- paste(cosmic.celline[duplix, "Zygosity"], collapse="///")
    # tt[iix, "CDS_MUT_SYNTAX"] <- paste(cosmic.celline[duplix, "CDS_MUT_SYNTAX"], collapse="///")
    # tt[iix, "AA_MUT_SYNTAX"] <- paste(cosmic.celline[duplix, "AA_MUT_SYNTAX"], collapse="///")
    # tt[iix, "NCBI36.genome.position"] <- paste(cosmic.celline[duplix, "NCBI36.genome.position"], collapse="///")
    # tt[iix, "GRCh37.genome.position"] <- paste(cosmic.celline[duplix, "GRCh37.genome.position"], collapse="///")
  }
  tt <- tt[-iix.rm, , drop=FALSE]
  tt <- tt[!is.na(tt[ , "Sample.name"]), , drop=FALSE]
  rownames(tt) <- tt[ , "Sample.name"]
  ## remove unnecessary annotations
  tt <- tt[ , c("Sample.name", "ID_sample", "ID_tumour", "Primary.site", "Site.subtype", "Primary.histology", "Histology.subtype", "Sample.source", "Tumour.origin", "Comments"), drop=FALSE]
  cosmic.celline <- tt
  save(list=c("cosmic.celline"), compress=TRUE, file=myfn)
} else { load(myfn) }

## annotations from GDSC (Genomics of Drug Sensitivity in Cancer)
myfn <- file.path(saveres, "gdsc_annotations.RData")
if(!file.exists(myfn)) {
  message("Download GDSC annotations for cell liness")
  myfn2 <- file.path(path.cell, "cgp_celline_collection.csv")
  if(!file.exists(myfn2)) {
    dir.create(file.path(path.cell, "dwl"), showWarnings=FALSE, recursive=TRUE)
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_cell_lines_w5.csv", destfile=file.path(path.cell, "dwl", "gdsc_cell_lines_w5.csv"), quiet=TRUE)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    file.copy(from=file.path(path.cell, "dwl", "gdsc_cell_lines_w5.csv"), to=myfn2)
  }
  gdsc.celline <- read.csv(file=file.path(path.cell, "cgp_celline_collection.csv"))
  gdsc.celline[gdsc.celline == "" | gdsc.celline == " " | gdsc.celline == "  "] <- NA
  gdsc.celline <- gdsc.celline[!is.na(gdsc.celline[ , "CELL_LINE_NAME"]), , drop=FALSE]
  dupln <- unique(gdsc.celline[ , "CELL_LINE_NAME"][duplicated(gdsc.celline[ , "CELL_LINE_NAME"])])
  gdsc.celline <- gdsc.celline[!duplicated(gdsc.celline[ , "CELL_LINE_NAME"]), , drop=FALSE]
  rownames(gdsc.celline) <- gdsc.celline[ , "CELL_LINE_NAME"]
  save(list=c("gdsc.celline"), compress=TRUE, file=myfn)
} else { load(myfn) }

## merge GDSC and COSMIC annotations through COSMIC_ID
message("Merge COSMIC and GDSC annotations for cell liness")
iix <- which(complete.cases(gdsc.celline[ , c("CELL_LINE_NAME", "COSMIC_ID")]) & !is.element(gdsc.celline[ , "COSMIC_ID"], cosmic.celline[ , "ID_sample"]) & !is.element(gdsc.celline[ , "CELL_LINE_NAME"], cosmic.celline[ , "Sample.name"]))
tt <- data.frame(matrix(NA, nrow=nrow(cosmic.celline) + length(iix), ncol=ncol(cosmic.celline), dimnames=list(c(rownames(cosmic.celline), rownames(gdsc.celline)[iix]), colnames(cosmic.celline))))
tt[rownames(cosmic.celline), ] <- cosmic.celline
tt[rownames(gdsc.celline)[iix], "Sample.name"] <- gdsc.celline[iix, "CELL_LINE_NAME"]
tt[rownames(gdsc.celline)[iix], "ID_sample"] <- gdsc.celline[iix, "COSMIC_ID"]
celline.cgp <- tt

## download drug information
message("Download drug information")
myfn <- file.path(path.drug, "cgp_drug_information.csv")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "dwl"), showWarnings=FALSE, recursive=TRUE)
  # dwl.status <- download.file(url="http://www.cancerrxgene.org/action/ExportJsonTable/CSV", destfile=file.path(path.drug, "dwl", "export-Automatically_generated_table_data.csv"), quiet=TRUE)
  # if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }  
  tables <- XML::readHTMLTable("http://www.cancerrxgene.org/translation/Drug")
  drugs <- tables[1][[1]]
  write.csv(drugs, row.names=FALSE, file=file.path(path.drug, "dwl", "export.csv"))
  file.copy(from=file.path(path.drug, "dwl", "export.csv"), to=myfn)
}
myfn <- file.path(path.drug, "nature_supplementary_information.xls")
if (!file.exists(myfn)) {
  dir.create(file.path(path.drug, "dwl"), showWarnings=FALSE, recursive=TRUE)
  dwl.status <- download.file(url="http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip", destfile=file.path(path.drug, "dwl", "nature11005-s2.zip"), quiet=TRUE)
  ff <- as.character(unzip(zipfile=file.path(path.drug, "dwl", "nature11005-s2.zip"), list=TRUE)[1, 1])
  unzip(zipfile=file.path(path.drug, "dwl", "nature11005-s2.zip"), exdir=file.path(path.drug, "dwl"))
  file.copy(from=file.path(path.drug, "dwl", ff), to=myfn)
}

########################
## normalize and format data
########################
myfn <- file.path(saveres, "cgp_data.RData")
if(!file.exists(myfn)) {

  require(affy) || stop("Library affy is not available!")
  require(Hmisc) || stop("Library Hmisc is not available!")
  require(genefu) || stop("Library genefu is not available!")
  require(frma) || stop("Library frma is not available!")
  require(hthgu133afrmavecs) || stop("Library hthgu133afrmavecs is not available!")
  data(hthgu133afrmavecs)
  require(hthgu133acdf) || stop("Library hthgu133acdf is not available!")
  data(hthgu133acdf)

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

  ## phenotype for the drugs
  message("Read drug sensitivity measurements")
  myfn2 <- file.path(saveres, "cgp_drug_sensitivity.RData")
  if(!file.exists(myfn2)) {
    drugpheno <- read.csv(file.path(path.drug, "cgp_drug_sensitivity.csv"))
    drugpheno[drugpheno == "" | drugpheno == " "] <- NA
    save(list="drugpheno", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  ## format column names
  coln2 <- unlist(drugpheno[1, ,drop=TRUE])
  coln2[coln2 == ""] <- NA
  drugpheno <- drugpheno[!is.na(drugpheno[ , "Cell.Line"]), ,drop=FALSE]
  coln <- colnames(drugpheno)
  coln2[is.na(coln2)] <- coln[is.na(coln2)]
  coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
  myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
  coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
  colnames(drugpheno) <- coln2
  ## drug identifiers and names
  dn <- toupper(gsub(badchars, "", sapply(strsplit(coln, "_"), function(x) { return(x[[1]]) })))
  ## manual curation for drug names starting with a figure
  dn[!is.na(dn) & dn == "X17AAG"] <- "17AAG"
  dn[!is.na(dn) & dn == "X681640"] <- "681640"
  did <- sapply(strsplit(coln2, "_"), function(x) { if(x[[1]] == "drugid") { return(x[[2]]) } else { return(NA) } })
  drugnid <- cbind("drug.name"=dn, "drug.id"=did)[!is.na(did) & !duplicated(did), ]
  rownames(drugnid) <- paste("drugid", drugnid[ , "drug.id"], sep="_")

  ## cell line identifiers
  dupln <- duplicated(drugpheno[ ,"Cell.Line"])
  if(sum(dupln) > 1) { warning("some cell lines are duplicated, only the first instance is kept") }
  drugpheno <- drugpheno[!dupln, , drop=FALSE]
  
  if(any(!is.element(drugpheno[ ,"Cell.Line"], celline.cgp[ , "Sample.name"]))) { stop("Some cell line names are not included in the COSMIC database") }
  celln <- drugpheno[ ,"Cell.Line"]
  drugpheno <- data.frame("cellid"=celln, drugpheno)
  rownames(drugpheno) <- celln
  
  ## protein coding variants
  ## Genetic mutation data for cancer genes. Includes MSI status (1 = unstable and 0 = stable) and gene-fusions. A binary code 'x::y' description is used for each gene where 'x' identifies a coding variant and 'y' indicates copy number information from SNP6.0 data. For gene fusions, cell lines are identified as fusion not-detected (0) or the identified fusion is given. The following abbreviations are used: not analysed (na), not detected or wild-type (wt), no copy number information (nci).
  ## we assume that AKT2 and WT1 are the first and last genes in the file
  rangeg <- which(colnames(drugpheno) == "AKT2"):which(colnames(drugpheno) == "MLL_AFF1")
  mutation <- as.matrix(drugpheno[ , rangeg, drop=FALSE])
  mutation <- apply(X=mutation, MARGIN=c(1, 2), FUN=function(x) {
    x <- unlist(strsplit(x, split="::"))
    if(length(x) == 2) {
      if(!is.na(x[[1]]) && (x[[1]] == "na")) {
        x <- NA
      } else {
        x <- x[[1]]
      }
    } else { x <- NA }
    return(x)
  })
  write.csv(mutation, file=file.path(path.mut, "gdsc_mutation.csv"))

  ## info about each experiment
  message("Read sample information")
  sampleinfo <- read.csv(file.path(path.ge, "cgp_ge_sampleinfo.txt"), sep="\t")
  sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
  ## curate cell line names
  sampleinfo[sampleinfo[ , "Source.Name"] == "MZ2-MEL.", "Source.Name"] <- "MZ2-MEL"
  iix <- which(!duplicated(sampleinfo[ , "Source.Name"]) & !is.element(sampleinfo[ , "Source.Name"], celline.cgp[ , "Sample.name"]))
  if(length(iix) > 0) {
    ## enrich the list of cell lines
    tt <- matrix(NA, nrow=length(iix), ncol=ncol(celline.cgp), dimnames=list(sampleinfo[iix, "Source.Name"], colnames(celline.cgp)))
    tt[ , "Sample.name"] <- sampleinfo[iix, "Source.Name"]
    celline.cgp <- rbind(celline.cgp, tt)
  }
  fn <- gsub(patter="[.]CEL", replacement="", x=sampleinfo[ ,"Array.Data.File"])
  if(any(!is.element(fn[!is.na(fn)], names(celfns)))) { stop("some CEL files are missing for the CGP project") }
  rownames(sampleinfo) <- fn
  sampleinfo <- sampleinfo[names(celfn), , drop=FALSE]
  sampleinfo <- data.frame("samplename"=names(celfns), "filename"=celfns, "chiptype"=chipt, "hybridization.date"=chipd[ ,"day"], "hybridization.hour"=chipd[ ,"hour"], "file.day"=celfile.timestamp[ ,"file.day"], "file.hour"=celfile.timestamp[ ,"file.hour"], "batch"=NA, "cellid"=sampleinfo[ , "Source.Name"], sampleinfo)
  sampleinfo2 <- sampleinfo
  ## remove duplcated cell line hybridization
  dupln <- duplicated(sampleinfo[ ,"cellid"])
  if(sum(dupln) > 1) { warning("some cell lines have been profiled for gene expression multiple times, only the first instance is kept") }
  sampleinfo <- sampleinfo[!dupln, , drop=FALSE]
  rownames(sampleinfo) <- sampleinfo[ ,"cellid"]

  ## update of cgp cell line collection
  celline.cgp <- data.frame("cellid"=as.character(celline.cgp[ , "Sample.name"]), celline.cgp)
  celline.cgp[ , "cellid"] <- as.character(celline.cgp[ , "cellid"])
  ## add url based on COSMIC IDs
  uurl <- paste("http://cancer.sanger.ac.uk/cell_lines/sample/overview?id=", celline.cgp[ , "ID_sample"], sep="")
  uurl[is.na(celline.cgp[ , "ID_sample"])] <- NA
  celline.cgp <- data.frame("cellid"=celline.cgp[ , "cellid"], "link"=uurl, celline.cgp[ , !is.element(colnames(celline.cgp), "cellid")])

  ## drugpheno
  cellnall <- sort(unique(c(as.character(sampleinfo[ ,"cellid"]), as.character(drugpheno[ ,"cellid"]))))
  dd <- data.frame(matrix(NA, ncol=ncol(drugpheno), nrow=length(cellnall), dimnames=list(cellnall, colnames(drugpheno))))
  newlev <- sapply(drugpheno, levels)
  newlev$cellid <- cellnall
  dd <- genefu::setcolclass.df(df=dd, colclass=sapply(drugpheno, class), factor.levels=newlev)
  dd[rownames(drugpheno),colnames(drugpheno)] <- drugpheno
  dd[ ,"cellid"] <- cellnall
  drugpheno <- dd
  
  ## mutation
  dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
  dd[rownames(mutation), colnames(mutation)] <- mutation
  rownames(dd) <- cellnall
  mutation <- dd

  ## reproducibility between different screening sites
  ## camptothecin was screened at MGH (drug id 195) and WTSI (drug id 1003)
  ## data only available in the supplementary infomration of the Nature website
  myfn2 <- file.path(saveres, "nature_supplinfo_drugpheno_cgp.RData")
  if(!file.exists(myfn2)) {
    drugpheno.nature <- gdata::read.xls(xls=file.path(path.drug, "nature_supplementary_information.xls"), sheet=2)
    drugpheno.nature[drugpheno.nature == "" | drugpheno.nature == " "] <- NA
    save(list="drugpheno.nature", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
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
  ## camptothecin
  pdf(file.path(saveres, "cgp_camptothecin_mgh_wtsi_paper.pdf"))
  yy <- -log10(ic50[ , "drugid_195", drop=FALSE])
  xx <- -log10(ic50[ , "drugid_1003", drop=FALSE])
  ccix <- complete.cases(xx, yy)
  nnn <- sum(ccix)
  cc <- cor.test(x=xx, y=yy, method="spearman", use="complete.obs", alternative="greater")
  cci <- spearmanCI(x=cc$estimate, n=sum(ccix))
  par(mar=c(4, 4, 3, 1) + 0.1)
  llim <- round(range(c(xx, yy), na.rm=TRUE) * 10) / 10
  myScatterPlot(x=xx, y=yy, xlab="-log10 IC50 (WTSI)", ylab="-log10 IC50 (MGH)", main="CAMPTOTHECIN", pch=16, method="transparent", transparency=0.75)
  legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2)
  dev.off()

  ## drug information
  message("Read drug information")
  druginfo <- read.csv(file.path(path.drug, "cgp_drug_information.csv"))
  druginfo[!is.na(druginfo) & (druginfo == " " | druginfo == " ")] <- NA
  druginfo <- data.frame("drug.name"=toupper(gsub(badchars, "", druginfo[ ,"Name"])), druginfo)
  myx <- match(druginfo[ , "drug.name"], drugnid[ , "drug.name"])
  if (any(is.na(myx))) { stop ("Some drugs have missing annotations") }
  ## correct ambiguity for AZD6482: drugid_156 corresponds to the first occurence of AZD6482 while drugid_1066 corresponds to the second
  ## table(!is.na(drugpheno[ , "drugid_156_AUC"]))
  ## table(!is.na(drugpheno[ , "drugid_1066_AUC"]))
  myx[druginfo[ , "drug.name"] == "AZD6482"][2] <- which(drugnid[ , "drug.name"] == "AZD6482")[2]
  druginfo <- data.frame("drugid"=rownames(drugnid)[myx], drugnid[myx, , drop=FALSE], druginfo)
  rownames(druginfo) <- as.character(druginfo[ ,"drugid"])
  ## complement drug infomration with the supplementary infomration from the Nature website
  myfn2 <- file.path(saveres, "nature_supplinfo_druginfo_cgp.RData")
  if(!file.exists(myfn2)) {
    druginfo.nature <- gdata::read.xls(xls=file.path(path.drug, "nature_supplementary_information.xls"), sheet=4)
    druginfo.nature[druginfo.nature == "" | druginfo.nature == " "] <- NA
    save(list="druginfo.nature", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  rownames(druginfo.nature) <- paste("drugid", druginfo.nature[ , "Drug.ID"], sep="_")
  druginfo <- data.frame(druginfo, druginfo.nature[rownames(druginfo), c("Brand.name", "Site.of.screening", "Drug.type", "Drug.class.I", "Drug.class.II", "Target.family", "Effector.pathway.biological.process", "Clinical.trials", "Source")])

  ## drug concentration
  message("Read drug concentration")
  drugconc <- read.csv(file.path(path.drug, "cgp_drug_concentration.csv"))
  drugconc[!is.na(drugconc) & (drugconc == "" | drugconc == " ")] <- NA
  drugconc <- data.frame("drug.name"=toupper(gsub(badchars, "", drugconc[ ,"Compound.Name"])), drugconc)
  if(all(!is.element(drugconc[ , "drug.name"], drugnid[ , "drug.name"]))) { stop("Screening concentration for drugs without identifiers!") }
  myx <- match(drugconc[ , "drug.name"], drugnid[ , "drug.name"])
  ## correct ambiguity for AZD6482: drugid_156 corresponds to the first occurence of AZD6482 while drugid_1066 corresponds to the second
  myx[drugconc[ , "drug.name"] == "AZD6482"][2] <- which(drugnid[ , "drug.name"] == "AZD6482")[2]
  rownames(drugconc) <- rownames(drugnid)[myx]
  drugconc <- data.frame("drugid"=rownames(drugconc), drugconc)

  ## combine all drugs
  dix <- sort(unique(c(rownames(druginfo), rownames(drugconc), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
  ## update druginfo
  druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
  newlev <- sapply(druginfo, levels)
  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  druginfo2 <- genefu::setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
  druginfo2[match(rownames(druginfo), dix), colnames(druginfo)] <- druginfo
  druginfo2[ , "drugid"] <- newlev$drugid
  druginfo <- druginfo2
  ## update drugconc
  drugconc2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(drugconc), dimnames=list(dix, colnames(drugconc))))
  newlev <- sapply(drugconc, levels)
  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  drugconc2 <- genefu::setcolclass.df(df=drugconc2, colclass=sapply(drugconc, class), factor.levels=newlev)
  drugconc2[match(rownames(drugconc), dix), colnames(drugconc)] <- drugconc
  drugconc2[ , "drugid"] <- newlev$drugid
  drugconc <- drugconc2

  ## report concentrations per cell line and per drug
  drugconc2 <- data.frame(matrix(NA, nrow=nrow(drugconc) * length(cellnall), ncol=6, dimnames=list(paste(rep(rownames(drugconc), times=length(cellnall)), rep(cellnall, each=nrow(drugconc)), sep="..."), c("cellid", "drugid", "drug.name", "nbr.conc.tested", "min.Dose.uM", "max.Dose.uM"))))
  drugconc2[ , "cellid"] <- rep(cellnall, times=nrow(drugconc))
  drugconc2[ , "drugid"] <- rep(rownames(drugconc), each=length(cellnall))
  drugconc2[ , "drug.name"] <- rep(as.character(drugconc[ ,"drug.name"]), each=length(cellnall))
  ## as mentioned in the supplementary information of Garnett et al., a single cell line is used on each plate and treated with 28 different drugs over a 9-pt, 256-fold concentration range
  drugconc2[ , "nbr.conc.tested"] <- 9
  drugconc2[ , "min.Dose.uM"] <- rep(drugconc[ , "Min.Concentration.micromolar."], each=length(cellnall))
  drugconc2[ , "max.Dose.uM"] <- rep(drugconc[ , "Max.Concentration.micromolar."], each=length(cellnall))
  drugconc <- drugconc2
  
  ## normalization of gene expression data
  message("Normalize gene expression data")
  myfn2 <- file.path(saveres, "cgp_ge_norm.RData")
  if(!file.exists(myfn2)) {
    # rr <- just.rma(filenames=celfn)
    ## frma normalization using parallel
    genexprs <- NULL
    ss <- parallel::splitIndices(nx=length(celfn), ncl=ceiling(length(celfn) / max.celfiles))
    for (ii in 1:length(ss)) {
      celfnt <- celfn[ss[[ii]]]
      splitix <- parallel::splitIndices(nx=length(celfnt), ncl=nbcore)
      splitix <- splitix[sapply(splitix, length) > 0]
      res <- parallel::mclapply(splitix, function(x, celfn) {
        ## fRMA
        tt <- celfn[x]
        names(tt) <- NULL
        abatch <- affy::read.affybatch(filenames=tt)
        rr <- frma::frma(object=abatch, summarize="robust_weighted_average", verbose=TRUE, input.vecs=hthgu133afrmavecs)
        return(exprs(rr))
      }, celfn=celfnt)
      genexprs <- rbind(genexprs, t(do.call(cbind, res)))
    }
    ## match the experiment labels
    rownames(genexprs) <- rownames(sampleinfo2)[match(rownames(genexprs), as.character(sampleinfo2[ ,"filename"]))]
    save(list=c("genexprs"), compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  
  ## build annotation matrix
  message("Build annotation matrix")
  myfn2 <- file.path(saveres, "cgp_ge_annot.RData")
  if(!file.exists(myfn2)) {
    ## select the best probe for a single gene
    require(jetset) || stop("Library jetset is not available!")
    js <- jetset::jscores(chip="hgu133a", probeset=colnames(genexprs))
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
  
  ## save the full dataset, with duplicates
  sampleinfo.full.cgp <- sampleinfo2
  data.ge.full.cgp <- genexprs
  rownames(data.ge.full.cgp) <- gsub(".CEL.gz", "", rownames(data.ge.full.cgp))
  data.ge.full.cgp <- data.ge.full.cgp[match(rownames(sampleinfo.full.cgp), rownames(data.ge.full.cgp)), , drop=FALSE]
  annot.ge.full.cgp <- annot
  save(list=c("data.ge.full.cgp", "annot.ge.full.cgp", "sampleinfo.full.cgp"), compress=TRUE, file=file.path(saveres, "cgp_full_frma.RData"))

  ## match the experiment labels
  myx <- rownames(sampleinfo)[match(rownames(genexprs), gsub(".CEL.gz", "", as.character(sampleinfo[ ,"filename"])))]
  genexprs <- genexprs[!is.na(myx), , drop=FALSE]
  myx <- myx[!is.na(myx)]
  rownames(genexprs) <- myx

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
  iix <- grep("_IC_50$", colnames(drugpheno))
  iixn <- gsub("_IC_50$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.ic50 <- exp(dd)
  ## IC25 in microM
  iix <- grep("_IC_25$", colnames(drugpheno))
  iixn <- gsub("_IC_25$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.ic25 <- exp(dd)
  ## IC75 in microM
  iix <- grep("_IC_75$", colnames(drugpheno))
  iixn <- gsub("_IC_75$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.ic75 <- exp(dd)
  ## IC90 in microM
  iix <- grep("_IC_90$", colnames(drugpheno))
  iixn <- gsub("_IC_90$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.ic90 <- exp(dd)
  ## AUC (higher value represents sensitivity)
  iix <- grep("_AUC$", colnames(drugpheno))
  iixn <- gsub("_AUC$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  drugpheno.auc <- 1 - dd
  ## slope
  iix <- grep("_BETA$", colnames(drugpheno))
  iixn <- gsub("_BETA$", "", colnames(drugpheno)[iix])
  dd <- matrix(NA, ncol=length(drugnall), nrow=length(cellnall), dimnames=list(cellnall, drugnall))
  dd[rownames(drugpheno), iixn] <- data.matrix(drugpheno[ , iix])
  dd[dd < 0] <- 0
  dd[dd > 1] <- 1
  drugpheno.slope <- dd
  
  ## compute drug sensitivity calls using waterfall plot
  ## IC50
  ic50.call2 <- ic50.call3 <- data.frame(matrix(NA, ncol=ncol(drugpheno.ic50), nrow=nrow(drugpheno.ic50), dimnames=dimnames(drugpheno.ic50)))
  ## IC50 sensitivity calling 2 levels
  pdf(file.path(saveres, "cgp_ic50_sensitivity_calling2_drugs.pdf"), width=10, height=10)
  for(i in 1:ncol(drugpheno.ic50)) {
    ic50.call2[ , i] <- PharmacoGx::callingWaterfall(x=drugpheno.ic50[ , i] / 10^6, type="IC50", intermediate.fold=0, cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CGP)", colnames(drugpheno.ic50)[i]))
  }
  dev.off()
  ## IC50 sensitivity calling 3 levels
  pdf(file.path(saveres, "cgp_ic50_sensitivity_calling3_drugs.pdf"), width=10, height=10)
  for(i in 1:ncol(drugpheno.ic50)) {
    ic50.call3[ , i] <- PharmacoGx::callingWaterfall(x=drugpheno.ic50[ , i] / 10^6, type="IC50", intermediate.fold=4, cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CGP)", colnames(drugpheno.ic50)[i]))
  }
  dev.off()
  dimnames(ic50.call2) <- dimnames(ic50.call3) <- dimnames(drugpheno.ic50)
  ## AUC
  auc.call2 <- auc.call3 <- data.frame(matrix(NA, ncol=ncol(drugpheno.auc), nrow=nrow(drugpheno.auc), dimnames=dimnames(drugpheno.auc)))
  ## AUC sensitivity calling 2 levels
  pdf(file.path(saveres, "cgp_auc_sensitivity_calling2_drugs.pdf"), width=10, height=10)
  for(i in 1:ncol(drugpheno.auc)) {
    auc.call2[ , i] <- PharmacoGx::callingWaterfall(x=drugpheno.auc[ , i] / 10^6, type="AUC", intermediate.fold=0, cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CGP)", colnames(drugpheno.auc)[i]))
  }
  dev.off()
  ## AUC sensitivity calling 3 levels
  pdf(file.path(saveres, "cgp_auc_sensitivity_calling3_drugs.pdf"), width=10, height=10)
  for(i in 1:ncol(drugpheno.auc)) {
    auc.call3[ , i] <- PharmacoGx::callingWaterfall(x=drugpheno.auc[ , i] / 10^6, type="AUC", intermediate.fold=1.2, cor.min.linear=0.95, plot=TRUE, name=sprintf("%s (CGP)", colnames(drugpheno.auc)[i]))
  }
  dev.off()
  dimnames(auc.call2) <- dimnames(auc.call3) <- dimnames(drugpheno.auc)

  ## save drug phenotypes
  drugpheno.orig <- drugpheno
  drugpheno <- list("AUC"=drugpheno.auc, "IC25"=drugpheno.ic25, "IC50"=drugpheno.ic50, "IC75"=drugpheno.ic75, "IC90"=drugpheno.ic90, "SLOPE"=drugpheno.slope, "IC50.CALL2"=ic50.call2, "IC50.CALL3"=ic50.call3, "AUC.CALL2"=auc.call2, "AUC.CALL3"=auc.call3)
  
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
  
  data.ge.cgp <- genexprs
  annot.ge.cgp <- annot
  sampleinfo.cgp <- sampleinfo
  drugpheno.cgp <- drugpheno
  druginfo.cgp <- druginfo
  drugconc.cgp <- drugconc
  mutation.cgp <- mutation

  ## make sure that cellid are not factors
  celline.cgp[, "cellid"] <- as.character(celline.cgp[, "cellid"])
  sampleinfo.cgp[, "cellid"] <- as.character(sampleinfo.cgp[, "cellid"])
  drugconc.cgp[, "cellid"] <- as.character(drugconc.cgp[, "cellid"])

  message("Save data")
  write.csv(celline.cgp, file=file.path(saveres, "cell_line_collection_cgp.csv"), row.names=FALSE)
  write.csv(annot.ge.cgp, file=file.path(saveres, "annot_ge_cgp.csv"), row.names=FALSE)
  write.csv(t(data.ge.cgp), file=file.path(saveres, "data_ge_cgp.csv"))
  write.csv(drugpheno.orig, file=file.path(saveres, "drugpheno_cgp.csv"), row.names=FALSE)
  write.csv(drugpheno.cgp$IC50, file=file.path(saveres, "drugpheno_ic50_cgp.csv"), row.names=TRUE)
  write.csv(drugpheno.cgp$IC25, file=file.path(saveres, "drugpheno_ic25_cgp.csv"), row.names=TRUE)
  write.csv(drugpheno.cgp$IC75, file=file.path(saveres, "drugpheno_ic75_cgp.csv"), row.names=TRUE)
  write.csv(drugpheno.cgp$IC90, file=file.path(saveres, "drugpheno_ic90_cgp.csv"), row.names=TRUE)
  write.csv(drugpheno.cgp$AUC, file=file.path(saveres, "drugpheno_auc_cgp.csv"), row.names=TRUE)
  write.csv(drugpheno.cgp$SLOPE, file=file.path(saveres, "drugpheno_slope_cgp.csv"), row.names=TRUE)
  write.csv(druginfo.cgp, file=file.path(saveres, "druginfo_cgp.csv"), row.names=FALSE)
  write.csv(drugconc.cgp, file=file.path(saveres, "drugconc_cgp.csv"), row.names=FALSE)
  write.csv(sampleinfo.cgp, file=file.path(saveres, "sampleinfo_cgp.csv"), row.names=FALSE)
  write.csv(mutation.cgp, file=file.path(saveres, "mutation_cgp.csv"), row.names=TRUE)
  save(list=c("data.ge.cgp", "annot.ge.cgp", "sampleinfo.cgp", "mutation.cgp", "drugpheno.cgp", "druginfo.cgp", "drugconc.cgp", "celline.cgp"), compress=TRUE, file=myfn)
}





## end




