########################
## Zhaleh Safikhani, Benjamin Haibe-Kains
## Code under License Artistic-2.0
## March 6, 2015
########################

require(stringr) || stop("Library stringr is not available!")
require(magicaxis) || stop("Library magicaxis is not available!")
require(gdata) || stop("Library gdata is not available!")
require(RColorBrewer) || stop("Library RColorBrewer is not available!")


## additional functions
source(file.path("code", "cdrug2_foo.R"))

path.data <- file.path("data")
if(!file.exists(path.data)) { dir.create(path.data, showWarnings=FALSE, recursive=TRUE) }
path.result <- file.path("result")
if(!file.exists(path.result)) { dir.create(path.result, showWarnings=FALSE, recursive=TRUE) }

## download ccle raw sensitivity data 
## all the analysis has been done using previous version: "CCLE_NP24.2009_Drug_data_2012.02.20.csv" which is not availbale for download anymore
## if you need previous version please contact benjamin.haibe.kains@utoronto.ca

archivn <- "CCLE_NP24.2009_Drug_data_2015.02.24.csv"
if(!file.exists(file.path(path.data, "dwl"))){dir.create(file.path(path.data, "dwl"), showWarnings=FALSE)}
myfn <- file.path(path.data, "dwl", archivn)
if (!file.exists(myfn)) {
  dwl.status <- download.file(url="http://www.broadinstitute.org/ccle/downloadFile/DefaultSystemRoot/exp_10/ds_27/CCLE_NP24.2009_Drug_data_2015.02.24.csv?downloadff=true&fileId=20777", destfile=file.path(path.data, "dwl", archivn), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
}

ccle.raw.drug.sensitivity <- read.csv(myfn, stringsAsFactors = FALSE)
ccle.raw.drug.sensitivity.list <- do.call(c, apply(ccle.raw.drug.sensitivity, 1, list))

## create the ccle.drug.response object including information viablilities and concentrations for each cell/drug pair
concentrations.no <- max(unlist(lapply(ccle.raw.drug.sensitivity$Doses..uM., function(x){length(unlist(strsplit(x, split = ",")))})))
obj <- array(NA, dim = c(length(unique(ccle.raw.drug.sensitivity$Primary.Cell.Line.Name)), length(unique(ccle.raw.drug.sensitivity$Compound)), 2, concentrations.no), dimnames = list(unique(ccle.raw.drug.sensitivity$Primary.Cell.Line.Name), unique(ccle.raw.drug.sensitivity$Compound), c("concentration","viability"), 1:concentrations.no))


ccle.raw.drug.sensitivity.list <- do.call(c, apply(ccle.raw.drug.sensitivity, 1, list))
ccle.raw.drug.sensitivity.res <- mapply(fnStandardizeSensitivity, values = ccle.raw.drug.sensitivity.list, method = "ccle")
ccle.drug.response <- obj

## create an ccle.AUC object including AUCs for each cell/drug pair
ccle.AUC <- matrix(NA, nrow = length(unique(ccle.raw.drug.sensitivity$Primary.Cell.Line.Name)), ncol = length(unique(ccle.raw.drug.sensitivity$Compound)))
rownames(ccle.AUC) <- unique(ccle.raw.drug.sensitivity$Primary.Cell.Line.Name)
colnames(ccle.AUC) <- unique(ccle.raw.drug.sensitivity$Compound)
test <- lapply(ccle.raw.drug.sensitivity.list, function(x){
	ccle.AUC[x["Primary.Cell.Line.Name"],x["Compound"]] <<- (as.numeric(x["ActArea"])/concentrations.no)
	})

## create an ccle.IC50 object including IC50s for each cell/drug pair
ccle.IC50 <- matrix(NA, nrow = length(unique(ccle.raw.drug.sensitivity$Primary.Cell.Line.Name)), ncol = length(unique(ccle.raw.drug.sensitivity$Compound)))
rownames(ccle.IC50) <- unique(ccle.raw.drug.sensitivity$Primary.Cell.Line.Name)
colnames(ccle.IC50) <- unique(ccle.raw.drug.sensitivity$Compound)
test <- lapply(ccle.raw.drug.sensitivity.list, function(x){
	ccle.IC50[x["Primary.Cell.Line.Name"],x["Compound"]] <<- as.numeric(x["IC50..uM."])
	})
ccle.IC50 <- -log10(exp(as.matrix(ccle.IC50))/10^6)

## download cgp raw sensitivity data
archivn <- "gdsc_drug_sensitivity_raw_data_w5"
if(!file.exists(file.path(path.data, "dwl"))){dir.create(file.path(path.data, "dwl"), showWarnings=FALSE)}
myfn <- file.path(path.data, "dwl", sprintf("%s.csv", archivn))
if (!file.exists(myfn)) {
  dwl.status <- download.file(url=sprintf("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/%s", sprintf("%s.zip", archivn)), destfile=file.path(path.data, "dwl", sprintf("%s.zip", archivn)), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
}
res <- unzip(zipfile=file.path(path.data, "dwl", sprintf("%s.zip", archivn)),exdir = file.path(path.data, "dwl"))
cgp.raw.drug.sensitivity <- read.csv(myfn, stringsAsFactors = FALSE)

## create the cgp.drug.response object including information viablilities and concentrations for each cell/drug pair
concentrations.no <- length(grep("raw", names(cgp.raw.drug.sensitivity)))
obj <- array(NA,  dim = c(length(unique(cgp.raw.drug.sensitivity$CELL_LINE_NAME)),  length(unique(cgp.raw.drug.sensitivity$DRUG_ID)),  2,  concentrations.no), dimnames = list(unique(cgp.raw.drug.sensitivity$CELL_LINE_NAME), unique(cgp.raw.drug.sensitivity$DRUG_ID), c("concentration","viability"), 1:concentrations.no))


cgp.raw.drug.sensitivity.list <- do.call(c, apply(cgp.raw.drug.sensitivity, 1, list))
cgp.raw.drug.sensitivity.res <- mapply(fnStandardizeSensitivity, values = cgp.raw.drug.sensitivity.list, method = "cgp")
cgp.drug.response <- obj

## download cgp drug sensitivity data
archivn <- "gdsc_manova_input_w5"
if(!file.exists(file.path(path.data, "dwl"))){dir.create(file.path(path.data, "dwl"), showWarnings=FALSE)}
myfn <- file.path(path.data, "dwl", sprintf("%s.csv", archivn))
if (!file.exists(myfn)) {
  dwl.status <- download.file(url=sprintf("ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/%s", sprintf("%s.csv", archivn)), destfile=file.path(path.data, "dwl", sprintf("%s.csv", archivn)), quiet=TRUE)
  if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
}
cgp.manova <- read.csv(myfn, stringsAsFactors = FALSE)
indices <- grep("IC_50$",colnames(cgp.manova))
colnames(cgp.manova)[indices] <- cgp.manova[1,indices]
indices <- grep("AUC$",colnames(cgp.manova))
colnames(cgp.manova)[indices] <- cgp.manova[1,indices]
cgp.manova <- cgp.manova[6:nrow(cgp.manova),]
rownames(cgp.manova)  <- cgp.manova$Cell.Line

## create an cgp.AUC object including AUCs for each cell/drug pair
cgp.AUC <- matrix(NA, nrow = length(unique(cgp.raw.drug.sensitivity$CELL_LINE_NAME)), ncol = length(unique(cgp.raw.drug.sensitivity$DRUG_ID)))
rownames(cgp.AUC) <- unique(cgp.raw.drug.sensitivity$CELL_LINE_NAME)
colnames(cgp.AUC) <- unique(cgp.raw.drug.sensitivity$DRUG_ID)
test <- lapply(cgp.raw.drug.sensitivity.list, function(x){ drug <- stringr::str_trim(x["DRUG_ID"]); y <- cgp.manova[x["CELL_LINE_NAME"], sprintf("%s_AUC", drug)]; if(length(y) > 0) {cgp.AUC[x["CELL_LINE_NAME"],drug] <<- 1 - as.numeric(y)}})


## create an cgp.IC50 object including IC50s for each cell/drug pair
cgp.IC50 <- matrix(NA, nrow = length(unique(cgp.raw.drug.sensitivity$CELL_LINE_NAME)), ncol = length(unique(cgp.raw.drug.sensitivity$DRUG_ID)))
rownames(cgp.IC50) <- unique(cgp.raw.drug.sensitivity$CELL_LINE_NAME)
colnames(cgp.IC50) <- unique(cgp.raw.drug.sensitivity$DRUG_ID)
test <- lapply(cgp.raw.drug.sensitivity.list, function(x){drug <- stringr::str_trim(x["DRUG_ID"]); y <- cgp.manova[x["CELL_LINE_NAME"], sprintf("%s_IC_50", drug)]; if(length(y) > 0) {cgp.IC50[x["CELL_LINE_NAME"],drug] <<- as.numeric(y)}})
###Original IC50 is in ln scale
cgp.IC50 <- -log10(exp(as.matrix(cgp.IC50))/10^6)


## drug mapping
drug.map <- read.csv(file.path(path.data,"matching_drug_CCLE_CGP.csv"), row.names=1, stringsAsFactors = FALSE)

## cell line mapping
cell.map <- read.csv(file.path(path.data,"matching_cell_line_CCLE_CGP.csv"), row.names=1, stringsAsFactors = FALSE)

## intersection between CCLE and CGP
cix <- rownames(cell.map)[is.element(cell.map[ , "CCLE"], rownames(ccle.AUC)) & is.element(cell.map[ , "CGP"], rownames(cgp.AUC))]
dix <- rownames(drug.map)[is.element(drug.map[ , "CCLE"], colnames(ccle.AUC)) & is.element(drug.map[ , "CGP"], colnames(cgp.AUC))]


## computation of slope and slope*(ST) for all the common experiments between ccle and cgp
cgp.slope <- cgp.slope.star <- ccle.slope <- ccle.slope.star <- matrix(NA, ncol = length(dix), nrow = length(cix))
dimnames(cgp.slope) <- dimnames(cgp.slope.star) <- dimnames(ccle.slope) <- dimnames(ccle.slope.star) <- list(cix, dix)

for(i in cix) {
  x <- cell.map[i,]
  for(j in dix) {
    y <- drug.map[j,]
    if(as.character(x["CCLE"]) %in% dimnames(ccle.drug.response)[[1]] &
         as.character(x["CGP"]) %in% dimnames(cgp.drug.response)[[1]]) {
      ccle.range <- ccle.drug.response[as.character(x["CCLE"]),as.character(y["CCLE"]),"concentration",]
      cgp.range <- cgp.drug.response[as.character(x["CGP"]),as.character(y["CGP"]),"concentration",]
      
      if((!is.na(table(is.na(cgp.range))["FALSE"])) & (!is.na(table(is.na(ccle.range))["FALSE"]))) {
        common.ranges <- fnCommonConcentrationRange(list(ccle.range,cgp.range))
        ccle.slope.star[i,j] <- fnComputeSlope(doses = common.ranges[[1]], responses = ccle.drug.response[as.character(x["CCLE"]),as.character(y["CCLE"]),"viability",1:length(common.ranges[[1]])])
        ccle.slope[i,j] <- fnComputeSlope(doses = ccle.range, responses = ccle.drug.response[as.character(x["CCLE"]),as.character(y["CCLE"]),"viability",])
        cgp.slope.star[i,j] <- fnComputeSlope(doses = common.ranges[[2]], responses = cgp.drug.response[as.character(x["CGP"]),as.character(y["CGP"]),"viability",1:length(common.ranges[[2]])])
        cgp.slope[i,j] <- fnComputeSlope(doses = cgp.range, responses = cgp.drug.response[as.character(x["CGP"]),as.character(y["CGP"]),"viability",])
      }
    }
  }
}

## arbitrary cutoff of slope = -16
cutoff.slope <- -16

## slope*(ST) drug sensitivity calling
cgp.slope.star.sensitivity.call <- cgp.slope.star <= cutoff.slope
cgp.slope.star.sensitivity.call <- data.frame(t(apply(cgp.slope.star.sensitivity.call, 1, function (x) {
  x <- as.factor(x)
  levels(x) <- c("insensitive", "sensitive")
  return (x)
})), stringsAsFactors=TRUE)

ccle.slope.star.sensitivity.call <- ccle.slope.star <= cutoff.slope
ccle.slope.star.sensitivity.call <- data.frame(t(apply(ccle.slope.star.sensitivity.call, 1, function (x) {
  x <- as.factor(x)
  levels(x) <- c("insensitive", "sensitive")
  return (x)
})), stringsAsFactors=TRUE)

## slope drug sensitivity calling
cgp.slope.sensitivity.call <- cgp.slope <= cutoff.slope
cgp.slope.sensitivity.call <- data.frame(t(apply(cgp.slope.sensitivity.call, 1, function (x) {
  x <- as.factor(x)
  levels(x) <- c("insensitive", "sensitive")
  return (x)
})), stringsAsFactors=TRUE)

ccle.slope.sensitivity.call <- ccle.slope <= cutoff.slope
ccle.slope.sensitivity.call <- data.frame(t(apply(ccle.slope.sensitivity.call, 1, function (x) {
  x <- as.factor(x)
  levels(x) <- c("insensitive", "sensitive")
  return (x)
})), stringsAsFactors=TRUE)


## arbitrary cutoff of AUC = 0.2
cutoff.AUC <- 0.2

### extract auc for all the common experiments between ccle and cgp
ccle2.AUC <- data.matrix(ccle.AUC[cell.map[cix, "CCLE"], drug.map[dix, "CCLE"]])
cgp2.AUC <- data.matrix(cgp.AUC[cell.map[cix, "CGP"], as.character(drug.map[dix, "CGP"])])
dimnames(cgp2.AUC) <- dimnames(ccle2.AUC) <- list(cix, dix)
nix <- is.na(cgp2.AUC) | is.na(ccle2.AUC)
cgp2.AUC[nix] <- ccle2.AUC[nix] <- NA


### auc drug sensitivity calling
cgp.auc.sensitivity.call <- cgp2.AUC >= cutoff.AUC
cgp.auc.sensitivity.call <- data.frame(t(apply(cgp.auc.sensitivity.call, 1, function (x) {
  x <- as.factor(x)
  levels(x) <- c("insensitive", "sensitive")
  return (x)
})), stringsAsFactors=TRUE)

ccle.auc.sensitivity.call <- ccle2.AUC >= cutoff.AUC
ccle.auc.sensitivity.call <- data.frame(t(apply(ccle.auc.sensitivity.call, 1, function (x) {
  x <- as.factor(x)
  levels(x) <- c("insensitive", "sensitive")
  return (x)
})), stringsAsFactors=TRUE)

## download bouhaddou et al manual curation file
## contact M Birtwistle (Mount Sinai Hospital, NY)
## myfn <- file.path("Manual_Curation_All.xlsx")
manual.curation.file <- gdata::read.xls(xls=myfn, sheet=1, stringsAsFactors = FALSE)

is.ccle.sensitive <- function(x) x == 1 | x == 3
is.cgp.sensitive <- function(x) x == 1 | x == 4

### manual curation sensitivity calling
mc.cgp.sensitivity.call <-  mc.ccle.sensitivity.call <- list()
for(i in c("EC","TT","ER"))
{
  mc.ccle.sensitivity.call[[i]] <- matrix(NA, ncol = length(unique(manual.curation.file$Drug.ID.number)), nrow = length(unique(manual.curation.file$Cell.Line.Name)))
  colnames(mc.ccle.sensitivity.call[[i]]) <- unique(manual.curation.file$Drug.ID.number); rownames(mc.ccle.sensitivity.call[[i]]) <- unique(manual.curation.file$Cell.Line.Name)
  mc.cgp.sensitivity.call[[i]] <- mc.ccle.sensitivity.call[[i]]
  
  x <- manual.curation.file[,i]
  mc.ccle <- sapply(1:length(x), function(idx){
    t <- ifelse(is.ccle.sensitive(x[idx]), "sensitive", "insensitive")
  })
  
  mc.cgp <- sapply(1:length(x), function(idx){
    t <- ifelse(is.cgp.sensitive(x[idx]), "sensitive", "insensitive")
  })
  mc.ccle <- cbind("drug" = manual.curation.file$Drug.ID.number, "cell" = manual.curation.file$Cell.Line.Name, "value" = mc.ccle)
  mc.cgp <- cbind("drug" = manual.curation.file$Drug.ID.number, "cell" = manual.curation.file$Cell.Line.Name, "value" = mc.cgp)
  
  test <- sapply(1:nrow(mc.ccle), function(idx){mc.ccle.sensitivity.call[[i]][mc.ccle[idx,"cell"],mc.ccle[idx,"drug"]] <<- mc.ccle[idx,"value"]})
  test <- sapply(1:nrow(mc.cgp), function(idx){mc.cgp.sensitivity.call[[i]][mc.cgp[idx,"cell"],mc.cgp[idx,"drug"]] <<- mc.cgp[idx,"value"]})
  colnames(mc.cgp.sensitivity.call[[i]]) <- colnames(mc.ccle.sensitivity.call[[i]]) <- unlist(lapply( colnames(mc.ccle.sensitivity.call[[i]]), function(x){rownames(drug.map)[which(drug.map$bouhaddou == as.character(x))]}))
}

require(epibasix) || stop("Library epibasix is not available!")
## AUC
## contingency tables
tt <- table("CGP"=unlist(cgp.auc.sensitivity.call), "CCLE"=unlist(ccle.auc.sensitivity.call))
print(tt/sum(tt))
tt <- list(tt)
for(i in 1:ncol(cgp.auc.sensitivity.call)) {
  tt <- c(tt, list(table("CGP"=cgp.auc.sensitivity.call[ , i], "CCLE"=ccle.auc.sensitivity.call[ , i])))
}
names(tt) <- c("All", colnames(cgp.auc.sensitivity.call))
## kappa
kapp <- sapply(tt, function (x) {
  err <- try(rr <- epibasix::epiKappa(x, k0=0), silent=TRUE)
  if(class(err) == "try-error") {
    rr <- NA
  } else {
    rr <- rr$kappa
  }
  return (rr)
})
auc.kapp <- kapp

## Slope*(ST)
## contingency tables
tt <- table("CGP"=unlist(cgp.slope.star.sensitivity.call), "CCLE"=unlist(ccle.slope.star.sensitivity.call))
print(tt/sum(tt))
tt <- list(tt)
for(i in 1:ncol(cgp.slope.star.sensitivity.call)) {
  tt <- c(tt, list(table("CGP"=cgp.slope.star.sensitivity.call[ , i], "CCLE"=ccle.slope.star.sensitivity.call[ , i])))
}
names(tt) <- c("All", colnames(cgp.slope.star.sensitivity.call))
## kappa
kapp <- sapply(tt, function (x) {
  err <- try(rr <- epibasix::epiKappa(x, k0=0), silent=TRUE)
  if(class(err) == "try-error") {
    rr <- NA
  } else {
    rr <- rr$kappa
  }
  return (rr)
})
slope.star.kapp <- kapp

## Slope
## contingency tables
tt <- table("CGP"=unlist(cgp.slope.sensitivity.call), "CCLE"=unlist(ccle.slope.sensitivity.call))
print(tt/sum(tt))
tt <- list(tt)
for(i in 1:ncol(cgp.slope.sensitivity.call)) {
  tt <- c(tt, list(table("CGP"=cgp.slope.sensitivity.call[ , i], "CCLE"=ccle.slope.sensitivity.call[ , i])))
}
names(tt) <- c("All", colnames(cgp.slope.sensitivity.call))
## kappa
kapp <- sapply(tt, function (x) {
  err <- try(rr <- epibasix::epiKappa(x, k0=0), silent=TRUE)
  if(class(err) == "try-error") {
    rr <- NA
  } else {
    rr <- rr$kappa
  }
  return (rr)
})
slope.kapp <- kapp

## MCs
## contingency tables

mc.kapp <- list()
for(mci in c("EC","TT","ER"))
{
  tt <- table("CGP"=unlist(mc.cgp.sensitivity.call[[mci]]), "CCLE"=unlist(mc.ccle.sensitivity.call[[mci]]))
  print(tt/sum(tt))
  tt <- list(tt)
  for(i in 1:ncol(mc.cgp.sensitivity.call[[mci]])) {
    tt <- c(tt, list(table("CGP"=mc.cgp.sensitivity.call[[mci]][ , i], "CCLE"=mc.ccle.sensitivity.call[[mci]][ , i])))
  }
  names(tt) <- c("All", colnames(mc.cgp.sensitivity.call[[mci]]))
  ## kappa
  kapp <- sapply(tt, function (x) {
    err <- try(rr <- epibasix::epiKappa(x, k0=0), silent=TRUE)
    if(class(err) == "try-error") {
      rr <- NA
    } else {
      rr <- rr$kappa
    }
    return (rr)
  })
  mc.kapp[[mci]] <- kapp
}

## Supplementary figure 1 and 2
mycol <- RColorBrewer::brewer.pal(n=6, name="Set3")
pdf(file.path(path.result, "kappa_consistency.pdf"), height=4, width=6)
par(xaxt="n", las=1)

for ( j in names(auc.kapp))
{
  if(j == "X17AAG"){i <- "17AAG"} else {i <- j}
  tt <- cbind(auc.kapp[j], slope.kapp[j], slope.star.kapp[j], mc.kapp[[1]][i], mc.kapp[[2]][i],mc.kapp[[3]][i])
  tt <- t(apply(tt,1, function(x){ifelse(x > 0 , x , 0)}))
  colnames(tt) <- c("AUC", "ST*", "ST", "MC1", "MC2", "MC3")
  mp <- barplot(t(tt), beside = T, space=0.3, col=mycol, ylab="Kappa", ylim=c(0, 1), main = i)
  axis(1, at=seq(1, ncol(tt), by=1), labels=FALSE)
  text(x=mp + (max(mp) * 0.02), y=par("usr")[3] - (par("usr")[4] * 0.05), pos=2, labels=colnames(tt), srt=45, xpd=NA, cex=0.8, font=1)
  abline(h=0.8, lty=2)
}
dev.off()

## Supplementary figue 3 and 4
pdf(file.path(path.result, sprintf("curves.pdf")), height=12, width=12)
par(mfrow = c(3,3))
for(cellline in cix)
{
  for(drug in dix)
  {
    if (!is.na(ccle.slope.star[cellline,drug]) & !is.na(cgp.slope.star[cellline,drug]))
    {
      ccle.avail.doses <- which(!is.na(ccle.drug.response[cell.map[cellline,"CCLE"],drug.map[drug,"CCLE"],"concentration",]))
      cgp.avail.doses <- which(!is.na(cgp.drug.response[cell.map[cellline,"CGP"],as.character(drug.map[drug,"CGP"]),"concentration",]))
      doses <- list("ccle" = ccle.drug.response[cell.map[cellline,"CCLE"],drug.map[drug,"CCLE"],"concentration",ccle.avail.doses], "cgp" = cgp.drug.response[cell.map[cellline,"CGP"],as.character(drug.map[drug,"CGP"]),"concentration",cgp.avail.doses])
      responses <- list("ccle" = ccle.drug.response[cell.map[cellline,"CCLE"],drug.map[drug,"CCLE"],"viability",ccle.avail.doses], "cgp" = cgp.drug.response[cell.map[cellline,"CGP"],as.character(drug.map[drug,"CGP"]),"viability",cgp.avail.doses])
      legend.values <- list("ccle"= ccle.slope.star[cellline,drug], "cgp" = cgp.slope.star[cellline,drug])
      fnPlotDrugResponseCurve(drug, cellline, doses, responses, legend.values)
    }
  } 
}
dev.off()

