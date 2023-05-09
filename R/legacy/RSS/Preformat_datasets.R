## Preparing the public datasets
setwd("/run/media/job/Noir/DATASETS/GEOQUERY")
library(GEOquery)
library(impute)
library(limma)

###############
## FUNCTIONS ##
###############

## Get and/or load GEO data
goGEO <- function(geoid, path=".", prefix="", suffix="", AnnotGPL=F, GSEMatrix=T, raw.include=T) {
  if (prefix != "") prefix <- paste(prefix, "_", sep="")
  if (suffix != "") suffix <- paste("_", suffix, sep="")
  dirname <- paste(path, "/", prefix, geoid, suffix, sep="")
  dir.create(dirname)
  myGEXdata <- getGEO(geoid, destdir=dirname, AnnotGPL=AnnotGPL, GSEMatrix=GSEMatrix)[[1]]
  if ((raw.include) & (length(dir(dirname, paste(geoid, "_RAW.tar", sep="")))==0)) getGEOSuppFiles(geoid, makeDirectory=F, baseDir=dirname)
  if(length(dir(path=dirname, pattern=paste(geoid, "_RAW.tar", sep="")))) save(myGEXdata, file=paste(dirname, "/", geoid, ".Rdata", sep=""))
  return(myGEXdata)
  gc()
}

## Merging mode functions
colMax <- function(X) apply(X, 2, function(x) {max(x, na.rm=T)})
colMin <- function(X) apply(X, 2, function(x) {min(x, na.rm=T)})
colMedian <- function(X) apply(X, 2, function(x) {median(x, na.rm=T)})
colNA <- function(X) apply(X, 2, function(x) {length(which(is.na(x)))})
colVar <- function(X) apply(X, 2, function(x) {var(x, na.rm=T)})
colSd <- function(X) apply(X, 2, function(x) {sd(x, na.rm=T)})
rowMax <- function(X) apply(X, 1, function(x) {max(x, na.rm=T)})
rowMin <- function(X) apply(X, 1, function(x) {min(x, na.rm=T)})
rowMedian <- function(X) apply(X, 1, function(x) {median(x, na.rm=T)})
rowNA <- function(X) apply(X, 1, function(x) {length(which(is.na(x)))})
rowVar <- function(X) apply(X, 1, function(x) {var(x, na.rm=T)})
rowSd <- function(X) apply(X, 1, function(x) {sd(x, na.rm=T)})

## Features merging
feat.merge <- function(fmat, mode="median") {
  usym <- unique(rownames(fmat))
  xmat <- t(sapply(1:length(usym), function(s) {
    symsel <- which(rownames(fmat) == usym[s])
    if (length(symsel) > 1) {
      if (mode == "median") return(as.numeric(apply(fmat[symsel, ], 2, function(x) {median(x, na.rm=T)})))
      if (mode == "min") return(as.numeric(apply(fmat[symsel, ], 2, function(x) {min(x, na.rm=T)})))
      if (mode == "max") return(as.numeric(apply(fmat[symsel, ], 2, function(x) {max(x, na.rm=T)})))
      if (mode == "max") return(as.numeric(apply(fmat[symsel, ], 2, function(x) {sum(x, na.rm=T)})))
      if (mode == "mean") return(as.numeric(colMeans(fmat[symsel, ])))
    } else return(fmat[symsel,])
  }))
  colnames(xmat) <- colnames(fmat)
  rownames(xmat) <- usym
  return(xmat)
}

## LOADING AND PREPROCESSING OF GEO DATASETS
preproc <- function(geoid, prefix="", suffix="", path, feat.mode="median", keep.sample=NULL, filter.nosymbol=T, NA.freqmax=.5, NA.impute=T, qn=T, AnnotGPL=F, GSEMatrix=T, raw.include=T) {
  ## Load data
  cat("Loading data ...\n")
  Xdata <- goGEO(geoid=geoid, prefix=prefix, path=path, AnnotGPL=AnnotGPL, GSEMatrix=GSEMatrix, raw.include=raw.include)
  ## Some basic infos regarding features annotation
  cat("Initial gex matrix dimensions :", dim(Xdata@assayData$exprs), "\n")
  ## Filtering out samples
  if (!is.null(keep.sample)) {
    cat("Restrained samples from", ncol(Xdata@assayData$exprs))
    Xdata <- Xdata[,keep.sample]
    cat(" to", ncol(Xdata@assayData$exprs),"\n")
  }
  ## Defining data format
  gex.type <- "int"
  if(max(Xdata@assayData$exprs, na.rm=T) < 1000) gex.type <- "log2"
  cat("Measure format found :", gex.type, "(", range(Xdata@assayData$exprs, na.rm=T), ")\n")
  ## Finding the Symbol column
  fcn <-colnames(Xdata@featureData@data)
  symb.cn <- fcn[grep("symbol", tolower(fcn))]
  cat("Symbol column found :", symb.cn, "\n")
  ## Checking the amount of unique symbols for the design
  symlen <- length(unique(Xdata@featureData@data[[symb.cn]]))
  cat("Found", symlen, "unique symbols.\n")
  ## Filtering out probes without a symbol
  Xgex <- exprs(Xdata)
  if(filter.nosymbol) {
    cat("Filtering out features without symbol :", nrow(Xgex))
    symsel <- which(Xdata@featureData@data[[symb.cn]] != "")
    Xdata <- Xdata[symsel,]
    Xgex <- Xgex[symsel,]
    cat(" to", nrow(Xgex),"\n")
  }
  ## Shifting data if necessary
  if( (gex.type == "int") & (min(Xgex, na.rm=T) < 1)) Xgex <- Xgex - min(Xgex, na.rm=T) + 1
  if( (gex.type == "log2") & (min(Xgex, na.rm=T) < 0)) Xgex <- Xgex - min(Xgex, na.rm=T)
  ## Filtering out probes with more than XX% of NA in the exprs
  if(NA.freqmax < 1) {
    cat("Filtering out features with more than", paste(NA.freqmax*100, "%", sep=""), "NA : from", nrow(Xgex))
    NAfreq <- rowNA(Xgex)/ncol(Xgex)
    Xgex <- Xgex[which(NAfreq < NA.freqmax),]
    Xdata <- Xdata[which(NAfreq < NA.freqmax),]
    cat(" to", nrow(Xgex),"\n")
  }
  ## Imputing missing values
  if((NA.impute) & (length(which(is.na(Xgex))) > 0)) {
    cat("Imputing (knn) remaining NAs ...\n")
    Xgex <- impute.knn(Xgex)$data
  }
  ## Transforming to log2
  if(gex.type == "int") {
    cat("Transforming data from int to log2 ...\n")
    Xgex <- log2(Xgex)
  }
  ## Merge features
  cat("Merging features to unique symbols ...\n")
  rownames(Xgex) <- Xdata@featureData@data[[symb.cn]]
  Xgex <- feat.merge(Xgex, mode=feat.mode)
  ## "raw" boxplots
  cat("Plotting boxplots ...\n")
  png(paste(path, "/", prefix, "_", geoid, "/", geoid, "_", suffix, "_box.png", sep=""), 1024, 768)
  boxplot(Xgex)
  dev.off()
  ## Dumping pre-QN data
  Xout <- data.frame(Symbol=rownames(Xgex), Xgex, stringsAsFactors=F)
  write.table(Xout, paste(path, "/", prefix, "_", geoid, "/", geoid, "_", suffix, "_gex_norm_", ncol(Xgex), ".", nrow(Xgex),"_l2i.txt", sep=""), sep="\t", quote=F, row.names=F)
  ## Quantiles-norm
  if(qn) {
    cat("Performing quantiles normalization ...\n")
    Xgex <- normalizeQuantiles(Xgex)
    cat("Plotting QN boxplots ...\n")
    png(paste(path, "/", prefix, "_", geoid, "/", geoid, "_", suffix, "_boxQN.png", sep=""), 1024, 768)
    boxplot(Xgex)
    dev.off()
    ## Dumping post-QN data
    Xout <- data.frame(Symbol=rownames(Xgex), Xgex, stringsAsFactors=F)
    write.table(Xout, paste(path, "/", prefix, "_", geoid, "/", geoid, "_", suffix, "_gex_normQN_", ncol(Xgex), ".", nrow(Xgex),"_l2i.txt", sep=""), sep="\t", quote=F, row.names=F)
  }
  ## Hierarchical clustering
  cat("Plotting hierarchical clustering ...\n")
  Xhc <- hclust(dist(t(Xgex), method="euclidean"), method="ward")
  png(paste(path, "/", prefix, "_", geoid, "/", geoid, "_", suffix, "_hclust_eu.wa_allgenes.png", sep=""), 1600, 1000)
  plot(Xhc, main=paste(prefix, " (", nrow(Xgex), ")", sep=""))
  dev.off()
  
  ## Dumping
  write.table(Xdata@phenoData@data, paste(path, "/", prefix, "_", geoid, "/", geoid, "_", suffix, "_annot_", ncol(Xgex), "_samples.annot.txt", sep=""), sep="\t", quote=F, row.names=F)
  
  return(list(exprs=Xgex, pheno=Xdata@phenoData@data))
}


##########
## CORE ##
##########

## MELANOMA
## AUTOMATIC LOADING (once samples to keep are known if dataset is heterogeneous)
## BOGUNOVIC  MM  44  21049   2branches
BOGUNOVIC.MM <- preproc(geoid="GSE19234", prefix="BOGUNOVIC", suffix="MM", path="GEX/MELANOMA/")
BOGUNOVIC.MM.SEL <- preproc(geoid="GSE19234", prefix="BOGUNOVIC", suffix="MM.sel", path="GEX/MELANOMA/", keep.sample=sample(1:44,30))
## JONSSON    MM  57  24614   2(3)branches
JONSSON.MM <- preproc(geoid="GSE22153", prefix="JONSSON", suffix="MM", path="GEX/MELANOMA/")
## MONTOYA    MM  65  21049   3branches
MONTOYA.MM <- preproc(geoid="GSE35640", prefix="MONTOYA", suffix="MM", path="GEX/MELANOMA/")  ## 21049 symbols
## KIEJDA     MM  82  18197   3branches
KIEJDA.MM <- preproc(geoid="GSE29377", prefix="KIEJDA", suffix="MM", path="GEX/MELANOMA/", keep.sample=c(9:90))
## DECECCO    MM  24  18402   2branches
DECECCO.MM <- preproc(geoid="GSE39945", prefix="DECECCO", "MM", path="GEX/MELANOMA/")
## HOLTAN     MM  22  21049   2branches             *WARNING* VERY LOW INTENSITIES!
HOLTAN.MM <- preproc(geoid="GSE23376", prefix="HOLTAN", suffix="MM", path="GEX/MELANOMA/")
## XU         MM  39  13211   2branches (moches)    *WARNING* FEW SYMBOLS!    Rq: Selection of TC>60%
XU.MM <- preproc(geoid="GSE8401", prefix="XU", suffix="MM", path="GEX/MELANOMA/", keep.sample=c(32:39,41,44:47,49,50,52,55,56,59,63:82))  ## 13211 symbols
## HARLIN     MM  44  13211   2(3)branches (moches) *WARNING* FEW SYMBOLS!
HARLIN.MM <- preproc(geoid="GSE12627", prefix="HARLIN", suffix="MM", path="GEX/MELANOMA/", keep.sample=c(1:12,14:18,21:24,29:41,43:52))
## RIKER      MM  40  21049   3branches (unbal)
RIKER.MM <- preproc(geoid="GSE7553", prefix="RIKER", suffix="MM", path="GEX/MELANOMA/", keep=c(9:16,35:40,56:81))
## RASKIN     MM  12  18585   
RASKIN.MM <- preproc(geoid="GSE15605", prefix="RASKIN", suffix="MM", path="GEX/MELANOMA/", keep.sample=c(63:74))
## SCATOLINI  MM  5   18841 (Agilent)
SCATOLINI.MM <- preproc(geoid="GSE12391", prefix="SCATOLINI", suffix="MM", path="GEX/MELANOMA/", keep.sample=c(42:46))
## JEWELL (ca en vaut meme pas la peine, 503 symbols!)

## HAYWARD    MCL 63  21049   2branches
HAYWARD.MCL <- preproc(geoid="GSE7127", prefix="HAYWARD", suffix="MCL", path="GEX/MELANOMA/")
## CARSON     MCL 39  12131   -branches          *WARNING* FEW SYMBOLS!
CARSON.MCL <- preproc(geoid="GSE40047", prefix="CARSON", suffix="MCL", path="GEX/MELANOMA/", keep.sample=c(1:5,21:54))
## ROSE       MCL 19  13211   2branches (unbal)     *WARNING* FEW SYMBOLS!
ROSE.MCL <- preproc(geoid="GSE22306", prefix="ROSE", suffix="MCL", path="GEX/MELANOMA/", keep.sample=c(4:22))
## HUANG      MCL 19  **requires a strsplit on "gene_assignment" to get the symbol**
HUANG.MCL <- preproc(geoid="GSE44660", prefix="HUANG", suffix="MCL", path="GEX/MELANOMA/", keep.sample=c(4:22))
## HOEK       MCL 45  21049   2 branches
HOEK.MCL <- preproc(geoid="GSE4843", prefix="HOEK", suffix="MCL", path="GEX/MELANOMA/", raw.include=F)
## CCLE       MCL 61  18901   2(3)branches (moches) **AUTOLOAD FAILED AS THERE IS NO SYMBOL FOR THIS F°CKING CUSTOM DESIGN. HOPEFULLY HAD ANOTHER NORMALIZED VERSION**
## NCI60      MCL 9   21049
NCI60.MCL <- preproc(geoid="GSE32474", prefix="NCI60", suffix="MCL", path="GEX/SPECIAL/", keep.sample=c(12,51:58), raw.include=F)
## SHAKHOVA   Melanoma cell line M010817 treated with a shRNA targetted against SOX10 profiled @ 48,96H (plus a control with scramble shRNA for each timepoint)  4  21049
SHAKHOVA.MCL <- preproc(geoid="GSE37059", prefix="SHAKHOVA", suffix="MCL", path="GEX/MELANOMA/", keep.sample=c(1,2,5,7), raw.include=F)
## PRINT      Melanoma cell line A375 treated with 45 different siRNA @ 48H (plus 3 control siRNA, and 3 non-treated replicates)  51  21049
PRINT.MCL <- preproc(geoid="GSE31534", prefix="PRINT", suffix="MCL", path="GEX/MELANOMA/", raw.include=F)

## RASKIN     PM  46  18585   
RASKIN.PM <- preproc(geoid="GSE15605", prefix="RASKIN", suffix="PM", path="GEX/MELANOMA/", keep.sample=c(17:62))
## SCATOLINI  PM  23  18841 (Agilent)
SCATOLINI.PM <- preproc(geoid="GSE12391", prefix="SCATOLINI", suffix="PM", path="GEX/MELANOMA/", keep.sample=c(19:41))
## JAEGER     ** MANUALLY LOADED **
## ECSEDI

## RASKIN     NORM  16  18585   
RASKIN.MEL.N <- preproc(geoid="GSE15605", prefix="RASKIN", suffix="MEL.NORM", path="GEX/MELANOMA/", keep.sample=c(1:16))
## SCATOLINI  NORM  18  18841 (Agilent)
SCATOLINI.NEV <- preproc(geoid="GSE12391", prefix="SCATOLINI", suffix="NEV", path="GEX/MELANOMA/", keep.sample=c(1:18))
## CARSON     NHEM  15  12131          *WARNING* FEW SYMBOLS!
CARSON.NHEM <- preproc(geoid="GSE40047", prefix="CARSON", suffix="NHEM", path="GEX/MELANOMA/", keep.sample=c(6:20))



## BREAST
## RICHARDSON   Breast tumors           40    GPL570  21049
RICHARDSON.BRE <- preproc(geoid="GSE3744", prefix="RICHARDSON", suffix="BRE", path="GEX/BREAST/", keep.sample=c(1:40))
## RICHARDSON   Breast normal           7     GPL570  21049
RICHARDSON.BRE.N <- preproc(geoid="GSE3744", prefix="RICHARDSON", suffix="BRE.NORM", path="GEX/BREAST/", keep.sample=c(41:47))
## DALFONSO     Invasive breast tumors  246   GPL570  21049       ** NOTE : FFPE SAMPLES ! **
DALFONSO.BRE <- preproc(geoid="GSE47109", prefix="DALFONSO", suffix="BRE", path="GEX/BREAST/")
## METZGER      Primary breast tumors   112   GPL570  21049
METZGER.BRE <- preproc(geoid="GSE43365", prefix="METZGER", suffix="BRE", path="GEX/BREAST/")
METZGER.BRE.SEL <- preproc(geoid="GSE43365", prefix="METZGER", suffix="BRE.sel", path="GEX/BREAST/", keep.sample=sample(1:111, 30))
## CLARKE       Breast tumors           104   GPL570  21049       ** NOTE : Well annotated, may be usefull for a breast analysis **
CLARKE.BRE <- preproc(geoid="GSE42568", prefix="CLARKE", suffix="BRE", path="GEX/BREAST/", keep.sample=c(18:121))
## CLARKE       Breast normal           17    GPL570  21049
CLARKE.BRE.N <- preproc(geoid="GSE42568", prefix="CLARKE", suffix="BRE.NORM", path="GEX/BREAST/", keep.sample=c(1:17))
## MIYAKE       ER- breast tumors       115   GPL570  21049       ** NOTE : biopsies taken before treatment **
MIYAKE.BRE <- preproc(geoid="GSE32646", prefix="MIYAKE", suffix="BRE", path="GEX/BREAST/")
## NCI60        Cell-lines              5     GPL570  21049
NCI60.BRE <- preproc(geoid="GSE32474", prefix="NCI60", suffix="BRECL", path="GEX/SPECIAL/", keep.sample=c(7:9,11,174), raw.include=F)


## LUNG
## KUNER        NSCLC AC+SCC    58    GPL570  21049
KUNER.LUN <- preproc(geoid="GSE10245", prefix="KUNER", suffix="LUN", path="GEX/LUNG/")
## TARCA        NSCLC AC+SCC    150   GPL570  21049
TARCA.LUN <- preproc(geoid="GSE43580", prefix="TARCA", suffix="LUN", path="GEX/LUNG/")
TARCA.LUN.SEL <- preproc(geoid="GSE43580", prefix="TARCA", suffix="LUN.sel", path="GEX/LUNG/", keep.sample=sample(1:150,30))
## ROUSSEAUX    Several types   293   GPL570  21049
ROUSSEAUX.LUN <- preproc(geoid="GSE30219", prefix="ROUSSEAUX", suffix="LUN", path="GEX/LUNG/", raw.include=F, keep.sample=c(1:7,9:17,19:53,55:63,65:74,76:98,100:103,105:113,115:121,123:129,131:170,175:307))
## NCI60        Cell lines    9             21049
NCI60.LUNCL <- preproc(geoid="GSE32474", prefix="NCI60", suffix="LUNCL", path="GEX/SPECIAL/", keep.sample=c(34:42), raw.include=F)


## COLON
## SCHLICKER    CRC           62    GPL570  21049
SCHLICKER.COL <- preproc(geoid="GSE35896", prefix="SCHLICKER", suffix="COL", path="GEX/COLON/")
SCHLICKER.COL.SEL <- preproc(geoid="GSE35896", prefix="SCHLICKER", suffix="COL.sel", path="GEX/COLON/", keep.sample=sample(1:62,30))
## MARISA       Primary CRC   566   GPL570  21049     ** NOTE : CIT dataset **
MARISA.COL <- preproc(geoid="GSE39582", prefix="MARISA", suffix="COL", path="GEX/COLON/", raw.include=F)
## NCI60        Cell lines    7             21049
NCI60.COLCL <- preproc(geoid="GSE32474", prefix="NCI60", suffix="COLCL", path="GEX/SPECIAL/", keep.sample=c(19:25), raw.include=F)


## SARCOMA
## DELATTRE     Primary Ewing sarcoma             117   GPL570  21049
DELATTRE.EWS <- preproc(geoid="GSE34620", prefix="DELATTRE", suffix="EWS", path="GEX/SARCOMA/")
DELATTRE.EWS.SEL <- preproc(geoid="GSE34620", prefix="DELATTRE", suffix="EWS.sel", path="GEX/SARCOMA/", keep.sample=sample(1:117, 30))
## SAVOLA       Sarcoma (Ewing, PNET, Askin)      44    GPL570  21049
SAVOLA.SAR <- preproc(geoid="GSE17618", prefix="SAVOLA", suffix="SAR", path="GEX/SARCOMA/", keep.sample=c(1:44))
## CHIBON       Various complex genetic sarcoma   182   GPL570  21049   ** Removed an outlier **
CHIBON.SAR <- preproc(geoid="GSE16382", prefix="CHIBON", suffix="SAR", path="GEX/SARCOMA/", keep.sample=c(1:34,36:183))


## OVARY
## MATEESCU     Ovarian cancer (several types)    107   GPL570  21049
MATEESCU.OVA <- preproc(geoid="GSE26193", prefix="MATEESCU", suffix="OVA", path="GEX/OVARY/")
MATEESCU.OVA.SEL <- preproc(geoid="GSE26193", prefix="MATEESCU", suffix="OVA.sel", path="GEX/OVARY/", keep.sample=sample(1:107,30))
## FERRISS      Ovarian cancer before treatment   58    GPL570  21049   ** NOTE : FFPE SAMPLES ! **
FERRISS.OVA <- preproc(geoid="GSE30161", prefix="FERRISS", suffix="OVA", path="GEX/OVARY/")
## MOK          Advanced ovarian tumors           52    GPL570  21049
MOK.OVA <- preproc(geoid="GSE18521", prefix="MOK", suffix="OVA", path="GEX/OVARY/", keep.sample=c(13:40,42:65))
## NCI60        Cell lines                        7     GPL570  21049
NCI60.OVACL <- preproc(geoid="GSE32474", prefix="NCI60", suffix="OVACL", path="GEX/SPECIAL/", keep.sample=c(10,43:48), raw.include=F)
## QIN          Ovarian cancer cell-line with TGFb treatment (0,3,6,12H) 4   GPL570  21049
QIN.OVACL <- preproc(geoid="GSE6653", prefix="QIN", suffix="OVACL", path="GEX/OVARY/", keep.sample=c(1,3,5,7), raw.include=F)
## TAN          Ovarian cancer cell-lines         34    GPL570  21049
TAN.OVACL <- preproc(geoid="GSE28724", prefix="TAN", suffix="OVACL", path="GEX/OVARY/", raw.include=F)


## BLADDER
## URQUIDI      Urothelial cancer cells           52    GPL570  21049
URQUIDI.BLA <- preproc(geoid="GSE31189", prefix="URQUIDI", suffix="BLA", path="GEX/BLADDER/", keep.sample=c(1:52))
## URQUIDI      Urothelial normal cells           40    GPL570  21049
URQUIDI.BLA.N <- preproc(geoid="GSE31189", prefix="URQUIDI", suffix="BLA.NORM", path="GEX/BLADDER/", keep.sample=c(53:92))
## RIESTER      Urothelial carcinomas             93    GPL570  21049   ** Intensities seem truncated on the lower side... ** 
RIESTER.BLA <- preproc(geoid="GSE31684", prefix="RIESTER", suffix="BLA", path="GEX/BLADDER/")


## NEUROBLASTOMA
## VONSTEDINGK  Neuroblastoma                     88    GPL570  21049
VONSTEDINGK.NBA <- preproc(geoid="GSE16476", prefix="VONSTEDINGK", suffix="NBL", path="GEX/NEUROBLASTOMA/")
VONSTEDINGK.NBA.SEL <- preproc(geoid="GSE16476", prefix="VONSTEDINGK", suffix="NBL.sel", path="GEX/NEUROBLASTOMA/", keep.sample=sample(1:88, 30))
## OHTAKI       Neuroblastoma                     51    GPL570  21049
OHTAKI.NBA <- preproc(geoid="GSE16237", prefix="OHTAKI", suffix="NBL", path="GEX/NEUROBLASTOMA/")
## VOICKMANN  Neurobastoma cell-lines             24    GPL570  21049
VOICKMANN.NBACL <- preproc(geoid="GSE28019", prefix="VOICKMANN", suffix="NBLCL", path="GEX/NEUROBLASTOMA/")


## LIVER
## CHIANG       Hepatocellular carcinoma          91    GPL570  21049
CHIANG.LIV <- preproc(geoid="GSE9843", prefix="CHIANG", suffix="LIV", path="GEX/LIVER/", raw.include=F)
CHIANG.LIV.SEL <- preproc(geoid="GSE9843", prefix="CHIANG", suffix="LIV.sel", path="GEX/LIVER/", raw.include=F, keep.sample=sample(1:91, 30))
## NOTAS        HePG2 cell-line treated with APRIL (0,2,6,12H)    4    GPL570  21049
NOTAS.LIVCL <- preproc(geoid="GSE29375", prefix="NOTAS", suffix="LIVCL", path="GEX/LIVER/", raw.include=F)


## GLIOBLASTOMA
## STURM        GBM of various ages               46    GPL570  21049     ** Removed the 3rd sample with bad boxplot and weird clustering position **
STURM.GBM <- preproc(geoid="GSE36245", prefix="STURM", suffix="GBM", path="GEX/GLIOBLASTOMA/", keep.sample=c(1,2,4:46))
## DONSON       GBM                               34    GPL570  21049
DONSON.GBM <- preproc(geoid="GSE50161", prefix="DONSON", suffix="ALL", path="GEX/GLIOBLASTOMA/", keep.sample=c(47:80))


## KIDNEY
## YANG         Papillary RCC                     34    GPL570  21049
YANG.KID <- preproc(geoid="GSE2748", prefix="YANG", suffix="KID", path="GEX/KIDNEY/", raw.include=F)
## YUSENKO      Various renal tumors              62    GPL570  21049
YUSENKO.KID <- preproc(geoid="GSE11151", prefix="YUSENKO", suffix="KID", path="GEX/KIDNEY/", keep.sample=c(1:33,36,40:67))


## THYROID
## TOMAS        PTC (Papillary carcinoma)         49    GPL570  21049
TOMAS.THY <- preproc(geoid="GSE33630", prefix="TOMAS", suffix="THY", path="GEX/THYROID/", keep.sample=c(11:59))


## EPENDYMOMA
## DONSON       Ependymoma                        46    GPL570  21049
# DONSON.EPM <- preproc(geoid="GSE50161", prefix="DONSON", suffix="EPM", path="GEX/EPENDYMOMA/", keep.sample=c(1:46))

## SPECIAL
## EXPO : A public dataset of 2158 tumors of several tissue types and origins, with constantly updated clinical annotations, on GPL570
EXPO.SPE <- preproc(geoid="GSE2109", prefix="EXPO", suffix="SPE", path="GEX/SPECIAL/", raw.include=F)
## NCI60        Cell-lines              174 (replicates)     GPL570  21049
NCI60.ALL <- preproc(geoid="GSE32474", prefix="NCI60", suffix="ALL", path="GEX/SPECIAL/", raw.include=F)



# ## MANUAL LOADING (to get samples indexes, data type, etc)




## MANUAL LOADING OF CCLE 61MCL
ccle <- read.table("MELANOMA/CCLE_GSE36133/OK_CCLE_2012_AffyCustom_CL61/CCLE_61MCL_int.txt", header=T, sep="\t", as.is=T)
cclem <- as.matrix(ccle[,-c(1)])
rownames(cclem) <- ccle$Symbol
length(which(is.na(cclem)))   ## OK, no NA
range(cclem)                  ## OK, ints and min>0
cclem <- log2(cclem)
cclem <- feat.merge(cclem, mode="median")
png("MELANOMA/CCLE_GSE36133/GSE36133_MCL_box.png", 1024, 768)
boxplot(cclem)
dev.off()
cclem <- normalizeQuantiles(cclem)
png("MELANOMA/CCLE_GSE36133/GSE36133_MCL_boxQN.png", 1024, 768)
boxplot(cclem)
dev.off()
Xhc <- hclust(dist(t(cclem), method="euclidean"), method="ward")
png("MELANOMA/CCLE_GSE36133/GSE36133_MCL_hclust_eu.wa_allgenes.png", 1600, 1000)
plot(Xhc, main=paste("CCLE_MCL", " (", nrow(cclem), ")", sep=""))
dev.off()

## MANUAL LOADING OF ALL CCLE
ccle <- goGEO(geoid="GSE36133", prefix="CCLE", suffix="", path="GEX/SPECIAL/", raw.include=F)
write.table(ccle@phenoData@data, paste("SPECIAL/CCLE_GSE36133/GSE36133_ALL_annot_samples.annot.txt", sep=""), sep="\t", quote=F, row.names=F)





## V3 banks
nselec <- 40
## MELANOMA         BOGUNOVIC     Melanoma Metastases       44    21049       2branches
BOGUNOVIC.MM.SEL    <- preproc(geoid="GSE19234", prefix="BOGUNOVIC",  suffix="MM.selv3", path="GEX/MELANOMA/", keep.sample=sample(1:44, nselec))
## BREAST           MIYAKE        ER- breast tumors         115   21049       ** NOTE : biopsies taken before treatment **
MIYAKE.BRE.SEL      <- preproc(geoid="GSE32646", prefix="MIYAKE",     suffix="BRE.selv3", path="GEX/BREAST/", keep.sample=sample(1:115, nselec))
## LUNG             ROUSSEAUX     Several subtypes          293   21049
ROUSSEAUX.LUN.SEL   <- preproc(geoid="GSE30219", prefix="ROUSSEAUX", suffix="LUN.selv3", path="GEX/LUNG/", raw.include=F, keep.sample=sample(c(1:7,9:17,19:53,55:63,65:74,76:98,100:103,105:113,115:121,123:129,131:170,175:307), nselec))
## COLON            SCHLICKER     CRC                       62    21049
SCHLICKER.COL.SEL   <- preproc(geoid="GSE35896", prefix="SCHLICKER", suffix="COL.selv3", path="GEX/COLON/", keep.sample=sample(1:62, nselec))
## EWING SARCOMA    DELATTRE      Primary Ewing sarcoma     117  21049
DELATTRE.EWS.SEL    <- preproc(geoid="GSE34620", prefix="DELATTRE", suffix="EWS.selv3", path="GEX/SARCOMA/", keep.sample=sample(1:117, nselec))
## OVARY            MOK           Advanced ovarian tumors   52    21049
MOK.OVA.SEL         <- preproc(geoid="GSE18521", prefix="MOK", suffix="OVA.selv3", path="GEX/OVARY/", keep.sample=sample(c(13:40,42:65), nselec))
## BLADDER          URQUIDI       Urothelial cancer cells   52    21049
URQUIDI.BLA.SEL     <- preproc(geoid="GSE31189", prefix="URQUIDI", suffix="BLA.selv3", path="GEX/BLADDER/", keep.sample=sample(c(1:52), nselec))
## NEUROBLASTOMA    VONSTEDINGK   Neuroblastoma             88    21049
VONSTEDINGK.NBA.SEL <- preproc(geoid="GSE16476", prefix="VONSTEDINGK", suffix="NBL.selv3", path="GEX/NEUROBLASTOMA/", keep.sample=sample(1:88, nselec))
## LIVER            CHIANG        Hepatocellular carcinoma  91    21049
CHIANG.LIV.SEL      <- preproc(geoid="GSE9843", prefix="CHIANG", suffix="LIV.selv3", path="GEX/LIVER/", raw.include=F, keep.sample=sample(1:91, nselec))
## KIDNEY           YUSENKO       Various renal tumors      62    21049
YUSENKO.KID.SEL     <- preproc(geoid="GSE11151", prefix="YUSENKO", suffix="KID.selv3", path="GEX/KIDNEY/", keep.sample=sample(c(1:33,36,40:67), nselec))


