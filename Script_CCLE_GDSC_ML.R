library(PharmacoGx)
library(mCI)
library(PharmacoGxML)
library(Biobase)
library(RadioGx)
library(CoreGx)
library(Xeva)
library(RJSONIO)
library(rjson)
library(affycoretools)
library(affy)
library(ggplot2)
library(reshape2)
library(pander)
library(SummarizedExperiment)
library(devtools)
devtools::install("wCI-master")
devtools::install_github("bhklab/mci", ref="light")
devtools::install_github("bhklab/PharmacoGx-ML")
getwd()
availablePSets()
availablePSets(canonical=FALSE)[, c(1, 3)]
GDSC <- downloadPSet("GDSC_2020(v1-8.2)", saveDir=file.path(".", "PSets"))
CCLE <- downloadPSet("CCLE_2015", saveDir=file.path(".", "PSets"))
CMAP <- downloadPSet("CMAP_2016", saveDir=file.path(".", "PSets"))
#GDSC <- get(load('/Users/akaushik/Epigenetics/DrugRespPrediction/analysis2/PSets/GDSC_2020(v1-8.2).rds'))
#CCLE <- get(load('/Users/akaushik/Epigenetics/DrugRespPrediction/analysis2/PSets/CCLE_2015.rds'))
GDSC <- get(load('/Users/akaushik/Epigenetics/DATA/GDSC.RData'))
CCLE <- get(load('/Users/akaushik/Epigenetics/DATA/CCLE.RData'))
knitr::kable(cellInfo(CCLE)[1:3,1:3])
knitr::kable(cellInfo(CCLE))
head(cellNames(CCLE))
mycol <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462",
           "#b3de69","#fccde5","#d9d9d9","#bc80bd","#ccebc5","#ffed6f",
           "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c",
           "#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99","#b15928")
pie(table(cellInfo(CCLE)[,"tissueid"]),
    col=mycol,
    main="Tissue Types",
    radius=1,
    cex=0.8)

CCLE.auc <- PharmacoGx::summarizeSensitivityProfiles(
  CCLE,
  sensitivity.measure="auc_published",
  summary.stat="median",
  verbose=FALSE)

Tanespimycin.aac <- CCLE.auc["Tanespimycin",]
cells <- names(Tanespimycin.aac)[
  c(which.min(Tanespimycin.aac),
    which((Tanespimycin.aac > 0.2) & (Tanespimycin.aac < 0.4))[1],
    which.max(Tanespimycin.aac))]
par(mfrow=c(2, 2))
drugDoseResponseCurve(drug="Tanespimycin", cellline=cells[1],
                      pSets=CCLE, plot.type="Fitted",
                      legends.label="aac_published")
drugDoseResponseCurve(drug="Tanespimycin", cellline=cells[2],
                      pSets=CCLE, plot.type="Fitted",
                      legends.label="aac_published")
drugDoseResponseCurve(drug="Tanespimycin", cellline=cells[3],
                      pSets=CCLE, plot.type="Fitted",
                      legends.label="aac_published")


###### Experiment By Aman
common <- intersectPSet(pSets = list("CCLE"=CCLE, "GDSC"=GDSC),
                        intersectOn = c("cell.lines", "drugs"),
                        strictIntersect = TRUE)
#> Intersecting large PSets may take a long time ...
drugNames(CCLE)
drugs <- drugNames(common$CCLE)
drugs
##Example of concordant and discordant drug curves
cases <- rbind(
  #c("BT-549", "Irinotecan"),
  #c("MDA-MB-435", "Irinotecan"),

  #c("BT-549", "Lapatinib"),

  ###c("BT-549", "Nilotinib"),
  ###c("HCC1806", "Nilotinib"),
  #c("MDA-MB-435", "Nilotinib"),

  #c("MDA-MB-435", "Nutlin-3"),

  #c("BT-549", "Paclitaxel"),
  ###c("HCC1187", "Paclitaxel"),
  #c("MDA-MB-435", "Paclitaxel"),
  #c("CAL-85-1", "Paclitaxel"),
  #c("BT-20", "Paclitaxel"),
  #c("HCC1395", "Paclitaxel"),
  #c("HCC1806", "Paclitaxel"),
  #c("MDA-MB-157", "Paclitaxel"),
  #c("MDA-MB-436", "Paclitaxel"),
  #c("MDA-MB-468", "Paclitaxel"),

  ###c("CAL-85-1", "PD-0325901"),
  ###c("HCC1395", "PD-0325901"),
  #c("MDA-MB-435", "PD-0325901"),

  #c("MDA-MB-435", "Crizotinib"),

  ###c("BT-549", "PLX4720"),
  #c("MDA-MB-435", "PLX4720"),

  #c("BT-549", "Sorafenib"),
  #c("MDA-MB-435", "Sorafenib"),

  #c("BT-549", "Topotecan"),
  #c("MDA-MB-435", "Topotecan"),

  #c("MDA-MB-435", "Vandetanib"),

  ###c("BT-549", "Palbociclib"),

  #c("BT-549", "Saracatinib"),

  ###c("CAL-85-1", "Selumetinib"),
  #c("MDA-MB-435", "Selumetinib"),

  ###c("BT-549", "Tanespimycin"),
  ###c("CAL-85-1", "Tanespimycin"),
  ###c("BT-20", "Tanespimycin"),
  ###c("HCC1395", "Tanespimycin"),
  ###c("HCC1806", "Tanespimycin"),
  #c("MDA-MB-435", "Tanespimycin"),
  ###c("MDA-MB-436", "Tanespimycin"),
  ###c("MDA-MB-468", "Tanespimycin")
  )


par(mfrow=c(2, 2))
for (i in seq_len(nrow(cases))) {
  drugDoseResponseCurve(pSets=common,
                        drug=cases[i,2],
                        cellline=cases[i,1],
                        legends.label="ic50_published",
                        plot.type="Fitted",
                        ylim=c(0,130))
}

cells <- c("CAL-85-1","BT-549", "BT-20", "HCC1187", "HCC1395", "HCC1806", "MDA-MB-157", "MDA-MB-435", "MDA-MB-436", "MDA-MB-468")
par(mfrow=c(1, 10), pty = "s")


drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[1],
                      pSets=list(CCLE,GDSC), plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))
drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[2],
                      pSets=CCLE, plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))
drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[3],
                      pSets=list(CCLE,GDSC), plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))
drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[4],
                      pSets=list(CCLE,GDSC), plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))
drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[5],
                      pSets=list(CCLE,GDSC), plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))
drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[6],
                      pSets=list(CCLE,GDSC), plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))
drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[7],
                      pSets=list(CCLE,GDSC), plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))
drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[8],
                      pSets=CCLE, plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))
drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[9],
                      pSets=list(CCLE,GDSC), plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))
drugDoseResponseCurve(drug="Trastuzumab", cellline=cells[10],
                      pSets=list(CCLE,GDSC), plot.type="Actual",
                      legends.label="aac_published",
                      ylim=c(0,130))


######################## END


library(ggplot2, verbose=FALSE)
library(reshape2, verbose=FALSE)
CCLE.aac <- PharmacoGx::summarizeSensitivityProfiles(CCLE, sensitivity.measure = "aac_recomputed")
melted_data <- melt(CCLE.aac)
NA_rows <- unique(which(is.na(melted_data), arr.ind=T)[,1])
melted_data <- melted_data[-NA_rows,]
ggplot(melted_data, aes(x=Var1,y=value)) +
  geom_boxplot(fill="red") +
  theme(axis.text.x=element_text(angle=90,hjust=1)) +
  xlab("Drugs") +
  ylab("AAC")
hist(CCLE.auc["Crizotinib",], xlab="Cells response to Crizotinib)",
     col="gray", main="")

common <- intersectPSet(pSets = list("CCLE"=CCLE, "GDSC"=GDSC),
                        intersectOn = c("cell.lines", "drugs"),
                        strictIntersect = TRUE)
#> Intersecting large PSets may take a long time ...
drugs <- drugNames(common$CCLE)

##Example of concordant and discordant drug curves
cases <- rbind(
  c("BT-549", "Nilotinib"),
  c("MDA−MB−435", "TAE684"),
  c("HCC1395", "Selumetinib"),
  c("BT-549", "Selumetinib"),
  c("MDA-MB-435", "Pha-665752"),
  c("MDA-MB-436", "Pha-665752"),
  c("BT-549", "Nilotinib"),
  c("MDA-MB-435", "LBW242"),
  c("BT-549", "PHA-665752"),
  c("CAL-85-1", "Pha-665752"))

par(mfrow=c(2, 2))
for (i in seq_len(nrow(cases))) {
  drugDoseResponseCurve(pSets=common,
                        drug=cases[i,2],
                        cellline=cases[i,1],
                        legends.label="ic50_published",
                        plot.type="Fitted",
                        ylim=c(0,130))
}


##AAC scatter plot
GDSC.aac <- PharmacoGx::summarizeSensitivityProfiles(
  common$GDSC,
  sensitivity.measure='aac_recomputed',
  summary.stat="median",
  verbose=FALSE)
CCLE.aac <- summarizeSensitivityProfiles(
  common$CCLE,
  sensitivity.measure='aac_recomputed',
  summary.stat="median",
  verbose=FALSE)

GDSC.ic50 <- summarizeSensitivityProfiles(
  common$GDSC,
  sensitivity.measure='ic50_recomputed',
  summary.stat="median",
  verbose=FALSE)
CCLE.ic50 <- summarizeSensitivityProfiles(
  common$CCLE,
  sensitivity.measure='ic50_recomputed',
  summary.stat="median",
  verbose=FALSE)

drug <- "Panobinostat"
#par(mfrow=c(1, 2))

myScatterPlot(x=GDSC.aac[drug,],
              y=CCLE.aac[drug,],
              method=c("transparent"),
              transparency=0.8, pch=16, minp=50,
              xlim=c(0, max(max(GDSC.aac[drug,], na.rm=T), max(CCLE.aac[drug,], na.rm=T))),
              ylim=c(0, max(max(GDSC.aac[drug,], na.rm=T), max(CCLE.aac[drug,], na.rm=T))),
              main="Cells Response to Vandetanib",
              cex.sub=0.7,
              xlab="AAC in GDSC",
              ylab="AAC in CCLE")
legend("topright",
       legend=sprintf("r=%s\nrs=%s\nCI=%s",
                      round(cor(GDSC.aac[drug,],
                                CCLE.aac[drug,],
                                method="pearson",
                                use="pairwise.complete.obs"),
                            digits=2),
                      round(cor(GDSC.aac[drug,],
                                CCLE.aac[drug,],
                                method="spearman",
                                use="pairwise.complete.obs"),
                            digits=2),
                      round(paired.concordance.index(GDSC.aac[drug,],
                                                     CCLE.aac[drug,],
                                                     delta.pred=0,
                                                     delta.obs=0)$cindex,
                            digits=2)),
       bty="n")


c_index <-  mc_index <- NULL
for(drug in drugs){
  tt <- mCI::paired.concordance.index(GDSC.aac[drug,], CCLE.aac[drug,], delta.pred=0, delta.obs=0, alternative="greater")
  c_index <- c(c_index, tt$cindex)
  tt <- mCI::paired.concordance.index(GDSC.aac[drug,], CCLE.aac[drug,], delta.pred=0.2, delta.obs=0.2, alternative="greater", logic.operator="or")
  mc_index <- c(mc_index, tt$cindex)
}
mp <- barplot(as.vector(rbind(c_index, mc_index)), beside=TRUE, col=c("black", "orange"), ylim=c(0, 1), ylab="Concordance Index", space=c(.15,.85), border=NA, main="rCI")
text(mp, par("usr")[3], labels=as.vector(rbind(drugs, rep("", 15))), srt=45, adj=c(1.1,1.1), xpd=TRUE, cex=.8)
abline(h=.7, lty=2)

############################## END

############################# Heat Map
library(Xeva)
downloadXevaSet(name = NULL, saveDir = file.path(".", "XevaSet"),
                XevaSetFileName = NULL, verbose = TRUE)
Xeva::downloadXevaSet()
##to download a dataset
brca_full = downloadXevaSet(name="PDXE_BRCA", saveDir="XevaSet")
brca <- brca_full
response(brca, model.id="X.1004.BG98", res.measure="mRECIST")
brca.mr <- summarizeResponse(brca, response.measure = "mRECIST",
                             group.by="patient.id")
plotmRECIST(brca.mr, control.name = "untreated")
############################ END

################## Biomarkers
## Extracting the first 100 gene names. featureInfo brings up the rowData of the "rna" SummarizedExperiment
chosen.genes <- sample(rownames(featureInfo(GDSC, "rna")))
#"Crizotinib", "Lapatinib", "Nilotinib", "Paclitaxel", "Palbociclib", "Panobinostat", "PD-0325901", "PLX4720", "Saracatinib", "Selumetinib", "Sorafenib",  "Tanespimycin"
sigs <- PharmacoGx::drugSensitivitySig(GDSC, "rna", sensitivity.measure = "aac_recomputed", drugs=c("Lapatinib", "Nilotinib", "Panobinostat", "PLX4720", "Selumetinib", "Tanespimycin"), features = chosen.genes)

knitr::kable(x = head(sigs[,1,]))
knitr::kable(x = head(as.data.frame(sigs[,2,])))
getwd()
View(sigs)
write.table(sigs, file = "/Users/akaushik/Downloads/PharmacoGxML/TNBC100BioMarkersNamesGDSC", sep="\t")


##### Second Part (Radiation Sensitivity Signatures)

protein_coding <- which(ft.info$GeneBioType %in% "protein_coding")
radSensSig <- radSensitivitySig(Cleveland, mDataType='rna', nthread=1, features=fNames(Cleveland, 'rna')[protein_coding])

# Convert to a data.frame
radSensSigDF <- data.frame(radSensSig@.Data[, , ], row.names=rownames(radSensSig))

# Order by estimated correlation
radSensSigDF <- radSensSigDF[order(radSensSigDF$estimate, decreasing=TRUE), ]

# Subset to only significant FDR
radSensSigDF <- radSensSigDF[radSensSigDF$fdr < 0.05, ]
radSensGeneSig <- head(radSensSigDF, )
radResistGeneSig <- tail(radSensSigDF[radSensSigDF$fdr < 0.05, ])
write.table(radSensSigDF, file = "/Users/akaushik/Downloads/PharmacoGxML/TNBCBioMarkersRadio-sensitivity", sep="\t")
write.table(radResistGeneSig, file = "/Users/akaushik/Downloads/PharmacoGxML/TNBCBioMarkersRadio-resistance", sep="\t")
knitr::kable(radSensGeneSig, caption="Top 10 Gene Markers of Radio-sensitivity in the Cleveland RSet")
knitr::kable(radResistGeneSig, caption="Top 10 Gene Markers of Radio-resistance in the Cleveland RSet")
###################### Comparing Sensitivity Signatures between Radiation and Drug Response

library(PharmacoGx)
library(RadioGx)

#GDSC1 <- downloadPSet("GDSC_2019(v1_8.0)")
Cleveland <- downloadRSet("Cleveland")
ft.rad <- rownames(featureInfo(Cleveland, "rna"))[which(featureInfo(Cleveland, "rna")$GeneBioType == "protein_coding")]

ft.info <- featureInfo(GDSC, "rna")
ft.drug <- rownames(ft.info)[which(ft.info$GeneBioType == "protein_coding")]
#"Lapatinib", "Nilotinib", "Panobinostat", "PLX4720", "Selumetinib", "Tanespimycin"
drugs.selected <- c("Crizotinib", "Lapatinib", "Nilotinib", "Paclitaxel", "Palbociclib", "Panobinostat", "PD-0325901", "PLX4720", "Saracatinib", "Selumetinib", "Sorafenib",  "Tanespimycin")
nthread <- 16
radSig <- radSensitivitySig(Cleveland, mDataType = "rna", features=ft.rad, nthread=nthread)
gdscSigs <- PharmacoGx::drugSensitivitySig(GDSC, mDataType="rna", sensitivity.measure="aac_recomputed", drugs=drugs.selected, features=ft.drug, nthread=nthread)

common.genes <- intersect(rownames(gdscSigs), rownames(radSig))
radSig <- radSig[common.genes,,,drop=FALSE]
gdscSigs <- gdscSigs[common.genes,,]
connectivity.res <- lapply(drugs.selected, function(drug){
  return(connectivityScore(x=radSig[,1,c("estimate", "pvalue")], y=gdscSigs[,drug, c("estimate","pvalue")], method = "gwc", nperm = 200))
})

names(connectivity.res) <- drugs.selected

print(connectivity.res)


##################  END

library(mRMRe, verbose=FALSE)
library(Biobase, verbose=FALSE)
library(Hmisc, verbose=FALSE)
library(glmnet, verbose=FALSE)
library(caret, verbose=FALSE)
library(randomForest, verbose=FALSE)

##Preparing trainig dataset
train_expr0 <- PharmacoGx::summarizeMolecularProfiles(GDSC, mDataType="rna", fill.missing=FALSE, verbose=FALSE)
train_expr <- t(assay(train_expr0))
aac <- summarizeSensitivityProfiles(GDSC, sensitivity.measure="aac_recomputed", drug="Nilotinib", fill.missing=FALSE, verbose=FALSE)
cells <- intersect(rownames(train_expr), names(aac))
df <- as.matrix(cbind(train_expr[cells,], "Nilotinib"=aac[cells]))


##Preparing validation dataset
validation_expr <- PharmacoGx::summarizeMolecularProfiles(GDSC, mDataType="rna", fill.missing=FALSE, verbose=FALSE)
actual_labels <- PharmacoGx::summarizeSensitivityProfiles(GDSC, sensitivity.measure="aac_recomputed", drug="Nilotinib", fill.missing=FALSE, verbose=FALSE)
######### Aman Start
#method=c("ridge", "lasso", "random_forest", "svm", "elastic_net"),
for(method in c("ridge")){
  #par(mfrow=c(1, 5))
  res <- optimization(train=df[, -ncol(df), drop=F],
                      labels=t(df[, ncol(df), drop=F]),
                      method=method,
                      folds.no=5,
                      sampling.no=1,
                      features.no=10,
                      feature.selection="mRMR",
                      assessment=c("corr", "CI", "rCI", "mCI"))

  validation_labels <- validation(model=res$model$Nilotinib,
                                  validation.set=t(assay(validation_expr)),
                                  validation.labels=actual_labels,
                                  method=method,
                                  assessment=("rCI"))

}
######### Aman End
########### svm

  par(mfrow=c(1, 5))
  method <- "svm"
  res <- optimization(train=df[, -ncol(df), drop=F],
                      labels=t(df[, ncol(df), drop=F]),
                      method="svm",
                      folds.no=5,
                      sampling.no=1,
                      features.no=10,
                      feature.selection="mRMR",
                      assessment=c("corr", "CI", "rCI", "mCI"))
  model <- res$model$Nilotinib
  validation_expr <- PharmacoGx::summarizeMolecularProfiles(GDSC, mDataType="rna", fill.missing=FALSE)
  actual_labels <- PharmacoGx::summarizeSensitivityProfiles(GDSC, sensitivity.measure="aac_recomputed", drug="Nilotinib", fill.missing=FALSE)
  validation.set=t(assay(validation_expr))
  validation_labels <-validation(model=res$model$Nilotinib,
                                 validation.set=validation.set,
                                 validation.labels=actual_labels,
                                 method = "svm",
                                 assessment= "rCI")

  ########### elastic_net

  par(mfrow=c(1, 5))
  method <- "elastic_net"
  res <- optimization(train=df[, -ncol(df), drop=F],
                      labels=t(df[, ncol(df), drop=F]),
                      method="elastic_net",
                      folds.no=5,
                      sampling.no=1,
                      features.no=10,
                      feature.selection="mRMR",
                      assessment=c("corr", "CI", "rCI", "mCI"))
  model <- res$model$Nilotinib
  validation_expr <- PharmacoGx::summarizeMolecularProfiles(CCLE, mDataType="rna", fill.missing=FALSE)
  actual_labels <- PharmacoGx::summarizeSensitivityProfiles(CCLE, sensitivity.measure="aac_recomputed", drug="Nilotinib", fill.missing=FALSE)
  validation.set=t(assay(validation_expr))
  validation_labels <-validation(model=model,
                                 validation.set=validation.set,
                                 validation.labels=actual_labels,
                                 method = "elastic_net",
                                 assessment= "rCI")

  ########### random_forest

  par(mfrow=c(1, 5))
  method <- "random_forest"
  res <- optimization(train=df[, -ncol(df), drop=F],
                      labels=t(df[, ncol(df), drop=F]),
                      method="random_forest",
                      folds.no=5,
                      sampling.no=1,
                      features.no=10,
                      feature.selection="mRMR",
                      assessment=c("corr", "CI", "rCI", "mCI"))
  model <- res$model$Nilotinib
  validation_expr <- PharmacoGx::summarizeMolecularProfiles(GDSC, mDataType="rna", fill.missing=FALSE)
  actual_labels <- PharmacoGx::summarizeSensitivityProfiles(GDSC, sensitivity.measure="aac_recomputed", drug="Nilotinib", fill.missing=FALSE)
  validation.set=t(assay(validation_expr))
  validation_labels <-validation(model=res$model$Nilotinib,
                                 validation.set=validation.set,
                                 validation.labels=actual_labels,
                                 method = "random_forest",
                                 assessment= "rCI")

  ###########  lasso

  par(mfrow=c(1, 3))
  method <- "lasso"
  res <- optimization(train=df[, -ncol(df), drop=F],
                      labels=t(df[, ncol(df), drop=F]),
                      method="lasso",
                      folds.no=5,
                      sampling.no=1,
                      features.no=10,
                      feature.selection="mRMR",
                      assessment=c("corr", "CI", "rCI", "mCI"))
  model <- res$model$Nilotinib
  validation_expr <- PharmacoGx::summarizeMolecularProfiles(CCLE, mDataType="rna", fill.missing=FALSE)
  actual_labels <- PharmacoGx::summarizeSensitivityProfiles(CCLE, sensitivity.measure="aac_recomputed", drug="Nilotinib", fill.missing=FALSE)
  validation.set=t(assay(validation_expr))
  validation_labels <-validation(model=res$model$Nilotinib,
                                 validation.set=validation.set,
                                 validation.labels=actual_labels,
                                 method = "lasso",
                                 assessment= "rCI")

  ########### ridge

  par(mfrow=c(1, 3))
  method <- "ridge"
  res <- optimization(train=df[, -ncol(df), drop=F],
                      labels=t(df[, ncol(df), drop=F]),
                      method="ridge",
                      folds.no=5,
                      sampling.no=1,
                      features.no=10,
                      feature.selection="mRMR",
                      assessment=c("corr", "CI", "rCI", "mCI"))
  model <- res$model$Nilotinib
  validation_expr <- PharmacoGx::summarizeMolecularProfiles(CCLE, mDataType="rna", fill.missing=FALSE)
  actual_labels <- PharmacoGx::summarizeSensitivityProfiles(CCLE, sensitivity.measure="aac_recomputed", drug="Nilotinib", fill.missing=FALSE)
  validation.set=t(assay(validation_expr))
  validation_labels <-validation(model=model,
                                 validation.set=validation.set,
                                 validation.labels=actual_labels,
                                 method = "ridge",
                                 assessment=("rCI"))


######################### Known Biomarkers
  #"Lapatinib", "Nilotinib", "Panobinostat", "PLX4720", "Selumetinib", "Tanespimycin"
  features <- PharmacoGx::fNames(CCLE, "rna")[
    which(featureInfo(CCLE,
                      "rna")$Symbol == "YIPF5")]
  ccle.sig.rna <- PharmacoGx::drugSensitivitySig(CCLE,
                                     mDataType="rna",
                                     drugs=c("Lapatinib"),
                                     features=features,
                                     sensitivity.measure="aac_published",
                                     molecular.summary.stat="median",
                                     sensitivity.summary.stat="median",
                                     verbose=FALSE)
  gdsc.sig.rna <- PharmacoGx::drugSensitivitySig(GDSC,
                                     mDataType="rna",
                                     drugs=c("Lapatinib"),
                                     features=features,
                                     sensitivity.measure="aac_recomputed",
                                     molecular.summary.stat="median",
                                     sensitivity.summary.stat="median",
                                     verbose=FALSE)
  ccle.sig.mut <- PharmacoGx::drugSensitivitySig(CCLE,
                                     mDataType="mutation",
                                     drugs=c("Lapatinib"),
                                     features="YIPF5",
                                     sensitivity.measure="aac_published",
                                     molecular.summary.stat="and",
                                     sensitivity.summary.stat="median",
                                     verbose=FALSE)
  gdsc.sig.mut <- PharmacoGx::drugSensitivitySig(GDSC,
                                     mDataType="mutation",
                                     drugs=c("Lapatinib"),
                                     features="YIPF5",
                                     sensitivity.measure="aac_recomputed",
                                     molecular.summary.stat="and",
                                     sensitivity.summary.stat="median",
                                     verbose=FALSE)
  ccle.sig <- rbind(ccle.sig.rna, ccle.sig.mut)
  gdsc.sig <- rbind(gdsc.sig.rna, gdsc.sig.mut)
  known.biomarkers <- cbind("GDSC effect size"=gdsc.sig[,1],
                            "GDSC pvalue"=gdsc.sig[,6],
                            "CCLE effect size"=ccle.sig[,1],
                            "CCLE pvalue"=ccle.sig[,6])
  rownames(known.biomarkers) <- c("Lapatinib + YIPF5")
  library(xtable, verbose=FALSE)
  xtable(known.biomarkers, digits=c(0, 2, -1, 2, -1), caption='Concordance of biomarkers across stuudies')
  par(mfrow=c(2, 2))
  CCLE_expr0 <- PharmacoGx::summarizeMolecularProfiles(CCLE, mDataType="rna", fill.missing=FALSE, verbose=FALSE)
  CCLE_expr <- t(assay(CCLE_expr0))
  CCLE_cells <- intersect(rownames(CCLE_expr), colnames(CCLE.aac))
  plot(CCLE.aac["Lapatinib", CCLE_cells], CCLE_expr[CCLE_cells, features],
       main="CCLE + Lapatinib + YIPF5",
       cex.main=1, ylab="Predictions", xlab="Drug Sensitivity", pch=20, col="red")

  GDSC_expr0 <- PharmacoGx::summarizeMolecularProfiles(GDSC, mDataType="rna", fill.missing=FALSE, verbose=FALSE)
  GDSC_expr <- t(assay(GDSC_expr0))
  #> Summarizing rna molecular data for:  GDSC
  GDSC_cells <- intersect(rownames(GDSC_expr), colnames(GDSC.aac))
  plot(GDSC.aac["Lapatinib", GDSC_cells], GDSC_expr[GDSC_cells, features],
       main="GDSC + Lapatinib + YIPF5",
       cex.main=1, ylab="Predictions", xlab="Drug Sensitivity", pch=20, col="blue")

  CCLE_mut0 <- PharmacoGx::summarizeMolecularProfiles(CCLE, mDataType="rna", fill.missing=FALSE, verbose=FALSE)
  CCLE_mut <- t(assay(CCLE_mut0))
  CCLE_cells <- intersect(rownames(CCLE_mut), colnames(CCLE.aac))
  boxplot(CCLE.aac["Lapatinib", CCLE_cells]~ CCLE_mut[CCLE_cells, "YIPF5"], col="red", pch=20, main="CCLE + Lapatinib + YIPF5",
          cex.main=1, xlab="Mutation", ylab="Drug Sensitivity")

  GDSC_mut0 <- PharmacoGx::summarizeMolecularProfiles(GDSC, mDataType="rna", fill.missing=FALSE, verbose=FALSE)
  GDSC_mut <- t(assay(GDSC_mut0))
  GDSC_cells <- intersect(rownames(GDSC_mut), colnames(GDSC.aac))
  boxplot(GDSC.aac["Lapatinib", GDSC_cells]~ GDSC_mut[GDSC_cells, "YIPF5"], col="blue", pch=20, main="GDSC + Lapatinib + YIPF5",
          cex.main=1, xlab="Mutation", ylab="Drug Sensitivity")



  ####################### END


  ############## Further analysis
  library(Biobase)
  library(SummarizedExperiment)
  library(S4Vectors)
  library(PharmacoGx)
  #data("GDSCsmall")
  #data("CCLEsmall")
  commonGenes <- intersect(fNames(GDSC, "rna"),
                           fNames(CCLE,"rna"))
  common <- intersectPSet(list('CCLE'=CCLE,
                               'GDSC'=GDSC),
                          intersectOn=c("cell.lines", "drugs"),
                          strictIntersect=TRUE)


  GDSC.aac <- summarizeSensitivityProfiles(
    common$GDSC,
    sensitivity.measure='aac_recomputed',
    summary.stat="median",
    verbose=FALSE)
  CCLE.aac <- summarizeSensitivityProfiles(
    common$CCLE,
    sensitivity.measure='aac_recomputed',
    summary.stat="median",
    verbose=FALSE)

  GDSC.ic50 <- summarizeSensitivityProfiles(
    common$GDSC,
    sensitivity.measure='ic50_recomputed',
    summary.stat="median",
    verbose=FALSE)
  CCLE.ic50 <- summarizeSensitivityProfiles(
    common$CCLE,
    sensitivity.measure='ic50_recomputed',
    summary.stat="median",
    verbose=FALSE)

  GDSCexpression <- PharmacoGx::summarizeMolecularProfiles(common$GDSC,
                                               cellNames(common$GDSC),
                                               mDataType="rna",
                                               features=commonGenes,
                                               verbose=FALSE)
  CCLEexpression <- PharmacoGx::summarizeMolecularProfiles(common$CCLE,
                                               cellNames(common$CCLE),
                                               mDataType="rna",
                                               features=commonGenes,
                                               verbose=FALSE)
  gg <- fNames(common[[1]], 'rna')
  cc <- cellNames(common[[1]])

  ge.cor <- sapply(cc, function (x, d1, d2) {
    stats::cor(d1[ , x], d2[ , x], method="spearman",
               use="pairwise.complete.obs")
    ## TO DO:: Ensure all assays are name so we can call by name instead of index
  }, d1=assay(GDSCexpression, 1), d2=assay(CCLEexpression, 1))
  ic50.cor <- sapply(cc, function (x, d1, d2) {
    stats::cor(d1[, x], d2[ , x], method="spearman",
               use="pairwise.complete.obs")
  }, d1=GDSC.ic50, d2=CCLE.ic50)
  auc.cor <- sapply(cc, function (x, d1, d2) {
    stats::cor(d1[ , x], d2[ , x], method="spearman",
               use="pairwise.complete.obs")
  }, d1=GDSC.aac, d2=CCLE.aac)

  w1 <- stats::wilcox.test(x=ge.cor, y=auc.cor,
                           conf.int=TRUE, exact=FALSE)
  w2 <- stats::wilcox.test(x=ge.cor, y=ic50.cor,
                           conf.int=TRUE, exact=FALSE)
  yylim <- c(-1, 1)
  ss <- sprintf("GENE vs. AAC = %.1E\nGENE vs. IC50 = %.1E",
                w1$p.value, w2$p.value)

  boxplot(list("GENE"=ge.cor,
               "AAC"=auc.cor,
               "IC50"=ic50.cor),
          main="Concordance b/w Cancer Cell Lines",
          ylab=expression(Expression-R[s]),
          sub=ss,
          ylim=yylim,
          col="#F2AA4CFF",
          pch=20,
          border="#101820FF")

  library(PharmacoGx)
  library(pander)
  #data(CMAPsmall)
  CMAP <- downloadPSet("CMAP_2016")
  drug.perturbation <- drugPerturbationSig(CMAP,
                                           mDataType="rna",
                                           verbose=FALSE)
  data(HDAC_genes)

  res <- apply(drug.perturbation[,,c("tstat", "fdr")],
               2, function(x, HDAC){
                 return(PharmacoGx::connectivityScore(x=x,
                                                      y=HDAC[,2,drop=FALSE],
                                                      method="fgsea", nperm=100))
               }, HDAC=HDAC_genes)

  rownames(res) <- c("Connectivity", "P Value")
  res <- t(res)
  res <- res[order(res[,1], decreasing=TRUE),]
  pander::pandoc.table(res,
                       caption='Connectivity Score results for HDAC inhibitor gene signature.',
                       style = 'rmarkdown')

  library(pander)
  #data(CCLE)
  features <- fNames(CCLE, "rna")[
    which(featureInfo(CCLE,
                      "rna")$Symbol == "YIPF5")]
  sig.rna <- drugSensitivitySig(object=CCLE,
                                mDataType="rna",
                                drugs=c("Lapatinib"),
                                features=features,
                                sensitivity.measure="auc_published",
                                molecular.summary.stat="median",
                                sensitivity.summary.stat="median",
                                verbose=FALSE)
  sig.mut <- drugSensitivitySig(object=CCLE,
                                mDataType="mutation",
                                drugs=c("Lapatinib"),
                                features="YIPF5",
                                sensitivity.measure="auc_published",
                                molecular.summary.stat="and",
                                sensitivity.summary.stat="median",
                                verbose=FALSE)
  sig <- rbind(sig.rna, sig.mut)
  rownames(sig) <- c("Lapatinib + YIPF5","Lapatinib + YIPF5")
  colnames(sig) <- dimnames(sig.mut)[[3]]
  pander::pandoc.table(t(sig), style = "rmarkdown", caption='P Value of Gene-Drug Association' )

  #############  END


## download and process the HDAC signature
mydir <- "1132939s"
downloader::download(paste("http://www.sciencemag.org/content/suppl/2006/09/29/313.5795.1929.DC1", mydir, ".zip", sep=""), destfile=paste(mydir,".zip",sep=""))
unzip(paste(mydir,".zip", sep=""))
getwd()
library(hgu133a.db)
library(PharmacoGx)

HDAC_up <- gdata::read.xls(paste(mydir, paste(mydir, "sigS1.xls", sep="_"),sep="/"), sheet=1, header=FALSE, as.is=TRUE)
HDAC_down <- gdata::read.xls(paste(mydir, paste(mydir, "sigS1.xls", sep="_"),sep="/"), sheet=2, header=FALSE, as.is=TRUE)
HDAC <- as.data.frame(matrix(NA, nrow=nrow(HDAC_down)+nrow(HDAC_up), ncol=2))
annot <- AnnotationDbi::select(hgu133a.db, keys = c(HDAC_up[[1]], HDAC_down[[1]]), columns=c("ENSEMBL"), keytype="PROBEID")
gene_up <- unique(annot[match(HDAC_up[[1]], annot[,1]),2])
gene_down <- na.omit(unique(annot[match(HDAC_down[[1]], annot[,1]),2]))
HDAC_genes <- as.data.frame(matrix(NA, nrow=length(gene_down)+length(gene_up), ncol=2))


HDAC_genes[ , 2] <- c(rep(1, times=length(gene_up)), rep(-1, times=length(gene_down)))
HDAC_genes[ , 1] <- c(gene_up, gene_down)
rownames(HDAC_genes) <- HDAC_genes[ , 1]
HDAC <- HDAC_genes[ , 2]
names(HDAC) <- rownames(HDAC_genes)

drug.perturbation <- PharmacoGx::downloadPertSig("CMAP_2016")
dimnames(drug.perturbation)[[1]] <- gsub("_at", "", dimnames(drug.perturbation)[[1]])
library(boot)
library(snow)
message("Be aware that computing sensitivity will take some time...")
cl <- parallel::makeCluster(2)
res <- parApply(drug.perturbation[ , , c("tstat", "fdr")], 2, function(x, HDAC){
  return(PharmacoGx::connectivityScore(x=x, y=HDAC, method="fgsea", nperm=100))
}, cl=cl, HDAC=HDAC)
stopCluster(cl)
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)


res <- apply(drug.perturbation[ , , c("tstat", "fdr")], 2, function(x, HDAC){
  return(PharmacoGx::connectivityScore(x=x, y=HDAC, method="fgsea", nperm=100))
}, HDAC=HDAC)
rownames(res) <- c("Connectivity", "P Value")
res <- t(res)


HDAC_inhibitors <- c("vorinostat")

#TRY# HDAC_inhibitors <- c("vorinostat", "trichostatin A", "HC toxin", "valproic acid")
res <- res[order(res[,1], decreasing=T), ]
HDAC_ranks <- which(rownames(res) %in% HDAC_inhibitors)## download and process the HDAC signature

