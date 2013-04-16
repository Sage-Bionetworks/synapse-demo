require(synapseClient)
require(rGithubClient)
require(Biobase)
require(survival)
require(illuminaHumanv3.db)
require(ggplot2)

exploratoryProjectId <- "syn1761555"

## GET OSLO VALIDATION CLINICAL DATA
ovClinEnt <- loadEntity("syn1710251")
ovClin <- ovClinEnt$objects$oslovalClinicalTable
ovClin$tripleNeg <- ovClin$Her2.Expr=="-" & ovClin$ER.Expr=="-" & ovClin$PR.Expr=="-"

## GET SURVIVAL DATA
ovSurvEnt <- loadEntity("syn1710257")
ovSurv <- ovSurvEnt$objects$oslovalSurvData
ovClin$osYears <- ovSurv[, 1] / 365
ovClin$osStatus <- ovSurv[, 2]

## READ IN THE OSLO VALIDATION EXPRESSION MATRIX
osloData <- loadEntity("syn1710255")
osloExpression <- exprs(osloData$objects$oslovalExprData)

## GET EXPRESSION ANNOTATIONS FROM BIOCONDUCTOR
mapDat <- as.list(illuminaHumanv3ALIAS2PROBE)

## READ IN THE LIST OF ATTRACTOR METAGENES
metagenesData <- loadEntity("syn1725748")
## GRAB chr8q24-3 METAGENE
chr8q24.3 <- readLines(file.path(metagenesData$cacheDir, "chr8q24-3.txt"))
theseProbes <- unlist(mapDat[chr8q24.3])
mgExpr <- colMeans(osloExpression[theseProbes[!is.na(theseProbes)], ])
mgExpr <- mgExpr - median(mgExpr)
ovClin$chr8q24.3 <- mgExpr


## CREATE KAPLAN MEIER PLOTS BASED ON TRIPLE NEGATIVE STATUS
analysisRepo <- getRepo("Sage-Bionetworks/synapse-demo")
sourceRepoFile(analysisRepo, "utilityFunctions/mgKmPlot.R")
plotCode <- getPermlink(analysisRepo, "utilityFunctions/mgKmPlot.R")

nonTnPlot <- mgKmPlot(ovClin, "chr8q24.3", subset=!ovClin$tripleNeg)
nonTnPlot

tnPlot <- mgKmPlot(ovClin, "chr8q24.3", subset=ovClin$tripleNeg)
tnPlot

#####
## CREATE ACTIVITY TO SHOW HOW PLOTS WERE PRODUCED
#####
tnActivity <- Activity(name="triple negative exploration", 
                       used=list(
                         list(url=plotCode, name=basename(plotCode), wasExecuted=T),
                         list(entity=osloData, wasExecuted=F),
                         list(entity=ovClinEnt, wasExecuted=F),
                         list(entity=metagenesData, wasExecuted=F),
                         list(entity=ovSurvEnt, wasExecuted=F)
                       ))
tnActivity <- storeEntity(tnActivity)


#####
## WRITE OUT AND PUSH PLOTS UP TO SYNAPSE
#####
## NON TRIPLE NEGATIVE
nonTnPlotFile <- file.path(tempdir(), "chr8q24.3-km-nonTripleNegative.png")
png(nonTnPlotFile, width=900, height=600)
nonTnPlot
dev.off()

nonTnSyn <- File(nonTnPlotFile, parentId=exploratoryProjectId)
generatedBy(nonTnSyn) <- tnActivity
nonTnSyn <- storeEntity(nonTnSyn)
tnActivity <- generatedBy(nonTnSyn)

## TRIPLE NEGATIVE
tnPlotFile <- file.path(tempdir(), "chr8q24.3-km-tripleNegative.png")
png(tnPlotFile, width=900, height=600)
tnPlot
dev.off()

tnSyn <- File(tnPlotFile, parentId=exploratoryProjectId)
generatedBy(tnSyn) <- tnActivity
tnSyn <- storeEntity(tnSyn)


onWeb(tnSyn)

