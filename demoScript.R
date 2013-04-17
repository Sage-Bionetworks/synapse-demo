require(synapseClient)
require(rGithubClient)
require(survival)
require(ggplot2)

exploratoryProjectId <- ""

## LOAD METABRIC AND OSLO VAL DATA
ovSurvEnt <- loadEntity("syn1710257")
ovSurv <- ovSurvEnt$objects$oslovalSurvData
mbSurvEnt <- loadEntity("syn1710277")
mbSurv <- mbSurvEnt$objects$metabricSurvData
mbClinEnt <- loadEntity("syn1710260")
mbClin <- mbClinEnt$objects$metabricClinicalTable
mbClin$tripleNeg <- mbClin$Her2.Expr=="-" & mbClin$ER.Expr=="-" & mbClin$PR.Expr=="-"
ovClinEnt <- loadEntity("syn1710251")
ovClin <- ovClinEnt$objects$oslovalClinicalTable
ovClin$tripleNeg <- ovClin$Her2.Expr=="-" & ovClin$ER.Expr=="-" & ovClin$PR.Expr=="-"
allDatEnt <- loadEntity("syn1752231")
mbImputedAll <- allDatEnt$objects$mbImputedAll
ovImputedAll <- allDatEnt$objects$ovImputedAll


## BUILD UP ANALYSIS DATASETS -- CENSOR ALL SURVIVAL TIMES AT 15 YEARS
mb <- mbImputedAll
mb$osYears <- mbSurv[, 1] / 365
mb$osStatus <- mbSurv[, 2]
mb$osYears[ mb$osYears >= 15 ] <- 15
mb$osStatus[ mb$osYears >= 15 ] <- 0
mb$chr8q24.3 <- mb$puf60
mb$tripleNeg <- mbClin$tripleNeg
ov <- ovImputedAll
ov$osYears <- ovSurv[, 1] / 365
ov$osStatus <- ovSurv[, 2]
ov$osYears[ ov$osYears >= 15 ] <- 15
ov$osStatus[ ov$osYears >= 15 ] <- 0
ov$chr8q24.3 <- ov$puf60
ov$tripleNeg <- ovClin$tripleNeg


## CREATE KAPLAN MEIER PLOTS BASED ON TRIPLE NEGATIVE STATUS
analysisRepo <- getRepo("Sage-Bionetworks/synapse-demo")
sourceRepoFile(analysisRepo, "utilityFunctions/mgKmPlot.R")
plotCode <- getPermlink(analysisRepo, "utilityFunctions/mgKmPlot.R")


## EXPLORE SURVIVAL WRT chr8q24.3
mbTn <- mgKmPlot(mb, "chr8q24.3", mb$tripleNeg)
mbTn
ovTn <- mgKmPlot(ov, "chr8q24.3", ov$tripleNeg)
ovTn


#####
## CREATE ACTIVITY TO SHOW HOW PLOTS WERE PRODUCED
#####
tnActivity <- Activity(name="triple negative exploration", 
                       used=list(
                         list(url=plotCode, name=basename(plotCode), wasExecuted=T),
                         list(entity=allDatEnt, wasExecuted=F),
                         list(entity=mbClinEnt, wasExecuted=F),
                         list(entity=mbSurvEnt, wasExecuted=F),
                         list(entity=ovClinEnt, wasExecuted=F),
                         list(entity=ovSurvEnt, wasExecuted=F)
                       ))
tnActivity <- storeEntity(tnActivity)


#####
## WRITE OUT AND PUSH PLOTS UP TO SYNAPSE
#####
## METABRIC
mbTnPlotFile <- file.path(tempdir(), "mb-tripleNegative-chr8q24.3-km.png")
png(mbTnPlotFile, width=900, height=600)
mbTn
dev.off()

mbTnSyn <- File(mbTnPlotFile, parentId=exploratoryProjectId)
generatedBy(mbTnSyn) <- tnActivity
mbTnSyn <- storeEntity(mbTnSyn)
tnActivity <- generatedBy(mbTnSyn)

## OSLO VALIDATION
ovTnPlotFile <- file.path(tempdir(), "ov-tripleNegative-chr8q24.3-km.png")
png(ovTnPlotFile, width=900, height=600)
ovTn
dev.off()

ovTnSyn <- File(ovTnPlotFile, parentId=exploratoryProjectId)
generatedBy(ovTnSyn) <- tnActivity
ovTnSyn <- storeEntity(ovTnSyn)

onWeb(mbTnSyn)

