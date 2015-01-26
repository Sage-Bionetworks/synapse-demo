require(synapseClient)
require(rGithubClient)
require(ggplot2)

folderId <- ""

myRepo <- getRepo("Sage-Bionetworks/synapse-demo")
thisCode <- getPermlink(myRepo, "demos/demo-201501.R")

## SUBTYPE CALLS
subSyn <- synGet("syn2755232")
subtypes <- read.csv(getFileLocation(subSyn), as.is=T, na.strings=c("NOLBL", "UNK"))
tcgaSub <- subtypes[ subtypes$dataset == "tcga_rnaseqAll", ]
tcgaSub <- tcgaSub[ !is.na(tcgaSub$CMS4network_plus_classifier_in_noncore_samples), ]
rownames(tcgaSub) <- tcgaSub$sample

## TCGA EXPRESSION DATA FROM COAD AND READ
crcExpr <- synGet("syn2326100")
crc <- read.delim(getFileLocation(crcExpr), as.is=T, row.names=1, check.names=FALSE)
crc <- crc[apply(crc, 1, sd) !=0, ]

## FIND INTERSECTION
these <- intersect(colnames(crc), rownames(tcgaSub))
tcgaSub <- tcgaSub[ these, ]
crc <- crc[, these]

## SPECIFICALLY INTERESTED IN EDN2
plotDF <- data.frame(expression=as.numeric(crc["EDN2",]), subtype=tcgaSub$CMS4network_plus_classifier_in_noncore_samples)

## PLOT EXPRESSION VALUES BY SUBTYPE
p <- ggplot(data=plotDF, aes(x=expression, fill=factor(subtype))) + 
  geom_density(alpha=0.3)
p
filePath <- file.path(tempdir(), "edn2-by-subtype.png")
png(filePath, height=400, width=600)
show(p)
dev.off()

f <- synStore(File(path=filePath, parentId=folderId),
              activityName="EDN2 in CRC Subtypes",
              used=list(
                list(url=thisCode, name=basename(thisCode), wasExecuted=T),
                list(entity=subSyn, wasExecuted=F),
                list(entity=crcExpr, wasExecuted=F)))

