require(devtools)
require(synapseClient)

## SOURCE IN KNITR INTEGRATION
source_gist("6117476")

folderId <- ""

helloWorld <- "~/workspace/repos/synapse-demo/demos/helloWorld.Rmd"
f <- synStore(File(path=helloWorld, parentId=folderId))

w <- knit2synapse(file=getFileLocation(f), owner=f@properties$id)

onWeb(f)
