require(synapseClient)
require(rGithubClient)
require(devtools)

## SOURCE IN KNITR INTEGRATION
source_gist("6117476")

folderId <- ""

helloWorld <- getPermlink('Sage-Bionetworks/synapse-demo', 'demos/helloWorld.Rmd', type="raw")
f <- synStore(File(path=helloWorld, parentId=folderId, synapseStore = F))
f <- synGet(f@properties$id)

w <- knit2synapse(file=getFileLocation(f), owner=f@properties$id)

onWeb(f)
