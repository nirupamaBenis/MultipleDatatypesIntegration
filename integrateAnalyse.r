#### this script gets the correlation matrices, builds networks and sends the networks to Cytoscape, also merges the 10 networks to one

## load library and neccessary functions
library(mixOmics)
source("~/allInfoDir/workspace/cyrest.r")
source("/home/user/allInfoDir/workspace/newFunctions/functionsFound.r")

## get all the names of the networks before hand, from allDatasetsSplsSwapped.r
allCorrMix <- c("microb.metabAS", "metabAS.microb", "microb.metabACU", "metabACU.microb", "gene.microb", "microb.gene", "cyto.microb", "microb.cyto", "gene.cyto", "cyto.gene", "gene.metabAS", "metabAS.gene", "gene.metabACU", "metabACU.gene", "cyto.metabAS", "metabAS.cyto", "cyto.metabACU", "metabACU.cyto", "metabAS.metabACU", "metabACU.metabAS")

## perform the correlation twice with each dataset taking turns at being an independent variable or a dependent one 
for (i in seq(1, length(allCorrMix), by = 2)) {
  corRes <- list()
  ## build the correlation matrix
  for (j in 0:1) {
    tmpSplsObj <- loadRData(paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/", allCorrMix[i+j], ".RData"))
    corRes[[j+1]] <- corrMixSpls(tmpSplsObj, ncomp = 5)
    range(corRes[[j+1]])
  }
  corRes[[j+1]] <- t(corRes[[j+1]])
  corResAvg <- do.call("+", corRes)/length(corRes)
  write.table(corResAvg, paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/meanSpls/", allCorrMix[i+j], "Mean.txt"), sep = "\t")
  write(allCorrMix[i+j], file = "/home/user/allInfoDir/samMiceWorkSpace/integration/fileNamesAllMixMean.txt", append = T)
}

## get all the average correlation matrices to build networks
## use all the 10 dataset combinations in the for loop
fileNamesMix <- scan(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/fileNamesAllMixMean.txt", what = "")
## create variables to store information from the for loop
tmpNw <- list()
nwInfo <- data.frame(NetworkNames = fileNamesMix, cutoffPos = NA, cutoffNeg = NA, Nodes = NA, NodesX = NA, NodesY = NA)

for (i in 1:length(fileNamesMix)){
  ## read in the averaged correlation matrix
  tmpMix <- read.delim(paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/meanSpls/", fileNamesMix[i], "Mean.txt"))
  ## find the top 5% of the correlation coefficient distribution to use as threshold
  cutoff1 <- unname(quantile(as.matrix(tmpMix), 0.95))
  cutoff2 <- unname(quantile(as.matrix(tmpMix), 0.05))
  ## store the information in a data.frame
  nwInfo$cutoffPos[i] <- round(cutoff1, digits = 2)
  nwInfo$cutoffNeg[i] <- round(cutoff2, digits = 2)
  ## apply the thresholds to the correlation matrix, anything that does not pass is made 0
  filtTmpMix <- tmpMix
  for (j in 1:dim(tmpMix)[1]) {
    tmpNeg <- tmpMix[j, ]
    tmpNeg[tmpMix[j, ] >= 0] <- NA
    tmpNeg[tmpNeg > cutoff2] <- 0
    tmpPos <- tmpMix[j, ]
    tmpPos[tmpMix[j, ] <= 0] <- NA
    tmpPos[tmpPos < cutoff1] <- 0
    tmpRow <- unlist(ifelse(is.na(tmpNeg), tmpPos, tmpNeg))
    filtTmpMix[j, ] <- tmpRow
  }
  ## write the filtered correlation matrices for later use (this filtration can take some time)
  write.table(filtTmpMix, file = paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/meanSpls/", fileNamesMix[i], "MeanFilt.txt"))
  ## build the network from the filtered correlation matrices
  tmpNw[[i]] <- network(as.matrix(filtTmpMix), comp = 1:5, shape.node = c("rectangle", "rectangle"), color.node = c("deeppink", "blue"), threshold= (min(abs(cutoff1), abs(cutoff2)) - 0.1))
  ## the function network() from mixOmics opens a graph and there could be errors if that graphics device is left open
  dev.off()
  ## get the network ready to send to Cytoscape using cyrest
  tmpNames <- V(tmpNw[[i]]$gR)$name
  tmpDataNames <- unlist(strsplit(fileNamesMix[i], split = "\\."))
  V(tmpNw[[i]]$gR)$data <- replace(tmpNames, c(grep("X", tmpNames), grep("Y", tmpNames)), c(rep(tmpDataNames[2], length(grep("X", tmpNames))), rep(tmpDataNames[1], length(grep("Y", tmpNames)))))
  ## load information into the data.frame
  nwInfo$Nodes[i] <- length(V(tmpNw[[i]]$gR)$name)
  nwInfo$NodesX[i] <- length(grep("X", V(tmpNw[[i]]$gR)$name))
  nwInfo$NodesY[i] <- length(grep("Y", V(tmpNw[[i]]$gR)$name))
  V(tmpNw[[i]]$gR)$name <- V(tmpNw[[i]]$gR)$label
  igraphObjGC <- set.graph.attribute(tmpNw[[i]]$gR, "name", paste0(fileNamesMix[i], "5percent"))
  E(igraphObjGC)$interaction <- rep("undirected", length(E(igraphObjGC)))
  E(igraphObjGC)$name <- apply(get.edgelist(igraphObjGC), 1, function(x) paste(x, collapse = " (undirected) "))
  ## send the network to an open Cytoscape session
  cygraph <- toCytoscape(igraphObjGC)
  test <- send2cy(cygraph, "default%white", 'force-directed')
} 
save(tmpNw, file = paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/AllCompMeanNw5perc.RData"))
write.table(nwInfo, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/nwInfo5.txt", sep = "\t", row.names = F)

## Analyse all the networks, get a look at the network properties
tmpNw1 <- list()
allNodesAllNws <- list()
for(i in 1:10) {
  tmpNw1[[i]] <- tmpNw[[i]]$gR
  ## change to get different network properties, given below are some examples
  #   print(transitivity(tmpNw1[[i]]))
  #   print(igraph_density(tmpNw1[[i]]))
  #   print(range(degree(tmpNw1[[i]])))
  allNodesAllNws[[i]] <- V(tmpNw[[i]]$gR)$label
}
x <- unique(unlist(allNodesAllNws))
xx <- unlist(allNodesAllNws)
numOfNetworks <- unlist(lapply(lapply(x, grep, xx), length))
names(numOfNetworks) <- x
str(intersect(x, V(mergedNw)$name))
onlyOneNw <- names(numOfNetworks[numOfNetworks == 1])
onlyOneNw <- gsub("^X([[:digit:]])", "\\1", onlyOneNw)

################## merge the 10 networks ##################
mergedNw <- Reduce("%u%", tmpNw1)
round(transitivity(mergedNw), 3)
namesMerged <- V(mergedNw)$name
## also merged in Cytoscape because it is easier to get datatype names
cytoscapeMerged <- read.csv(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/mergedInCyto.csv", as.is = T)
cytoscapeMerged <- cytoscapeMerged[unlist(lapply(namesMerged, grep, cytoscapeMerged$label)), ]
V(mergedNw)$data <- cytoscapeMerged$data
igraphObj <- set.graph.attribute(mergedNw, "name", "Merged5percent")
cygraph <- toCytoscape(igraphObj)
test <- send2cy(cygraph, "default%white", 'force-directed')
edgeTableCytoscape <- read.csv(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/mergedEdgeTableChg.csv", as.is = T)
E(mergedNw)$name <- apply(get.edgelist(mergedNw), 1, function(x) paste(x, collapse = " (undirected) "))

## assign weights to the merged igraph object
mergedWeight <- NULL
for (i in E(mergedNw)$name) {
  tmpPos <- grep(i, edgeTableCytoscape$name, fixed = T)
  if (length(tmpPos) == 0) {
    i = gsub("^(.*) (\\(undirected\\)) (.*)$", "\\3 \\2 \\1", i)  
    tmpPos <- grep(i, edgeTableCytoscape$name, fixed = T)
  }
  mergedWeight[i] <- edgeTableCytoscape$weight[tmpPos]
}
E(mergedNw)$weight <-  mergedWeight

## remove few nodes that do not have meaning and change names for clarity, also remove nodes that are in only one network
mergedNwCur <- mergedNw
mergedNwCur <- delete_vertices(mergedNwCur, unlist(lapply(onlyOneNw, grep, V(mergedNwCur)$name))) ## only one network
mergedNwCur <- delete_vertices(mergedNwCur, unlist(lapply("Rik$", grep, V(mergedNwCur)$name))) ## genes that don't have information
V(mergedNwCur)$name <- gsub("^X([[:digit:]])", "\\1", V(mergedNwCur)$name)
V(mergedNwCur)$name <- gsub(".g$", "", V(mergedNwCur)$name) ## change bacterial node names for clarity
save(mergedNwCur, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/merged5NwsCur.RData")

## look at merged network
degreeNw <- degree(mergedNwCur)
hisMerg <- hist(degreeNw)
allData <- unique(V(mergedNwCur)$data)
namesMerged <- V(mergedNwCur)$name

## find connectivity hubs based on type of data of the neighbours
levelHubs <- NULL
for (i in 1:length(namesMerged)) {
  neighboursData <- unique(neighbors(graph = mergedNwCur, v = namesMerged[i])$data)
  ## make sure the neighbours are from 4 other types of data
  if (length(neighboursData) == 4 && !(neighboursData %in% V(mergedNwCur)$data[i])) {
    levelHubs <- c(levelHubs, namesMerged[i])
  }
}

## prepare the network for cytoscape
V(mergedNwCur)$levelHubs <- as.numeric(factor(V(mergedNwCur)$name %in% levelHubsDF$Names)) ## add as an attribute to the network, a column in the cytoscape node table
V(mergedNwCur)$Degree <-  degree(mergedNwCur)
V(mergedNwCur)$label <- V(mergedNwCur)$name
mergedNwCur <- set.graph.attribute(mergedNwCur, "name", "MergedCuratedLit")
E(mergedNwCur)$interaction <- rep("undirected", length(E(mergedNwCur)))
E(mergedNwCur)$name <- apply(get.edgelist(mergedNwCur), 1, function(x) paste(x, collapse = " (undirected) "))
connectionsDataPids <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/connectionsDataPids.txt", as.is = T)
E(mergedNwCur)$literatureResult <-  ifelse(apply(get.edgelist(mergedNwCur), 1, function(x) paste(x, collapse = " AND ")) %in% connectionsDataPids$connections, 1, 0)
save(mergedNwCur, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/merged5NwsCur.RData")
## send to Cytoscape
cygraph <- toCytoscape(mergedNwCur)
test <- send2cy(cygraph, "default%white", 'force-directed')

## look into the connectivity hubs in more detail
levelHubsDF <- cbind.data.frame(Names = levelHubs, DataType = V(mergedNwCur)$data[V(mergedNwCur)$name %in% levelHubs], Degree = degreeNw[names(degreeNw) %in% levelHubs])
write.table(levelHubsDF, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/mergedLevelHubsDFCur.txt", sep = "\t", row.names = F)
levelHubsDF <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/mergedLevelHubsDFCur.txt", as.is = T)
levelHubsDF$DataType <- gsubSeveral(pattern = unique(levelHubsDF$DataType), c("Microbiota", "MetabolomicsSerum", "MetabolomicsUrine", "Transcriptomics", "Cytokines"), levelHubsDF$DataType)
dataTypes <- c("MetabolomicsUrine", "MetabolomicsSerum", "Cytokines", "Transcriptomics", "Microbiota")
newDF <- data.frame(Names = character(), DataType = character(), Degree = numeric())
for (i in 1:length(dataTypes)) {
  tmpRows <- levelHubsDF[levelHubsDF$DataType == dataTypes[i], ]
  tmpRows <- tmpRows[order(tmpRows$Degree, decreasing = T), ]
  newDF <- rbind.data.frame(newDF, tmpRows)
}
write.table(newDF, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/mergedLevelHubsDFCurProp.txt", sep = "\t", row.names = F, quote = F)
