#### This script uses the names of all the nodes in the merged network to search for co-occurence of connected nodes in literature

## get edges for literature search
library(igraph)
load(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/merged5Nws.RData")
mergedNwDF <- get.data.frame(mergedNw)
mergedNwDF <- mergedNwDF[, c(1,2)]
write.table(mergedNwDF, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/merged5NwsEdgesDF.txt", sep = "\t")

## novel connections
s24Neighbors <- neighbors(mergedNwCur, v= 3)$name
s24Genes <- s24Neighbors[s24Neighbors %in% V(mergedNwCur)$name[V(mergedNwCur)$data == "gene"]]
write(s24Genes, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/s24GeneNeighbors.txt", sep = "\t")
geneID2GO <- readMappings(file="geneID2GOAll.txt")
refHumanMouse <- read.delim(file = "/home/user/allInfoDir/workspace/conversion/refHumanMouse.txt", as.is = T)
s24GenesH <- refHumanMouse$symbHuman[refHumanMouse$symbMouse %in% s24Genes]
write(s24GenesH, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/s24GeneNeighborsHumanSymb.txt", sep = "\t")
s24GenesEnt <- refHumanMouse$entrezHuman[refHumanMouse$symbMouse %in% s24Genes]
goInputGenes <- factor(as.integer(names(geneID2GO) %in% s24GenesEnt))
names(goInputGenes) <- names(geneID2GO)
GOdataBP <- new("topGOdata", ontology = "BP", allGenes = goInputGenes, annot = annFUN.gene2GO, gene2GO = geneID2GO)
anno.genes <- genesInTerm(GOdataBP) #opp to geneID2GO, for each term the list of genes
goNames <- names(anno.genes)
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")#enrichment is done here, there are other statistical tests you can choose
resultFisher <- getSigGroups(GOdataBP, test.stat)
pvalFis <- score(resultFisher)#gives named vector of pvalues
pvalFis <- sort(pvalFis)
pval.01 <- sum(pvalFis < 0.01)#filter for 0.01 cutoff
sigGO <- names(pvalFis[1:pval.01])

sigGeneTerm <- list()
enrichTable <- as.data.frame(matrix(rep(0, 6 * pval.01), nrow = pval.01))
colNames=c("RowNum", "GO.ID", "GO.Term","p.value", "GeneCount", "GenePerc")
colnames(enrichTable) <- colNames                             
for(i in 1:pval.01) {
  tmpGO <- grep(sigGO[i], goNames)
  genesTerm <- anno.genes[[tmpGO]]
  tGenes <- intersect(s24GenesEnt, genesTerm)
  sigGeneTerm[[i]] <- genesTerm
  enrichTable$RowNum[i] <-i
  enrichTable$GO.ID[i] <-sigGO[i]
  enrichTable$GO.Term[i] <- GOTERM[[sigGO[i]]]@Term
  enrichTable$p.value[i] <- unname(pvalFis[i])
  enrichTable$GeneCount[i] <- length(tGenes)
  enrichTable$GenePerc[i] <- length(tGenes)/length(anno.genes[[tmpGO]])
}
enrichTable1 <- enrichTable[enrichTable$GeneCount >= 3, ]

# gene bp
geneBP <- names(geneID2GO)[unlist(lapply(s24GenesEnt, grep, geneID2GO))]
geneTerms <- unlist(lapply(geneBP, function(x) GOTERM[[x]]@Term))

## bifido and s24-7
commonNeighbors <- intersect(s24Neighbors, bifidoNeighbors)

#shared neighbours
finalSharedNeighb <- data.frame(V(mergedNwCur)$name)
for (i in 1:length(V(mergedNwCur)$name)) {
  tmpSharedNeighb <- NULL
  tmpPerc <- NULL
  tmpNeigh1 <- ego(mergedNwCur, order = 1, nodes = i)[[1]]$name
  for (j in 1:length(V(mergedNwCur)$name)) {
    tmpNeigh2 <- ego(mergedNwCur, order = 1, nodes = j)[[1]]$name
    tmpSharedNeighb[j] <- length(intersect(tmpNeigh1, tmpNeigh2))
    tmpPerc[j] <- length(intersect(tmpNeigh1, tmpNeigh2))/length(tmpNeigh2)
  }
  finalSharedNeighb <- cbind.data.frame(finalSharedNeighb, tmpSharedNeighb, tmpPerc)
}

library(rentrez)
library(igraph)
x <- entrez_search(db="pubmed", term="Akkermansia AND S24-7")

mergedNwDF <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/merged5NwsEdgesDF.txt", as.is = T)
load(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/merged5Nws.RData")
verticesData <- get.data.frame(mergedNw, what = "vertices")
verticesData <- verticesData[, c(21, 22)]
verticesData <- verticesData[grep("Rik", verticesData$name, invert = T), ]
geneNames <- verticesData$name[verticesData$data == "gene"]
refHumanMouse <- read.delim(file = "/home/user/allInfoDir/workspace/conversion/refHumanMouse.txt", as.is = T)
microbNames <- verticesData$name[verticesData$data == "microb"]
metabNames <- verticesData$name[c(which(verticesData$data == "metabAS"), which(verticesData$data == "metabACU"))]
metabNames <- gsub("^X", "", metabNames)
cytoNames <- verticesData$name[verticesData$data == "cyto"]

removeDF <- mergedNwDF[grep("Rik", mergedNwDF$from), ]
removeDF <- rbind.data.frame(removeDF, mergedNwDF[grep("Rik", mergedNwDF$to), ])
tmpmergedNwDF <- mergedNwDF[intersect(grep("Rik", mergedNwDF$to, invert = T), grep("Rik", mergedNwDF$from, invert = T)), ]
tmpmergedNwDF$from <- gsub("^X", "", tmpmergedNwDF$from)
tmpmergedNwDF$to <- gsub("^X", "", tmpmergedNwDF$to)
write.table(tmpmergedNwDF, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/mergedNwDFCurated.txt", sep = "\t")
tmpmergedNwDF <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/mergedNwDFCurated.txt", as.is = T)
tmpmergedNwDF$from <- gsub("\\.g$", "", tmpmergedNwDF$from)
tmpmergedNwDF$to <- gsub("\\.g$", "", tmpmergedNwDF$to)
tmpmergedNwDF$from <- gsub("S24.7", "S24-7", tmpmergedNwDF$from)
tmpmergedNwDF$to <- gsub("S24.7", "S24-7", tmpmergedNwDF$to)


mergedNwDF <- tmpmergedNwDF
# mergedN
pubmedRes <- list()
absCount <- list()
pids <- list()
for (i in 1:dim(mergedNwDF)[1]) {
  tmpPair <- paste(mergedNwDF[i, ], collapse = " AND ")
  pubmedRes[[i]] <- entrez_search(db="pubmed", term = tmpPair)
  if (pubmedRes[[i]]$count != 0 && pubmedRes[[i]]$count <= 20) {
    pids[[i]] <- pubmedRes[[i]]$ids
    absCount[[i]] <- pubmedRes[[i]]$count
  }
  if (pubmedRes[[i]]$count != 0 && pubmedRes[[i]]$count > 20) {
    absCount[[i]] <- pubmedRes[[i]]$count
    pubmedRes[[i]] <- entrez_search(db = "pubmed", term = tmpPair, retmax = absCount[[i]])
    pids[[i]] <- pubmedRes[[i]]$ids
  } else {
    pids[[i]] <- NA
  }
}
names(pids) <- apply(mergedNwDF, 1, function(x) paste(x, collapse = " AND "))
save(pids, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/pidsMergedNwEntrezBactCor.RData") ## nothing new

save(pids, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/pidsMergedNwEntrez.RData")
load(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/pidsMergedNwEntrez.RData")
save(pubmedRes, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/pubmedRes.RData")
sum(is.na(pids))
connectionsPids <- data.frame(connections = names(pids)[!is.na(pids)], numPids = unname(unlist(lapply(pids[!is.na(pids)], length))), stringsAsFactors = F)
write.table(connectionsPids, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/pidsDF.txt", sep = "\t")

connectionsDataPids <- data.frame(connections = names(pids)[!is.na(pids)], data = NA, numPids = unname(unlist(lapply(pids[!is.na(pids)], length))), stringsAsFactors = F)
for (i in 1:dim(connectionsPids)[1]) {
  tmpConn <- connectionsPids$connections[i]
  tmpData <- NULL
  for (j in 1:2) {
    tmpData[j] <- V(mergedNw)$data[V(mergedNw)$name == unlist(strsplit(tmpConn, split = " AND "))[j]]
  }
  connectionsDataPids$data[i] <- paste(tmpData, collapse = " AND ")
}
write.table(connectionsDataPids, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/connectionsDataPids.txt", sep = "\t")
dataPids <- aggregate.data.frame(connectionsDataPids$numPids, by = list(connectionsDataPids$data), FUN = sum)
colnames(dataPids) <- c("LevelConnetions", "LiteratureHits")
write.table(dataPids, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/dataPids.txt", sep = "\t")

connectionsDataPids <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/connectionsDataPids.txt", as.is = T)
nodesFound <- unique(unlist(strsplit(split = " AND ", connectionsDataPids$connections)))
str(intersect(nodesFound, tableLevelHubs$Names)) ## from integrate.analyse

load(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/merged5Nws.RData")
nodesConn <- data.frame(nodeName = unique(sepConnectionsFound), numIds = rep(0, length(unique(sepConnectionsFound))), data = NA, stringsAsFactors = F)
for (i in 1:dim(nodesConn)[1]) {
  tmpNode <- nodesConn$nodeName[i]
  sumConn <- sum(connectionsPids$numPids[grep(tmpNode, connectionsPids$connections)])
  nodesConn$numIds[i] <- sumConn
  nodesConn$data[i] <- V(mergedNw)$data[V(mergedNw)$name == tmpNode]
}
write.table(nodesConn, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/nodesConnectionsData.txt", sep = "\t")

## bar plot of the data and connections
barplot(as.numeric(nodesConn[, 2]), names.arg = nodesConn$nodeName)
barplot(as.numeric(dataPids$LiteratureHits), names.arg = dataPids$LevelConnetions)
#
connectionsFound <- names(pids)[!is.na(pids)]
sepConnectionsFound <- unlist(strsplit(connectionsFound, split = " AND "))
pidsFound <- pids[!is.na(pids)]
save(pidsFound, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/pidsFound.RData")
load(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/pidsFound.RData")
write(connectionsFound, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/pairsFoundEntrez.txt", sep = "\t")
connectionsFound <- scan(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/pairsFoundEntrez.txt", sep = "\t", what = "", quote = "")

## for mapping to the edge table 
cytoscapeMergedNetworkEdge <- read.csv(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/mergedEdgeTableChg.csv", as.is = T)
edgeTableFound <- gsub("AND", "(undirected)", connectionsFound)
revEdge <- setdiff(edgeTableFound, cytoscapeMergedNetworkEdge$shared.name)
edgeTableFound[edgeTableFound %in% revEdge] <- gsub("(.*) (.*) (.*)", "\\3 \\2 \\1", revEdge)
dfEdgeFound <- data.frame(edgeTableFound, litFound = rep(1, length(edgeTableFound)))
write.table(dfEdgeFound, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/edgeTableConnctionsFound.txt", sep = "\t", quote = F, row.names = F)
## only glutathione pids
glutathionePids <- pidsFound[grep("Glutathione", names(pidsFound))]
for (i in 1:length(glutathionePids)) {
  write(glutathionePids[[i]], file = paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/", gsub(" ", "", names(glutathionePids)), ".txt")[i])
}
venn.diagram(glutathionePids[-6], category.names = gsub(" ", "", names(glutathionePids)[-6]), file = "/home/user/allInfoDir/samMiceWorkSpace/integration/glutathionePids.tiff")
mergedNwLit <- mergedNwDF
litConn <- rep(0, dim(mergedNwLit)[1])
for (i in 1:length(mergedNwLit$from)) {
  if (paste(mergedNwDF[i, ], collapse = " AND ") %in% connectionsFound) {
    litConn[i] <- 1
  }
}
mergedNwLit <- cbind.data.frame(mergedNwLit, litConn)
colnames(mergedNwLit) <- c("Source", "Target", "litConn")
write.table(mergedNwLit, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/literatureConnections.txt", sep = "\t", quote = F, row.names = F)
mergedNwLit <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/literatureConnections.txt", as.is = T)
tmpNw <- graph_from_data_frame(mergedNwLit)
V(tmpNw)$name <- V(tmpNw)$label
igraphObjGC <- set.graph.attribute(tmpNw, "name", paste0("litConn"))
cygraph <- toCytoscape(igraphObjGC)
test <- send2cy(cygraph, "default%white", 'force-directed') #MAKE A PROPER STYLE

## TRY FULLTEXT PACKAGE
library(fulltext)

dois <- ft_search(query = 'S24-7', from = c('plos', 'entrez', 'bmc'), bmc_key = "141bd19da7f82b05cf5b26fd01d7992e")
maxRet <- max(unlist(lapply(dois, function(x) x$found)))
dois <- ft_search(query = 'S24-7', from = c('plos', 'entrez', 'bmc'), bmc_key = "141bd19da7f82b05cf5b26fd01d7992e", limit = maxRet)
allDois <- c(dois$entrez$data$doi, dois$plot$data$id, dois$bmc$data$doi)
allDois <- unique(allDois[!is.na(allDois)])
ftOfDois1 <- ft_get(dois$entrez$data$doi[!is.na(dois$entrez$data$doi)], from = 'entrez')$entrez$data$data
ftOfDois2 <- ft_get(dois$plos$data$id[!is.na(dois$plos$data$id)], from = 'plos')$plos$data$data
ftOfDois3 <- ft_get(dois$bmc$data$doi)$bmc$data$data)
ftOfDois <- ft_get(allDois)

for (i in 1:length(allDois)) { 
  tmpFT <- ft_get(allDois[i])
xmlParsed <- xmlTreeParse(ftOfDois$plos$data$data[[i]])
xmlRootDoc <- xmlRoot(xmlParsed)
allDocData <- xmlSApply(xmlRootDoc, function(x) xmlSApply(x, xmlValue))
x <- grepl("microbe|Tumor Necrosis Factor", allDocData, ignore.case = T)
if (any(x))
  print(i)
}
#


library(bmc)
library(rplos)
library(rcrossref)
pubmedRes <- list()
absCount <- list()
pids <- list()
for (i in 1:dim(mergedNwDF)[1]) {
  tmpPair <- paste(mergedNwDF[i, ], collapse = " AND ")
  pubmedRes[[i]] <- entrez_search(db="pubmed", term = tmpPair)
  if (pubmedRes[[i]]$count != 0 && pubmedRes[[i]]$count <= 20) {
    pids[[i]] <- pubmedRes[[i]]$ids
    absCount[[i]] <- pubmedRes[[i]]$count
  }
  if (pubmedRes[[i]]$count != 0 && pubmedRes[[i]]$count > 20) {
    absCount[[i]] <- pubmedRes[[i]]$count
    pubmedRes[[i]] <- entrez_search(db = "pubmed", term = tmpPair, retmax = absCount[[i]])
    pids[[i]] <- pubmedRes[[i]]$ids
  } else {
    pids[[i]] <- NA
  }
}