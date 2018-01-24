## This script scrambles the data to distort the biological information and then performs the same type of integration as the correlation networks. This correlation is done 1000 times to gauge how much the correlation networks are based on biological information. The top 5% of the correlation coefficients (thresholds in the correlation networks) in the 1000 iterations are plotted with the corresponding correlation network values.

## load libraries
library(mixOmics)
library(igraph)
## load function to create the correlation matrix
source("~/corrMixSpls.r")

## read in information on the experimental design
animalNumbers <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/animalNamesOrderedProp.txt", as.is = TRUE, header = TRUE)
color.edge <- colorRampPalette(c("darkgreen", "green", "yellow", "red", "darkred"))

## read in data
microb.data <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/javiAnalysis/processedData/diffMicrobDataShorter.txt", as.is = T, na.strings = "")
microb.data <- microb.data[-(grep(remR, rownames(microb.data))), ] # remR = "Chloroplast", has too many zeros
geneExpr.data <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/selGeneExpr.txt", as.is = T)
cytokine.data <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/cytokineData/selCytoLimma1.txt", as.is = T)
# metab.dataAU <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/Metabolomics/selAULimma1.txt", as.is = T)
metab.dataAS <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/Metabolomics/selASLimma1.txt", as.is = T)
metab.dataACU <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/Metabolomics/selACULimma1.txt", as.is = T)

allDataList <- list(microb.data, metab.dataAS, metab.dataACU, geneExpr.data, cytokine.data)
names(allDataList) <- c("microb", "metabAS", "metabACU", "gene", "cyto")
fileNamesMix <- scan(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/fileNamesAllMixMean.txt", what = "")
##### the list object allDataList has 5 datasets that need to be integrated. The names of the list members must correspond to the names in the vector fileNamesMix. eg., for integrating two datasets, 'microb' and 'metabAS', the variable fileNamesMix must have the string 'microb.metabAS' to get the for loop working. #####

## run the integration 1k times for each of the data set pairs
nwInfo <- data.frame(NetworkNames = fileNamesMix, cutoffPosRange = NA, cutoffNegRange = NA, NodesXRange = NA, NodesYRange = NA)
for (i in 1:length(fileNamesMix)) { 
  ## get the right members of the list of datasets
  data1 <- data.frame(t(allDataList[[grep(unlist(strsplit(fileNamesMix[i], split = "\\."))[1], names(allDataList))]]))
  data2 <- data.frame(t(allDataList[[grep(unlist(strsplit(fileNamesMix[i], split = "\\."))[2], names(allDataList))]]))
  ## make sure the columns are consistent in both datasets
  data1 <- data1[rownames(data1) %in% rownames(data2), ]
  data2 <- data2[rownames(data2) %in% rownames(data1), ]
  ## make an empty data.frame to hold important information about the networks in the different permutations
  networkInfoPermute <- data.frame(NetworkNames = rep(fileNamesMix[i], 1000), cutoffPos = NA, cutoffNeg = NA, NodesX = NA, NodesY = NA)
  ## get variables ready for the permutations
  tmpNw <- NULL
  cutoff1 <- NULL
  cutoff2 <- NULL
  for (k in 1:1000) { 
    ## permute by rows
    data1p <- randomPermuteByCols(data1)
    data2p <- randomPermuteByCols(data2)
    ## run the spls function from mixOmics
    tmpSpls <- spls(X = data1p, Y = data2p, ncomp = 5)
    ## correlate with an adaptation of a function in mixOmics
    tmpMix <- corrMixSpls(tmpSpls) # in new functions, taken from network function in mixomics
    ## look at correlations that are in top 5% on either side of the distribution
    cutoff1[k] <- unname(quantile(as.matrix(tmpMix), 0.95))
    cutoff2[k] <- unname(quantile(as.matrix(tmpMix), 0.05))
    ## use the cutoffs to limit the correlation matrix, anything that does not pass the threshold is made 0
    networkInfoPermute$cutoffPos[k] <- cutoff1
    networkInfoPermute$cutoffNeg[k] <- cutoff2
    filtTmpMix <- tmpMix
    for (j in 1:dim(tmpMix)[1]) {
      tmpNeg <- tmpMix[j, ]
      tmpNeg[tmpMix[j, ] >= 0] <- NA
      tmpNeg[tmpNeg > cutoff2[k]] <- 0
      tmpPos <- tmpMix[j, ]
      tmpPos[tmpMix[j, ] <= 0] <- NA
      tmpPos[tmpPos < cutoff1[k]] <- 0
      tmpRow <- unlist(ifelse(is.na(tmpNeg), tmpPos, tmpNeg))
      filtTmpMix[j, ] <- tmpRow
    }
    ## this correlation matrix is used to build a network
    tmpNw <- network(as.matrix(filtTmpMix), comp = 1:5, shape.node = c("rectangle", "rectangle"), color.node = c("deeppink", "blue"), threshold= (min(abs(cutoff1), abs(cutoff2)) - 0.1))
    ## set aside some information on the resulting network
    networkInfoPermute$NodesX[k] <- length(grep("X", V(tmpNw$gR)$name))
    networkInfoPermute$NodesY[k] <- length(grep("Y", V(tmpNw$gR)$name))
  }
  ## write down all the thresholds used in the 1000 iterations
  write(cutoff1, file = paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/", fileNamesMix[i], "PosCutoff.txt"), sep = "\t")
  write(cutoff2, file = paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/", fileNamesMix[i], "NegCutoff.txt"), sep = "\t")
  ## plot all the cutoff values
  jpeg(file = "/home/user/allInfoDir/samMiceWorkSpace/integration/randomCutoffs.jpg")
  par(mfrow = c(2,2))
  hist(cutoff2, freq = T, main = paste(fileNamesMix[i], "Negative Cutoff"))
  hist(cutoff1, freq = T, main = paste(fileNamesMix[i], "Positive Cutoff"))
  dev.off()
  ## add all the information to a data.frame to have ranges for all the dataset combinations
  nwInfo$cutoffPosRange[i] <- paste(range(networkInfoPermute$cutoffPos), collapse = " to ")
  nwInfo$cutoffNegRange[i] <- paste(range(networkInfoPermute$cutoffNeg), collapse = " to ")
  nwInfo$NodesXRange[i] <- paste(range(networkInfoPermute$NodesX), collapse = " to ")
  nwInfo$NodesYRange[i] <- paste(range(networkInfoPermute$NodesY), collapse = " to ")
}

## make the figures used in the paper
## names to be used in the figure
properNames <- c("Metabolomics Serum and Microbiota", "Metabolomics Urine and Microbiota", "Microbiota and Transcriptomics", "Microbiota and Cytokines", "Transcriptomics and Cytokines", "Metabolomics Serum and Transcriptomics", "Metabolomics Urine and Transcriptomics", "Metabolomics Serum and Cytokines", "Metabolomics Urine and Cytokines", "Metabolomics Serum and Metabolomics Urine")
## iterate through each of the 10 combinations
for (i in 1:length(fileNamesMix)) {
  ## read in the thresholds of all the random iterations
  cutoff <- c(scan(file = paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/", fileNamesMix[i], "PosCutoff.txt"), sep = "\t"), scan(file = paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/", fileNamesMix[i], "NegCutoff.txt"), sep = "\t"))
  ## read in the network with the actual data and make a network
  filtMean <- read.table(file = paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/meanSpls/", fileNamesMix[i], "MeanFilt.txt"), as.is = T)
  vectorFiltMean <- unlist(filtMean)
  vectorFiltMean <- vectorFiltMean[vectorFiltMean != 0]
  ## get variables for the plot
  random <- cutoff
  real <- vectorFiltMean
  ## range of the x-axis for the plot
  xrange <- c(min(vectorFiltMean, cutoff), max(vectorFiltMean, cutoff))
  ## number of breaks in histogram
  lengthOutRnd <- 100
  lengthOutReal <- 50
  ## to scale the histogram counts
  Factor <- 10
  ## create and modify histograms and prepare arguments for the plot
  rr <- hist(random, xlim = xrange, breaks = seq(xrange[1], xrange[2], length.out = lengthOutRnd))
  xx <- hist(real, breaks = seq(xrange[1], xrange[2], length.out = lengthOutRnd))
  axisxx <- axis(2) #store axis info 
  xx$counts <- Factor*xx$counts  #scale reads
  
  ## final plot
  svg(paste0("/home/user/allInfoDir/samMiceWorkSpace/integration/nwImages/", fileNamesMix[i], "randomNetworkPosNeg.svg"))
  par(mar = c(5,5,2,5))
  plot(rr, col = "grey", xlim = xrange, main = properNames[i], ylab = "Frequecy", xlab = "Correlation Values")
  plot(xx, add = TRUE, col = "red")
  axis(4, at = axisxx*Factor, labels = axisxx)
  mtext(side = 4, line = 2, "Frequency", at = mean(axisxx * Factor))
  dev.off()
}