#### build and save all the spls objects neccessary
## each spls calculation is done twice to allow each dataset to be both dependent and independent variable, just to be sure

## load library
library(mixOmics)

## read in experiment information
animalNumbers <- read.delim(file="/home/user/allInfoDir/samMiceWorkSpace/animalNamesOrderedProp.txt", as.is=T, header=T)

################ READ IN ALL DATA TYPES ###############################
microb.data <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/javiAnalysis/processedData/diffMicrobDataShorter.txt", as.is = T, na.strings = "")
microb.data <- microb.data[-(grep(remR, rownames(microb.data))), ] # remR = "Chloroplast", has too many zeros
geneExpr.data <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/selGeneExpr.txt", as.is = T)
cytokine.data <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/cytokineData/selCytoLimma1.txt", as.is = T)
# metab.dataAU <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/Metabolomics/selAULimma1.txt", as.is = T) ## decided not to use
metab.dataAS <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/Metabolomics/selASLimma1.txt", as.is = T)
metab.dataACU <- read.delim(file = "/home/user/allInfoDir/samMiceWorkSpace/Metabolomics/selACULimma1.txt", as.is = T)

####################### microb * amine serum #########################
microb.dataT <- data.frame(t((microb.data)))
metab.dataAST <- data.frame(t(metab.dataAS))
metab.dataAST <- metab.dataAST[rownames(metab.dataAST) %in% rownames(microb.dataT), ]
microb.dataT <- microb.dataT[rownames(microb.dataT) %in% rownames(metab.dataAST), ]
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(metab.dataAST)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(microb.dataT)))
## spls
microb.metabAS <- spls(X = microb.dataT, Y = metab.dataAST, ncomp = 5)
metabAS.microb <- spls(X = metab.dataAST, Y = microb.dataT, ncomp = 5)
## save results
save(microb.metabAS, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/microb.metabAS.RData")
save(metabAS.microb, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/metabAS.microb.RData")

####################### microb * acrn urine #########################
microb.dataT <- data.frame(t((microb.data)))
metab.dataACUT <- data.frame(t(metab.dataACU))
metab.dataACUT <- metab.dataACUT[rownames(metab.dataACUT) %in% rownames(microb.dataT), ]
microb.dataT <- microb.dataT[rownames(microb.dataT) %in% rownames(metab.dataACUT), ]
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(metab.dataACUT)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(microb.dataT)))
## spls
microb.metabACU <- spls(X = microb.dataT, Y = metab.dataACUT, ncomp = 5)
metabACU.microb <- spls(X = metab.dataACUT, Y = microb.dataT, ncomp = 5)
## save results
save(microb.metabACU, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/microb.metabACU.RData")
save(metabACU.microb, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/metabACU.microb.RData")

######################## microb * genes ##########################
microb.dataT <- data.frame(t(microb.data))
geneExpr.dataT <- data.frame(t(geneExpr.data))
geneExpr.dataT <- geneExpr.dataT[rownames(geneExpr.dataT) %in% rownames(microb.dataT),]
microb.dataT <- microb.dataT[rownames(microb.dataT) %in% rownames(geneExpr.dataT),]
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(microb.dataT)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(geneExpr.dataT)))
## spls
gene.microb <- spls(X = geneExpr.dataT, Y = microb.dataT, ncomp = 5) 
microb.gene <- spls(X = microb.dataT, Y = geneExpr.dataT, ncomp = 5)
## save results
save(gene.microb, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/gene.microb.RData")
save(microb.gene, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/microb.gene.RData")

######################## microb * cytokine ############################
microb.dataT <- data.frame(t(microb.data))
cytokine.dataT <- data.frame(t(cytokine.data))
cytokine.dataT <- as.matrix(cytokine.dataT[rownames(cytokine.dataT) %in% rownames(microb.dataT),])
microb.dataT <- as.matrix(microb.dataT[rownames(microb.dataT) %in% rownames(cytokine.dataT),])
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(microb.dataT)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(cytokine.dataT)))
## spls
cyto.microb <- spls(X = cytokine.dataT, Y = microb.dataT, ncomp = 5) 
microb.cyto <- spls(X = microb.dataT, Y = cytokine.dataT, ncomp = 5) 
## save results
save(cyto.microb, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/cyto.microb.RData")
save(microb.cyto, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/microb.cyto.RData")

######################## gene * cytokine ############################
geneExpr.dataT <- data.frame(t(geneExpr.data))
cytokine.dataT <- data.frame(t(cytokine.data))
geneExpr.dataT <- geneExpr.dataT[rownames(geneExpr.dataT) %in% rownames(cytokine.dataT),]
cytokine.dataT <- cytokine.dataT[rownames(cytokine.dataT) %in% rownames(geneExpr.dataT),]
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(geneExpr.dataT)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(cytokine.dataT)))
## spls
gene.cyto <- spls(X = geneExpr.dataT, Y = cytokine.dataT, ncomp = 5) 
cyto.gene <- spls(X = cytokine.dataT, Y = geneExpr.dataT, ncomp = 5) 
## save results
save(gene.cyto, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/gene.cyto.RData")
save(cyto.gene, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/cyto.gene.RData")

######################## gene * amine serum ############################
geneExpr.dataT <- data.frame(t(geneExpr.data))
metab.dataAST <- data.frame(t(metab.dataAS))
geneExpr.dataT <- geneExpr.dataT[rownames(geneExpr.dataT) %in% rownames(metab.dataAST),]
metab.dataAST <- metab.dataAST[rownames(metab.dataAST) %in% rownames(geneExpr.dataT),]
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(geneExpr.dataT)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(metab.dataAST)))
## spls
gene.metabAS <- spls(X = geneExpr.dataT, Y = metab.dataAST, ncomp = 5) 
metabAS.gene <- spls(X = metab.dataAST, Y = geneExpr.dataT, ncomp = 5) 
## save results
save(gene.metabAS, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/gene.metabAS.RData")
save(metabAS.gene, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/metabAS.gene.RData")

######################## gene * acrn serum ############################
geneExpr.dataT <- data.frame(t(geneExpr.data))
metab.dataACUT <- data.frame(t(metab.dataACU))
geneExpr.dataT <- geneExpr.dataT[rownames(geneExpr.dataT) %in% rownames(metab.dataACUT),]
metab.dataACUT <- metab.dataACUT[rownames(metab.dataACUT) %in% rownames(geneExpr.dataT),]
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(geneExpr.dataT)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(metab.dataACUT)))
## spls
gene.metabACU <- spls(X = geneExpr.dataT, Y = metab.dataACUT, ncomp = 5) 
metabACU.gene <- spls(X = metab.dataACUT, Y = geneExpr.dataT, ncomp = 5) 
## save results
save(gene.metabACU, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/gene.metabACU.RData")
save(metabACU.gene, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/metabACU.gene.RData")

######################## cytokine * amine serum ############################
cytokine.dataT <- data.frame(t(cytokine.data))
metab.dataAST <- data.frame(t(metab.dataAS))
cytokine.dataT <- cytokine.dataT[rownames(cytokine.dataT) %in% rownames(metab.dataAST),]
metab.dataAST <- metab.dataAST[rownames(metab.dataAST) %in% rownames(cytokine.dataT),]
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(cytokine.dataT)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(metab.dataAST)))
## spls
metabAS.cyto <- spls(X = metab.dataAST, Y = cytokine.dataT, ncomp = 5)
cyto.metabAS <- spls(X = cytokine.dataT, Y = metab.dataAST, ncomp = 5) #, keepX = rep(200, 5)
## save results
save(cyto.metabAS, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/cyto.metabAS.RData")
save(metabAS.cyto, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/metabAS.cyto.RData")

######################## cytokine * acrn urine ############################
cytokine.dataT <- data.frame(t(cytokine.data))
metab.dataACUT <- data.frame(t(metab.dataACU))
cytokine.dataT <- cytokine.dataT[rownames(cytokine.dataT) %in% rownames(metab.dataACUT),]
metab.dataACUT <- metab.dataACUT[rownames(metab.dataACUT) %in% rownames(cytokine.dataT),]
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(cytokine.dataT)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(metab.dataACUT)))
## spls
cyto.metabACU <- spls(X = cytokine.dataT, Y = metab.dataACUT, ncomp = 5)
metabACU.cyto <- spls(X = metab.dataACUT, Y = cytokine.dataT, ncomp = 5)
## save results
save(cyto.metabACU, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/cyto.metabACU.RData")
save(metabACU.cyto, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/metabACU.cyto.RData")

######################## amine serum * acrn urine ############################
metab.dataAST <- data.frame(t(metab.dataAS))
metab.dataACUT <- data.frame(t(metab.dataACU))
metab.dataAST <- metab.dataAST[rownames(metab.dataAST) %in% rownames(metab.dataACUT),]
metab.dataACUT <- metab.dataACUT[rownames(metab.dataACUT) %in% rownames(metab.dataAST),]
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(metab.dataAST)))
table(gsub("([[:alpha:]]+)[[:digit:]]{3}", "\\1", rownames(metab.dataACUT)))
## spls
metabAS.metabACU <- spls(X = metab.dataAST, Y = metab.dataACUT, ncomp = 5)
metabACU.metabAS <- spls(X = metab.dataACUT, Y = metab.dataAST, ncomp = 5)
## save results
save(metabAS.metabACU, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/metabAS.metabACU.RData")
save(metabACU.metabAS, file = "/home/user/allInfoDir/samMiceWorkSpace/integration/metabACU.metabAS.RData")
