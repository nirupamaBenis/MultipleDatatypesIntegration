#### this file contains some random useful functions, mostly from stackexchange

## a variation of load that allows you to set a name for the object being loaded
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

## opposite of intersect
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

##permute values in each row
randomPermuteByRows=function(xx){  
  yy <- NULL  ;for( i in 1:dim(xx)[1])  yy <- rbind(yy, xx[i,sample(1:dim(xx)[2])])
  colnames(yy) <- colnames(xx)
  rownames(yy) <- rownames(xx)
  return(yy)
}

##permute values in each col
randomPermuteByCols=function(xx){
  yy <- NULL  ;for( i in 1:dim(xx)[2])  yy <- cbind(yy, xx[sample(1:dim(xx)[1]),i])
  colnames(yy) <- colnames(xx)
  rownames(yy) <- rownames(xx)
  return(yy)
}

## takes several patterns and replacements to execute at once, length of pattern should be equal to length of replacement
gsubSeveral <- function (pattern, replacement, x) {
  tmpReplacement <- list()
 for (i in 1:length(pattern)) {
   tmpReplacement[[i]] = suppressWarnings(gsub(pattern[i], replacement[i], x[x == pattern[i]]))
 }
 res = unlist(tmpReplacement)
  return(res)
}

## use table function and order function simultaneously
orderedTable <- function (x, decreasing = T) {
  decreasing = T
  tableForm <- table(x)
  orderedTableForm <- tableForm[order(tableForm, decreasing = decreasing)]
  return(orderedTableForm)
}