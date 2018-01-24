## function to generate a correlation matrix based on the a splsObject (output of function spls from mixOmics)
## requires splsObj and number of components used when generating the spls object
## performs correlation using pairwise Pearson correlation
## returns a correlation matrix

corrMixSpls <- function (splsObj, ncomp) {
  mat <- splsObj
  row.names = mat$names$X
  col.names = mat$names$Y
  comp = 1:ncomp
  keep.X = apply(abs(mat$loadings$X), 1, sum) > 0
  keep.Y = apply(abs(mat$loadings$Y), 1, sum) > 0
  row.names = row.names[keep.X]
  col.names = col.names[keep.Y]
  cord.X = cor(mat$X[, keep.X], mat$variates$X[, comp], use = "pairwise")
  cord.Y = cor(mat$Y[, keep.Y], mat$variates$X[, comp], use = "pairwise")
  mat = cord.X %*% t(cord.Y)
  return(mat)
}
