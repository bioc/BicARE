plot.bicluster <- function(x, ...){
    op <- par(mar = c(10,4,4,2))
    matplot(t(x), type = "b", pch = 1, ylab = "expression level", axes = FALSE, ...)
    axis(1, at = 1:ncol(x), labels = colnames(x), las = 2, cex.axis = 0.75)
    axis(2)
    par(op)
}

