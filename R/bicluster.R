bicluster <- function (biclustering, k, graph=TRUE)  
{
    if (class(biclustering) != "biclustering")
        stop ("invalid object class, must be of class biclustering")
    idRow <- which(biclustering$bicRow[k,] == 1)
    idCol <- which(biclustering$bicCol[k,] == 1)
    bic <- exprs(biclustering$ExpressionSet)[idRow, idCol]
    rownames(bic) <- featureNames(biclustering$ExpressionSet[idRow, idCol])
    colnames(bic) <- sampleNames(biclustering$ExpressionSet[idRow, idCol])
    class(bic) <- "bicluster"
    if (graph) plot(bic)
    return(bic)
}
