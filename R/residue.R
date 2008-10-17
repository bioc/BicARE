residue <- function(Data)
{
    if (class(Data)!= "ExpressionSet") Data <- new("ExpressionSet", exprs=Data)
    expData <- exprs(Data)
    vecData <- as.vector(t(expData))
    bic.row <- rep(1, nrow(expData))
    bic.col <- rep(1, ncol(expData))
    resid <- .C("printres", as.integer(nrow(expData)), as.integer(ncol(expData)), as.double(vecData), as.integer(bic.row), as.integer(bic.col), res=double(1),PACKAGE="BicARE")$res
    return(resid)
}
