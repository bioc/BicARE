print.biclustering <-  function(x, ...){
    param <- x$param
    resVol <- x$mat.resvol.bic
    cl <- x$Call
    cat("\nCall : \n")
    print(cl)
    cat("\nParameters :")
    cat("\n number of biclusters : ", param[1,2])
    cat("\n residu threshold : ", param[2,2])
    cat("\n gene initial probability : ", param[3,2])
    cat("\n sample initial probability : ", param[4,2])
    cat("\n number of iterations : ", param[5,2])
    cat("\n date : ", param[6,2])    
    cat("\nBiclusters : \n")
    print(resVol, quote=FALSE)
    if (!is.null(x$covar) || !is.null(x$geneSet)) cat("\nAdditionnal analyses :\n")
    if (!is.null(x$covar)) cat("\t Sample covariates over-representation\n")
    if(!is.null(x$geneSet)) cat("\t Gene sets over-representation\n")
}
