FLOC <- function(Data, k=20, pGene=0.5, pSample=pGene, r=NULL, N=8, M=6, t=500, blocGene=NULL, blocSample=NULL)
{
    if (class(Data) != "ExpressionSet") Data <- new("ExpressionSet", exprs=Data)
    expData <- exprs(Data)
    
    if (is.null(r)) r <- residue(expData)/10
    
    cl <- match.call()
    
    param <- matrix(nrow=6, ncol=2)
    param[,1] <- c("number of bicluster", "residu threshold", "genes initial probability", "samples initial probability","number of iterations", "date")
    param[,2] <- c(k, r, pGene, pSample, t, date())
    
    vecData <- as.vector(t(expData))
    nbGenes <- nrow(expData)
    nbSamples <- ncol(expData)
    
    if (!is.null(blocGene) || !is.null(blocSample)) k <- max(ncol(blocGene), ncol(blocSample), k)
    vecBicRow <- matrix(data=0, nrow=nbGenes, ncol=k)
    vecBicCol <- matrix(data=0, nrow=nbSamples, ncol=k)
    
    if (is.null(blocGene)) blocGene <- matrix(data=0, nrow=nbGenes, ncol=k)
    else vecBicRow[,1:ncol(blocGene)] <- blocGene
    vecBlocGene <- vecBicRow <- as.vector(vecBicRow)
    
    if (is.null(blocSample)) blocSample <- matrix(data=0, nrow=nbSamples, ncol=k)
    else vecBicCol[,1:ncol(blocSample)] <- blocSample
    vecBlocSample <- vecBicCol <- as.vector(vecBicCol)                 
    
    rand1 <- runif(k * nbGenes)
    rand2 <- runif(k * nbSamples)
    
    vecBicRow[which(rand1 < pGene)] <- 1
    vecBicCol[which(rand2 < pSample)] <- 1
    
    vec.resvol.bic <- vector(mode="double", length = k * 4)
    
    res.floc <- .C("floc", vecData = as.double(vecData), as.integer(nbGenes), as.integer(nbSamples), vecBicRow = as.integer(vecBicRow), vecBicCol = as.integer(vecBicCol), vec.resvol.bic = as.double(vec.resvol.bic), as.double(r), as.integer(k),as.integer(N), as.integer(M), as.integer(t), as.integer(vecBlocGene), as.integer(vecBlocSample), PACKAGE="BicARE")
    
    bicRow <- matrix(data = res.floc$vecBicRow, nrow = k, ncol = nbGenes, byrow = TRUE)
    bicCol <- matrix(data = res.floc$vecBicCol, nrow = k, ncol = nbSamples, byrow = TRUE)
    
    mat.resvol.bic <- matrix(data = 0, nrow = k, ncol = 5)
    mat.resvol.bic[,1:4] <- matrix(data=res.floc$vec.resvol.bic, nrow = k, ncol = 4, byrow = TRUE)[,1:4]
    for (i in 1:k) mat.resvol.bic[i,5] <- mean(apply(expData[which(bicRow[i,] == 1), which(bicCol[i,] == 1)], 1, var))
    
    colnames(mat.resvol.bic) <- c("residue", "volume", "genes", "conditions", "rowvar")
    
    resBic <- list(cl, Data, param, bicRow, bicCol, mat.resvol.bic)
    names(resBic) <- c("Call", "ExpressionSet", "param", "bicRow", "bicCol", "mat.resvol.bic")
    
    class(resBic) <- "biclustering"
    return(resBic)
}
