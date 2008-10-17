testAnnot <- function(resBic, annot=NULL, covariates="all"){
    
    if (class(resBic)!="biclustering") stop("Not a biclustering object")
    nbBic <- as.integer(resBic$param[1, 2])
    ow <- options("warn")
    bicCol <- resBic$bicCol
    
    Data <- resBic$ExpressionSet
    if (is.null(annot)){
        annot <- pData(phenoData(Data))
        if (is.null(annot)) stop("No covariates provided !")
    }
    annot <- as.data.frame(annot)
    
    if ((length(covariates) == 1) && (covariates == "all")){
        whichVar <- 1:ncol(annot)
    } else whichVar <- match(covariates, colnames(annot))
    whichVar <- whichVar[!is.na(whichVar)]
    nbVar <- length(whichVar)
    
    pval <- matrix(nrow=nbBic, ncol=nbVar)
    colnames(pval) <- colnames(annot)[whichVar]
    index <- vector(length=nbVar, mode="list")
    resid <- vector(length=nbVar, mode="list")
    adjpvalue <- pval
    
    listVar <- vector(length=nbVar, mode="list")
    
    for (i in 1:nbVar){
        covar <- annot[, whichVar[i]]
        listVar[[i]] <- as.vector(summary.factor(covar))
        names(listVar[[i]]) <- levels(covar)
        names(listVar)[i] <- colnames(annot)[i]
        resid[[i]] <- matrix(nrow=nbBic, ncol=length(listVar[[i]]))
        index[[i]] <- matrix(nrow=nbBic, ncol=length(listVar[[i]]))
        names(index)[i] <- names(resid)[i] <- colnames(annot)[i]
        p <- listVar[[i]]
        idCol <- whichVar[i]
        colnames(resid[[i]]) <- colnames(index[[i]]) <- levels(covar)
        
        
        for (j in 1:nbBic){
            samples <- which(bicCol[j,]==1)
            x <- as.vector(summary(annot[samples, idCol]))
            
            options(warn=-1)
            test <- chisq.test(x, p=p, rescale.p=TRUE, simulate.p.value=FALSE)
            options(ow)
            
            pval[j,i] <- test$p.value
            resid[[i]][j,] <- test$residuals
            index[[i]][j,] <- x
            
        }
        adjpval <- mt.rawp2adjp(pval[,i] , proc="BY")
        adjpvalue[,i] <- matrix(data=adjpval$adjp[order(adjpval$index),2], nrow=nrow(pval), ncol=1)
    }
    colnames(adjpvalue) <- colnames(annot)[whichVar]
    
    res <- list(listVar, pval, adjpvalue, index, resid)
    names(res) <- c("covar", "pvalues", "adjpvalues", "index", "residuals")
    resBic$covar <- res
    return(resBic)
}
