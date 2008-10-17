testSet <- function(resBic, geneSetCol){
    
    nbBic <- as.integer(resBic$param[1, 2])
    nbSet <- length(geneSetCol)
    
    Data <- resBic$ExpressionSet
    expData <- exprs(Data)
    
    pval <- matrix(data = NA, nrow = nbBic, ncol = nbSet)
    inter <- matrix(data = NA, nrow = nbBic, ncol = nbSet)  
    setSize <- vector(length=nbSet)
    nbTest <- 0
    
    setPop <- GeneSet(Data, setName="pop")
    idTypePop <- geneIdType(setPop)
    dimPop <- length(geneIds(setPop))
    
    for (i in 1:nbSet){
        set <- geneSetCol[[i]]
        idTypeSet <- geneIdType(set)
        
        if (attributes(idTypeSet)$type != attributes(idTypePop)$type) setMapped <- mapIdentifiers(set, idTypePop)
        else setMapped <- set
        adjSet <- setPop & setMapped
        dimSet <- length(geneIds(adjSet))
        dimNonSet <- dimPop - dimSet
        
        for (j in 1:nbBic){
            bic <- bicluster(resBic, j, graph=FALSE)
            bic <- Data[rownames(bic),]
            setBic <- GeneSet(bic, setName="bic")
            dimBic <- length(geneIds(setBic))
            adjBic <- adjSet & setBic
            dimAdjBic <- length(geneIds(adjBic))
            
            if (dimAdjBic > 0){     ## mettre un seuil 
                pval[j, i] <- phyper(dimAdjBic, dimSet, dimNonSet, dimBic, lower.tail=FALSE)
            }
        }
    }
    
    adjpval <- mt.rawp2adjp(pval , proc="BY")
    adjpval <- matrix(data=adjpval$adjp[order(adjpval$index),2], nrow=nrow(pval), ncol=ncol(pval))
        
    res <- list(geneSetCol,pval,adjpval)
    names(res) <- c("geneSetCollection", "pvalues", "adjpvalues")
    class(res) <- "testSet"
    resBic$geneSet <- res
    return(resBic)
}
