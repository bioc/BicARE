makeReport <- function(dirPath, dirName, resBic, browse=TRUE){
    
    if(substr(dirPath, start=nchar(dirPath), stop=nchar(dirPath)) == "/") dirPath <- substr(dirPath, start=1, stop=nchar(dirPath)-1)
    
    outDir <- file.path(dirPath, dirName)
    dir.create(outDir, showWarnings=FALSE)
    BicAREPath <- file.path(attr(as.environment(match("package:BicARE", search())), "path"), "report")
    file.copy(c(list.files(BicAREPath, full.names=T)), outDir, overwrite=T)
    
    nbBic <- nrow(resBic$mat.resvol.bic)
    
    fileHome <- file.path(outDir, "home.html")
    
    contentHome <- readLines(fileHome, warn=FALSE)
    contentHome <- sub("DATE1", resBic$param[5,2], contentHome)
    contentHome <- sub("DATE2", date(), contentHome)
    contentHome <- sub("ANALYSISNAME", dirName, contentHome)
    contentHome <- sub("NBBIC", resBic$param[1,2], contentHome)
    contentHome <- sub("RESIDU", resBic$param[2,2], contentHome)
    contentHome <- sub("PROBAINIT", resBic$param[3,2], contentHome)
    contentHome <- sub("NBITE", resBic$param[4,2], contentHome)
    
    ## array type, if present
    Data <- resBic$ExpressionSet
    modPuce <- FALSE
    if ((length(annotation(Data))!=0) && (require(paste(annotation(Data), ".db", sep=""), character.only=TRUE))) modPuce <- TRUE
    
    ## biclusters
    fileListBic <- file.path(outDir, "listBic.html")
    contentListBic <- readLines(fileListBic, warn=FALSE)
    tab <- ""
    for (i in 1:nbBic){
        tabListGenes <- "<hr><br>\n<table>\n \t<tr><th></th><th>Symbol</th></tr>\n"
        tabListSamples <- "<hr><br>\n<table>\n \t<tr><th>Sample</th></tr>\n"
        
        pngName <- paste("bic", i, ".png", sep="")
        pngFile <- file.path(outDir, pngName)
        png(pngFile, width=600, height=600)
        bic <- bicluster(resBic, i, graph=TRUE)
        dev.off()
        fileBic <- file.path(outDir, paste("bic", i, ".html",sep=""))
        file.copy(file.path(outDir,"bic.html"), fileBic, overwrite=TRUE)
        contentBic <- readLines(fileBic)
        contentBic <- sub("biclusterName", paste("Bicluster", i, sep=" "), contentBic)
        contentBic <- sub("Graph", pngName, contentBic)
        
        listGenes <- rownames(bic)
        if (modPuce) symbol <- unlist(getSYMBOL(listGenes, annotation(Data)))  
        else symbol <- rep(NA, nrow(bic))
        for (j in 1:length(symbol)){
            tabListGenes <- paste(tabListGenes, "\t<tr> <td>", listGenes[j], "</td> <td>", symbol[j], "</td> </tr>\n", sep="")
        }
        tabListGenes <- paste(tabListGenes, "</table>\n", sep="")                                       
        contentBic <- sub("ListGenes", tabListGenes, contentBic)
        listSamples <- colnames(bic)
        
        for (j in 1:length(listSamples)){
            tabListSamples <- paste(tabListSamples, "\t<tr> <td>",listSamples[j], "</td> </tr>\n",sep="")
        }
        tabListSamples <- paste(tabListSamples, "</table>\n", sep="")                                       
        contentBic <- sub("ListSamples", tabListSamples, contentBic)
        
        write(contentBic, file=fileBic, append=FALSE)
        tab <- paste(tab, "\t<tr>\n\t\t<td><a href=\"bic", i, ".html\">bicluster ", i, "</a></td>\n",sep="")
        for (j in 1:5){
            tab <- paste(tab, "\t\t<td>", round(resBic$mat.resvol.bic[i,j], 3), "</td>\n", sep="")
        }
        tab <- paste(tab, "\t</tr>\n", sep="")
    }
    contentListBic <- sub("TAB", tab, contentListBic)
    write(contentListBic, file=fileListBic, append=FALSE)
    
    geneSet <- resBic$geneSet
    geneSetCol <- geneSet$geneSetCol
    fileGeneSet <- file.path(outDir, "geneSets.html")
    contentGeneSet <- readLines(fileGeneSet)
    if (is.null(geneSet)){
        contentGeneSet <- sub("content", "Analysis not performed", contentGeneSet)
        nbSets <- 0
    } else{
        nbSets <- length(geneSetCol)
        txt <- "\n<table>\n\t<tr>\n\t<th>Gene set</th>\n\t<th>Bicluster</th>\n\t<th>adj p-value</th>\n\t<th>p-value</th></tr>"
        index <- order(geneSet$adjpvalue, na.last=NA)
        
        index <- index[1:100]
        for (i in index){
            n <- floor((i-1)/nbBic)
            m <- round(nbBic * (i/nbBic - n),digits=0)
            nameSet <- setName(geneSetCol[[n+1]])
            
            txt2 <- paste("\n\t<tr>\n\t<td>",nameSet,"</td>\n\t<td><a href=\"bic", m,".html\">Bicluster ", m, "</a></td>\n\t<td>", geneSet$adjpvalue[i], "</td>\n\t<td>", geneSet$pvalue[i],"</td>\n\t</tr>", sep="") 
            txt <- paste(txt, txt2, sep="")
        }
        txt <- paste(txt, "\n</table>\n", sep="")
        contentGeneSet <- sub("content", txt, contentGeneSet)
    }
    contentHome <- sub("NBSETS", nbSets, contentHome)
    contentHome <- sub("ENRICHEDSET", sum(geneSet$adjpvalue<=0.2, na.rm=TRUE), contentHome)
    write(contentGeneSet, fileGeneSet, append=FALSE)
    
    covar <- resBic$covar
    fileCovar <- file.path(outDir, "covar.html")
    contentCovar <- readLines(fileCovar)
    if (is.null(covar)){
        contentCovar <- sub("content", "Analysis not performed", contentCovar)
        nbVar <- 0
    } else {
        nbVar <- length(covar$covar)
        txt <- ""
        txt2 <- "<ul>\n"
        for (i in 1:nbVar) {
            nameCov <-  names(covar$covar)[i]
            nbModa <- length(covar$covar[[i]])
            txt2 <- paste(txt2, "<li><a href=\"#", nameCov, "\">", nameCov,"</a></li>\n", sep="")
            txt <- paste(txt,"<hr><h3><a name=\"", nameCov, "\">", nameCov, "</a></h3>\n", sep="")
            txt <- paste(txt, "\n<table>\n\t<tr>\n\t<th></th>\n\t<th>adjp-value<br><i>(p-value)</i></th>\n", sep="")
            txt3 <- "\t<tr>\n\t<th></th>\n\t<th></th>\n"
            for (j in 1:nbModa){
                txt <- paste(txt, "\t<th>", names(covar$covar[[i]])[j], "</th>\n")
                txt3 <- paste(txt3, "\t<th>", covar$covar[[i]][j], "</th>\n",sep="")
            }
            txt <- paste(txt, "\t</tr>\n")
            txt <- paste(txt, txt3, "\t</tr>\n",sep="")
            pvalOrder <- order(covar$pvalue[,i])
            for (j in pvalOrder){
                txt <- paste(txt, "\t<tr>\n\t<td><a href=\"bic", j,".html\">Bicluster ", j, "</a></td>\n", sep="")
                txt <- paste(txt, "\t<td>", round(covar$adjpvalue[j,i], digits=4), "<br><i>(", round(covar$pvalue[j,i], digits=4), ")</i></td>\n", sep="")
                for (k in 1:nbModa){
                    txt <- paste(txt, "\t<td>", covar$index[[i]][j,k], "<br><i>",round(covar$residuals[[i]][j,k], digits=3),"</i></td>\n", sep="")
                }
                txt <- paste(txt, "\t</tr>\n", sep="")
            }
            txt <- paste(txt, "</table>\n", sep="")
        }
        txt2 <- paste(txt2, "</ul>\n")
        txt <- paste(txt2, txt)
        contentCovar <- sub("content", txt, contentCovar)
    }
    write(contentCovar, file=fileCovar, append=FALSE)
    contentHome <- sub("NBANNOT", nbVar, contentHome)
    contentHome <- sub("ENRICHEDANNOT", sum(covar$adjpvalues<=0.2, na.rm=TRUE), contentHome)

    write(contentHome, file=fileHome, append=FALSE)

    if (browse) browseURL(fileHome, browser = getOption("browser"))
    
}
