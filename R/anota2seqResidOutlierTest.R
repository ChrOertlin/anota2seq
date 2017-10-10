anota2seqResidOutlierTest <- function(Anota2seqDataSet = NULL, confInt=0.01, iter=5, 
                                      generateSingleGenePlots=FALSE, nGraphs=200, 
                                      generateSummaryPlot=TRUE, residFitPlot=TRUE, useProgBar=TRUE){
    ##########################
    #########Input param checks
    if(is.null(Anota2seqDataSet)){
        stop("Please provide an Anota2seqDataSet.\n")
    }
    if(class(Anota2seqDataSet)!= "Anota2seqDataSet"){
        stop("Please provide an Anota2seqDataSet.\n")
    }
    if(is.null(Anota2seqDataSet@qualityControl)){
        stop("No qualityControl slot found in Anota2seqDataSet.\nPlease run the anota2seqPerformQC function on the Anota2seqDataSet before running anota2seqResidOutlierTest.\n")
    }
    if(is.null(confInt)){
        stop("Please provide confInt parameter. Must a numeric value between 0 and 1.\n")
    }
    if(!class(confInt)%in%c("numeric")){
        stop("confInt parameter must be a numeric value between 0 and 1.\n")
    }
    if(is.null(iter)){
        stop("Please provide iter parameter. Must a numeric value of 1 and greater.\n")
    }
    if(!class(iter)%in%c("numeric")){
        stop("iter parameter must be a numeric value of 1 and greater.\n")
    }
    if(is.null(generateSingleGenePlots)){
        stop("Please provide generateSingleGenePlots parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!generateSingleGenePlots%in%c(TRUE,FALSE)){
        stop("generateSingleGenePlots parameter must be set to TRUE or FALSE.\n")
    }
    if(is.null(nGraphs)){
        stop("Please provide nGraphs parameter. Must a numeric value of 1 and greater.\n")
    }
    if(!class(nGraphs)%in%c("numeric")){
        stop("nGraphs parameter must be a numeric value of 1 and greater.\n")
    }
    if(is.null(generateSummaryPlot)){
        stop("Please provide generateSummaryPlot parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!generateSummaryPlot%in%c(TRUE,FALSE)){
        stop("generateSummaryPlot parameter must be set to TRUE or FALSE.\n")
    }
    if(is.null(residFitPlot)){
        stop("Please provide residFitPlot parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!residFitPlot%in%c(TRUE,FALSE)){
        stop("residFitPlot parameter must be set to TRUE or FALSE.\n")
    }
    if(is.null(useProgBar)){
        stop("Please provide useProgBar parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!useProgBar%in%c(TRUE,FALSE)){
        stop("useProgBar parameter must be set to TRUE or FALSE.\n")
    }
    
    ##########################
    ##Get data
    residualMatrix <- Anota2seqDataSet@qualityControl@residuals
    nData <- dim(residualMatrix)[1]
    geneNames <- rownames(residualMatrix)
    ######################################################
    ##Warnings
    if(is.null(residualMatrix)){
        message("ERROR:No residuals as input\n")
        stop()
    }
    ######################################################
    ##Single gene plotting
    if(generateSingleGenePlots==TRUE){
        pdf("ANOTA2SEQ_residual_distributions_single.pdf", width=9, height=12)
        par(mfrow=c(4,3))
    }
    ######################################################
    ##Data structures
    residualMatrixOutlierSum <- matrix(ncol=iter, nrow=nData, dimnames=list("rownames"=rownames(residualMatrix)))
    residualMatrixOutlier <- matrix(ncol=dim(residualMatrix)[2], nrow=nData, dimnames=list("rownames"=rownames(residualMatrix), "colnames"=colnames(residualMatrix)))
    allSortScaleTrue <- allXs <- matrix(nrow=dim(residualMatrix)[2], ncol=nData)
    residualMatrixOutlierSumP <- rep(NA, iter)
    ##############################################
    ##Run analysis using the set n iterations
    message("Running anota2seqResidOutlierTest\n")
    total <- iter
    if(useProgBar==TRUE){
        pb <- txtProgressBar(min=0, max=total, style=3)
    }
    for(j in 1:iter){
        if(useProgBar==TRUE){
            setTxtProgressBar(pb, (j))
        }
        ##generate the random normally distributed data. 
        rnormMat <- matrix(data=rnorm(dim(residualMatrix)[2]*((1/confInt)-1)),
                           nrow=dim(residualMatrix)[2], ncol=((1/confInt)-1))
        
        ##Calculate the upper and lower limits of the data and scale
        rnormMat <- apply(scale(rnormMat),2,sort)
        env <- t(apply(rnormMat, 1, range))
        for(i in 1:nData){            
            ##Scale true data per gene sort and cbind to rnorm data set
            trueVec <- sort(scale(residualMatrix[i,]))
            sampMat <- cbind(trueVec, rnormMat)
            ##get real data set qq
            rs <- sampMat[,1]
            xs <- qqnorm(rs, plot=FALSE)$x
            ##get range of the sampled distribution sort position
            ##Calculate if obtained residuals falls outside expected from rnorm
            rsLog <- rs<env[,1] | rs>env[,2]
            residualMatrixOutlierSum[i,j] <- sum(rsLog)
            residualMatrixOutlier[i,] <- rsLog
            ##save true data
            allSortScaleTrue[,i] <- trueVec
            allXs[,i] <- xs
            if(iter==j & generateSingleGenePlots==TRUE & nGraphs>=i){
                anota2seqResidOutlierPlot(xs=xs, rs=rs, env=env, geneName=geneNames[i])
            }
            
        }
        ##Collect data for single
        residualMatrixOutlierSumP[j] <- sum(residualMatrixOutlier>0)
    }
    message("\n\n")
    if(generateSingleGenePlots==1){
        dev.off()
    }
    ################################################
    ##calcaulte obtained expected
    ##only create full summary for last iteration
    residualMatrixOutlierLog <- residualMatrixOutlier>0
    residualMatrixOutlierSumP <- residualMatrixOutlierSumP/(nData*dim(residualMatrixOutlier)[2])
    obtVsExpected <- residualMatrixOutlierSumP[j]/confInt
    expected <- nData *dim(residualMatrixOutlier)[2] *confInt
    obtained <- sum(residualMatrixOutlierLog)
    ################################################
    outputList <- new("Anota2seqResidOutlierTest",
                      confInt = confInt,
                      inputResiduals = residualMatrix,
                      rnormIter = iter,
                      outlierMatrixLog = residualMatrixOutlierLog,
                      meanOutlierPerIteration = residualMatrixOutlierSumP,
                      obtainedComparedToExpected = obtVsExpected,
                      nExpected = expected,
                      nObtained = obtained)
    #################################################
    ##Plotting summary
    if(generateSummaryPlot==TRUE){
        jpeg("ANOTA2SEQ_residual_distribution_summary.jpeg", width=800, height=800, quality=100)
        anota2seqResidOutlierPlotAll(all=allSortScaleTrue, xsAll=allXs, env=env, obtained=obtained, expected=expected, obtRelExpected=obtVsExpected, confInt=confInt)
        dev.off()
    }
    ##plot fitted vs residuals
    if(residFitPlot==TRUE){
        jpeg("ANOTA2SEQ_residual_vs_fitted.jpeg", width=900, height=900, quality=100)
        par(mfrow=c(2,1))
        plot(x=as.vector(Anota2seqDataSet@qualityControl@fittedValues), y=as.vector(Anota2seqDataSet@qualityControl@residuals), ylab="residuals", xlab="Fitted values", main="Residual vs fitted values")
        dev.off()
    }
    
    Anota2seqDataSet@residOutlierTest <- outputList
    return(Anota2seqDataSet)   
}

