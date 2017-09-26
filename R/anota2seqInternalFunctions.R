
##################################### 
####################################
#Function used to display selected output stats in the show() methods
showSelectedOutput <- function(anota2seqDataSet,analysis){
    
    # get the selectedOutput class ...
    if(is.null(anota2seq.get.output.class(anota2seqDataSet,analysis,"selected")) == FALSE){
        outSelClass <- anota2seq.get.output.class(anota2seqDataSet,analysis,"selected")
        cat(paste("Selected output of ",analysis, " contains:\n",sep=""))
        if(outSelClass@useRVM == TRUE){
            for(cont in 1:dim(anota2seqDataSet@contrasts)[2]){
                cat(paste("\tfor contrast",cont,"\n",sep=" "))
                
                if(is.null(anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE))== FALSE){
                    outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE)
                    if(analysis == "mRNA abundance"){
                        if(anota2seqDataSet@mRNAAbundance@mRNASelect[[1]] ==TRUE &anota2seqDataSet@mRNAAbundance@mRNASelect[[2]] ==TRUE){
                            outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE)[[1]]
                        }
                        else{
                            outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE)
                        }
                    }
                    if(length(outSel) > 2){
                        up <- outSel[outSel[,"apvEff"] > 0,,drop = FALSE]
                        down <- outSel[outSel[,"apvEff"] <0,,drop = FALSE]
                        if(analysis =="buffering"){
                            up <- outSel[outSel[,"apvEff"] < 0,,drop = FALSE]
                            down <- outSel[outSel[,"apvEff"] >0,,drop = FALSE]  
                        }
                        cat(paste("\t\tTotal: ",nrow(outSel)," genes.\n",sep=""))
                        cat(paste("\t\tPositive log2FC: ",nrow(up)," genes.\n",sep=""))
                        cat(paste("\t\tNegative log2FC: ",nrow(down)," genes.\n",sep=""))
                    }
                    if(length(outSel)<2){
                        message("No significant output with set filtering criteria, try to change the thresholds...")
                    }
                }
                
                if(is.null(anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE))== TRUE){
                    outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE)
                    message("No significant output with set filtering criteria, try to change the thresholds...")
                }
            }
        }
        if(outSelClass@useRVM == FALSE){
            for(cont in 1:dim(anota2seqDataSet@contrasts)[2]){
                cat(paste("\tfor contrast",cont,"\n",sep=" "))
                if(is.null(anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE))== FALSE){
                    outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE)
                    if(analysis == "mRNA abundance"){
                        if(anota2seqDataSet@mRNAAbundance@mRNASelect[[1]] ==TRUE &anota2seqDataSet@mRNAAbundance@mRNASelect[[2]] ==TRUE){
                            outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE)[[1]]
                        }
                        else{
                            outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE)
                        }
                    }
                    if(length(outSel) > 2){
                        up <- outSel[outSel[,"apvEff"] > 0,,drop = FALSE]
                        down <- outSel[outSel[,"apvEff"] < 0,,drop = FALSE]
                        cat(paste("\t\tTotal: ",nrow(outSel)," genes.\n",sep=""))
                        cat(paste("\t\tPostive log2FC: ",nrow(up)," genes.\n",sep=""))
                        cat(paste("\t\tNegative log2FC: ",nrow(down)," genes.\n",sep=""))
                    }
                    if(length(outSel) <2){
                        message("No significant output with set filtering criteria, try to change the thresholds...")
                    }
                }
                if(is.null(anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE))== TRUE){
                    outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE)
                    message("No significant output with set filtering criteria, try to change the thresholds...")
                }
            }
        }
    }
}
###################################
###################################
#Function used to display regulatory mode output stats in the show() methods
showRegModeOutput <- function(anota2seqDataSet,regMode,analysis){
    # Check whether any analysis selSigGenes has been performed
    if(is.null(anota2seq.get.output.class(anota2seqDataSet,"translation","selected"))|
       is.null(anota2seq.get.output.class(anota2seqDataSet,"translated mRNA","selected"))|
       is.null(anota2seq.get.output.class(anota2seqDataSet,"buffering","selected"))|
       is.null(anota2seq.get.output.class(anota2seqDataSet,"total mRNA","selected"))){
        cat("No regulatory modes specified.\n")
    }
    # Display regModes gene number for transation, buffering and mRNA abundance for each contrast...
    if(anota2seqDataSet@selectedTranslation@regModes == TRUE){
        # get the selectedOutput class ...
        if(is.null(anota2seq.get.output.class(anota2seqDataSet,analysis,"selected")) == FALSE){
            outSelClass <- anota2seq.get.output.class(anota2seqDataSet,analysis,"selected")
            cat(paste("\nRegulatory mode output of ",analysis, " contains:\n",sep=""))
            if(outSelClass@useRVM == TRUE){
                for(cont in 1:dim(anota2seqDataSet@contrasts)[2]){
                    cat(paste("\tfor contrast",cont,"\n",sep=" "))
                    
                    if(is.null(anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE))== FALSE){
                        outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE)
                        
                        if(analysis == "mRNA abundance"){
                            if(anota2seqDataSet@mRNAAbundance@mRNASelect[[1]] ==TRUE &anota2seqDataSet@mRNAAbundance@mRNASelect[[2]] ==TRUE){
                                outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE)[[1]]
                            }
                            else{
                                outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,TRUE)
                            }
                        }
                        
                        outSel <- outSel[which(outSel[,"singleRegMode"] == regMode),]
                        if(length(outSel) > 2){
                            up <- outSel[outSel[,"apvEff"] > 0,,drop = FALSE]
                            down <- outSel[outSel[,"apvEff"] <0,,drop = FALSE]
                            if(analysis =="buffering"){
                                up <- outSel[outSel[,"apvEff"] < 0,,drop = FALSE]
                                down <- outSel[outSel[,"apvEff"] >0,,drop = FALSE]  
                            }
                            cat(paste("\t\tTotal: ",nrow(outSel)," genes.\n",sep=""))
                            cat(paste("\t\tPositive log2FC: ",nrow(up)," genes.\n",sep=""))
                            cat(paste("\t\tNegative log2FC: ",nrow(down)," genes.\n",sep=""))
                        }
                        if(length(outSel)<2){
                            cat("Total: 0 genes.\n")
                        }
                    }
                    
                }
            }
            if(outSelClass@useRVM == FALSE){
                for(cont in 1:dim(anota2seqDataSet@contrasts)[2]){
                    cat(paste("for contrast",cont,"\n",sep=" "))
                    if(is.null(anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE))== FALSE){
                        outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE)
                        if(analysis == "mRNA abundance"){
                            if(anota2seqDataSet@mRNAAbundance@mRNASelect[[1]] == TRUE &anota2seqDataSet@mRNAAbundance@mRNASelect[[2]] == TRUE){
                                outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE)[[1]]
                            }
                            else{
                                outSel <- anota2seq.get.output(anota2seqDataSet,analysis,"selected",cont,FALSE)
                            }
                        }
                        outSel <- outSel[which(outSel[,"singleRegMode"] == regMode),]
                        if(length(outSel) > 2){
                            up <- outSel[outSel[,"apvEff"] > 0,,drop = FALSE]
                            down <- outSel[outSel[,"apvEff"] < 0,,drop = FALSE]
                            cat(paste("Total: ",nrow(outSel)," genes.\n",sep=""))
                            cat(paste("Postive log2FC: ",nrow(up)," genes.\n",sep=""))
                            cat(paste("Negative log2FC: ",nrow(down)," genes.\n",sep=""))
                        }
                        if(length(outSel) <2){
                            cat("Total: 0 genes.\n")
                        }
                    }
                }
            }
        } 
    }
}
#################################
#################################
# S4methods input checks ...
s4MethodChecks <- function(object,selContrast,output,analysis,useRVM,getRVM,visualizeRegModes,plotToFile,myBw,inFunc){
    if(is.null(object)){
        stop("Please provide an anota2seqDataSet.\n")
    }
    if(inFunc%in%c("output","delta","thresholds")){
        if(is.null(selContrast)){
            stop("Please provide one contrast using the selContrast parameter.\n")
        }
        if(length(selContrast) > 1){
            stop("Please provide only one contrast.\n")
        }
        if(max(selContrast) > dim(object@contrasts)[2]){
            stop("The selected contrasts does not exist. selContrast must be a numeric vector with 1 or more contrasts.\n The values cannot be higher than the number of columns in the contrast matrix.")
        }
        if(selContrast%in%c(TRUE,FALSE) & class(selContrast) == "logical"){
            stop("The selected contrasts does not exist. selContrast must be a numeric vector with 1 or more contrasts.\n The values cannot be higher than the number of columns in the contrast matrix.")
        }
    }
    if(inFunc == "delta"){
        if(is.null(output) == TRUE){
            stop("Please provide the output parameter.\n Must be one of the following: full or selected.\n")
        }
        if(length(output) > 1){
            stop("output parameter must be set to either full or selected.\n")
        }
        if(!output %in% c("full","selected")){
            stop("output parameter must be set to either full or selected.\n")
        }
    }
    if(inFunc%in%c("anota2seqPlotFC","anota2seqPlotPvalues","anota2seqPlotGenes")){
        if(is.null(selContrast)){
            stop("Please provide one or more contrasts using the selContrast parameter.\n")
        }
        if(max(selContrast) > dim(object@contrasts)[2]){
            stop("One of the selected contrasts does not exist. selContrast must be a numeric vector with 1 or more contrasts.\n The values cannot be higher than the number of columns in the contrast matrix.")
        }
    }
    if(inFunc%in%c("delta","thresholds")){
        if(is.null(analysis) == TRUE){
            stop("Please provide the analysis parameter.\nMust be one of the following: translated mRNA, total mRNA, translation, buffering or mRNA abundance.")
        }
        if(!analysis %in% c("translated mRNA","total mRNA","translation","buffering")){
            stop("analysis parameter wrong.\nMust be one of the following: translated mRNA, total mRNA, translation or buffering.")
        }
    }
    if(inFunc == "anota2seqPlotGenes"){
        if(is.null(analysis) == TRUE){
            stop("Please provide the analysis parameter.\nMust be one of the following: translation or buffering. ")
        }
        if(!analysis %in% c("translated mRNA","total mRNA","translation","buffering", "mRNA abundance")){
            stop("analysis parameter wrong.\nMust be one of the following: translation or buffering. ")
        } 
    }
    if(inFunc == "anota2seqPlotPvalues"){
        if(is.null(myBw)){
            stop("Please provide the myBw parameter. This parameter is used in the density() function and corresponds to the bw (bandwidth) parameter.\n")
        }
        if(is.null(useRVM)){
            stop("Please provide useRVM parameter. Must be set to TRUE or FALSE.\n")
        }
        if(!useRVM%in%c(TRUE,FALSE)){
            stop("useRVM parameter must be set to TRUE or FALSE.\n")
        }
    }
    if(inFunc == "output"){
        if(is.null(getRVM)){
            stop("Please provide getRVM parameter. Must be set to TRUE or FALSE.\n")
        }
        if(!getRVM%in%c(TRUE,FALSE)){
            stop("getRVM parameter must be set to TRUE or FALSE.\n")
        }
        if(is.null(output) == TRUE){
            stop("Please provide the output parameter.\n Must be one of the following: full, selected or regModes.\n")
        }
        if(length(output) > 1){
            stop("output parameter must be set to either full, selected or regModes.\n")
        }
        if(!output %in% c("full","selected","regModes")){
            stop("output parameter wrong ... must be either full or selected.\n")
        }
        if(is.null(analysis) == TRUE){
            stop("Please provide the analysis parameter.\nMust be one of the following: translated mRNA, total mRNA, translation, buffering or mRNA abundance.")
        }
        if(!analysis %in% c("translated mRNA","total mRNA","translation","buffering", "mRNA abundance")){
            stop("analysis parameter wrong.\nMust be one of the following: translated mRNA, total mRNA, translation, buffering or mRNA abundance.")
        }
    }
    if(inFunc == "anota2seqPlotFC"){
        
        if(is.null(visualizeRegModes)){
            stop("Please provide a visualizeRegModes parameter.\nMust be one of the following: all, none, translation or buffering.\n")
        }
        if(!visualizeRegModes%in%c("all","none","translation","buffering")){
            stop("visualizeRegModes parameter wrong.\nMust be one of the following: all, none, translation or buffering.\n")
        }
    }
    if(inFunc%in%c("anota2seqPlotFC","anota2seqPlotGenes","anota2seqPlotPvalues")){
        if(is.null(plotToFile)){
            stop("Please provide the plotToFile parameter.Must be set to TRUE or FALSE.\n")
        }
        if(!plotToFile%in%c(TRUE,FALSE)){
            stop("plotToFile parameter must be set to TRUE or FALSE.\n")
        }
    }
}
#################################
#################################
# This function performs input parameter checks so that no non-sense parameters can be used
anota2seqCheckParameter <- function(normalize,dataType,transformation,filterZeroGenes,varCutOff,analysis,inFunc){
    if(inFunc == "dataset"){
        # normalize parater
        if(is.null(normalize)){
            stop("normalize parameter is NULL. normalize parameter must be set to TRUE or FALSE. \n")
        }
        
        if(!is.null(normalize)){
            if(normalize %in% c(TRUE,FALSE) == FALSE){
                stop("normalize parameter must be set to TRUE or FALSE. \n")
            }
        }
        # datatype parameter
        if(is.null(dataType)){
            stop("dataType parameter is NULL. dataType parameter must be set to RNAseq or microarray. \n")
        }
        
        if(!is.null(dataType)){
            if(dataType %in% c("RNAseq","microarray") == FALSE){
                stop("dataType parameter must be set to either RNAseq or microarray. \n " )
            }
        }
        if(dataType == "microarray" & normalize == TRUE){
            stop("Data coming from DNA-microarrays should be preprocessed by the user.\nWhile using DNA-microarray data set normalize to FALSE.\n")
        }
        
        #filterZeroGenes
        if(is.null(filterZeroGenes)){
            stop("filterZeroGenes parameter must be set to TRUE or FALSE.\n")
        }
        if(!is.null(filterZeroGenes)){
            if(filterZeroGenes %in% c(TRUE,FALSE) == FALSE){
                stop("filterZeroGenes parameter must be set to TRUE or FALSE.\n")
            }
        }
        
        if(dataType == "microarray" & filterZeroGenes == TRUE){
            stop("Data coming from DNA-microarrays should be preprocessed by the user. While using DNA-microarray data set filterZeroGenes to FALSE.\n")
        }
        
        #transformation
        if(is.null(transformation)){
            if(dataType == "RNAseq" & normalize == TRUE){
                stop("transformation parameter is NULL while dataType is RNAseq and normalize is TRUE.\n Please set transformation parameter to either rlog or TMM-log2. \n")
            }
        }
        if(!is.null(transformation)){
            if(transformation %in% c("TMM-log2","rlog") == FALSE){
                stop("transformation parameter must be either rlog or TMM-log2. \n")
            }
        }
        if(is.null(varCutOff) == FALSE){
            if(class(varCutOff)!= "numeric"){
                stop("Please provde a numeric varCutOff parameter.\n")
            }
        }
    }
    
    if(inFunc == "analysis"){
        if(is.null(analysis)){
            stop("analysis parameter is NULL, must be set a vector containing one or more of the following strings translated mRNA, total mRNA, translation or buffering.\n")
        }
        
        if(!is.null(analysis)){
            if(length(analysis)>4){
                stop("analysis parameter has length > 4, can only be 4 or less but at least 1.")
            }
            for(reg in 1:length(analysis)){
                if(analysis[reg] %in% c("translated mRNA","total mRNA","translation","buffering") == FALSE){
                    stop(paste("position ",reg," of the analysis parameter is wrong.\n",
                               "analysis parameter must be a vector containing one or more of the following strings translated mRNA, total mRNA, translation or buffering.\n",sep=""))
                }
            }
        }
    }
}
######################################################
######################################################
######################################################
#Used in anota2seqDataset constructors, anota2seqRun and anota2seqAnalyze this functions performs several input data checks
anota2seqCheckInput <- function(dataP=NULL,dataT=NULL,phenoVec=NULL,batchVec=NULL,contrasts=NULL,correctionMethod = NULL,inFunc="none"){
    ### Put in function CheckData and input
    if(is.null(dataT)){
        stop("No data for total mRNA specified; check dataT input\n")
    }
    if(is.null(dataP)){
        stop("No data for translated mRNA specified; check dataP input\n")
    }
    # Check dataP and dataT
    if(is.null(rownames(dataP)) | is.null(rownames(dataP))){
        stop("data P and dataT need non-numeric rownames.\n")
    }
    if(is.numeric(rownames( dataP)) == TRUE | is.numeric(rownames( dataT)) == TRUE){
        stop("dataP and dataT need non-numeric rownames.\n")
    }
    
    if(nrow(dataP) != nrow(dataT)){
        stop("Number of rows for dataT and dataP must be identical.\n")
    }
    if(ncol(dataP) != ncol(dataT)){
        stop("Number of columns for dataT and dataP must be identical.\n")
    }
    ### Check the range of the input data, get a clue wheter continuous or count data is being used as input...
    
    if (is.null(phenoVec)) {
        stop("No phenotypes specified. Please supply a phenoVec.\n")
    }
    
    if(length( phenoVec) != ncol( dataP) | length( phenoVec) != ncol( dataT)){
        
        stop("length(phenoVec) must correspond to number of columns in dataT or dataP.\n")
    }
    
    nPheno <- length(levels(as.factor(phenoVec)))
    if(nPheno < 2){
        stop("Only one sample class provided in phenoVec.\nanota2seq needs at least 2 sample classes to perform the analysis.\n")
    }
    phenoLev <- levels(as.factor(phenoVec))
    for(s in 1:nPheno){
        if(nPheno == 2){
            if(length(phenoVec[phenoVec== phenoLev[s]]) < 3){
                stop(paste("Sample class ",phenoLev[s]," has less than three samples.\nanota2seq needs at least 3 samples per sample class if there are only 2 sample classes.\n" ))
            }
        }
        if(nPheno >2){
            if(length(phenoVec[phenoVec== phenoLev[s]]) < 2){
                stop(paste("Sample class ",phenoLev[s]," has less than 2 samples.\nanota2seq needs at least 2 samples per sample class if there are more than 2 sample classes.\n" ))
            }
        }
    }
    
    if(inFunc == "fromMatrix"){
        if (is.null(batchVec) == FALSE) {
            if(length(batchVec) < length(phenoVec)){
                stop("Not all samples assigned to batchVec. length(batchVec) must correspond to number of column in dataP or dataT.\n")
            }
            if(length(batchVec) > length(phenoVec)){
                stop("More batches than samples provided. length(batchVec) must correspond to number of column in dataP or dataT.\n")
            }
            nBatch <- length(levels(as.factor(batchVec)))
            if(nBatch <2){
                stop("Less than 2 batch classes provided in batchVec. Provide more batches or remove batchVec.\n")
            }
            
        }
    }
    if(inFunc == "fromSE"){
        if (is.null(batchVec) == FALSE) {
            if(length(batchVec) < length(phenoVec)){
                stop("Not all samples assigned to batchVec. length(batchVec) must correspond to number of column in dataP or dataT.\n")
            }
            if(length(batchVec) > length(phenoVec)){
                stop("More batches than samples provided. length(batchVec) must correspond to number of column in dataP or dataT.\n")
            }
            nBatch <- length(levels(as.factor(batchVec)))
            if(nBatch <2){
                stop("Less than 2 batch classes provided in the SummarizedExperiment batch column. Please specify more batches or remove the batch column.\n")
            }
            
        }
    }
    if(!identical(rownames(dataP),rownames(dataT))){
        
        stop("Rownames of dataP and dataT do not follow the same order.\n")
    }
    if (is.null(contrasts) == FALSE) {
        if(!class(contrasts) == "matrix"){
            stop("custom contrasts need to be provided as a matrix.\nPlease check your contrast matrix.\n")
        }
        if (dim(contrasts)[2] != (nPheno - 1)) {
            if (dim(contrasts)[2] > (nPheno - 1)) {
                
                stop("Too many custom contrasts supplied.\nPlease check your contrast matrix.\n")
            }
            if (dim(contrasts)[2] < (nPheno - 1)) {
                
                stop("Too few custom contrasts supplied.\nPlease check your contrast matrix.\n")
            }
        }
        if(sum(apply(contrasts,2,sum)) != 0){
            stop("Sum for each column in the contrast matrix must be 0. Please check your contrast matrix.\n")
        }
        if(identical(rownames( contrasts),levels(as.factor( phenoVec))) == FALSE){
            stop("Contrast matrix rownames are wrong.\nCheck the anota2seqDataSet help for an example to build a custom contrast matrix.\n")
        }
        if(is.null(limma::nonEstimable(contrasts)) == FALSE){
            stop("contrast matrix is not full rank. Please check your contrast matrix.\n")
        }
    }
    # check phenoVec
    
    #### End contrasts
    if(is.null(correctionMethod)){
        stop("Please provide a correctionMethod.\ncorrectionMethod can be set to:\nBonferroni, Holm, Hochberg, SidakSS, SidakSD, BH, BY, ABH, TSBH or qvalue\n")
    }
    if(length(correctionMethod) > 1){
        stop("Please select only one correctionMethod. correctionMethod can be set to:\nBonferroni, Holm, Hochberg, SidakSS, SidakSD, BH, BY, ABH, TSBH or qvalue\n")
    }
    if(!correctionMethod%in% c("Bonferroni", "Holm", "Hochberg", "SidakSS", "SidakSD", "BH", "BY", "ABH", "TSBH","qvalue")){
        stop("correctionMethod not recognized. correctionMethod can be set to:\nBonferroni, Holm, Hochberg, SidakSS, SidakSD, BH, BY, ABH, TSBH or qvalue\n")
    }
    ### end put in function checkdata
}
######################################################################
######################################################################
anota2seqSimDfbs <- function(nData=2000, phenoVec=phenoVec, mode=mode, useProgBar=useProgBar){
    nSamples <- length(phenoVec)
    ##Create matrix for simulated data
    dataPMatrix <- dataTMatrix <- matrix(ncol=nSamples, nrow=nData, data=NA)
    rownames(dataPMatrix) <- rownames(dataTMatrix) <- c(1:nData)
    ##The analysis is performed over a range of realistic correlations which are defined by the sd of the covariate.
    ##We have shown that this setting does not matter for the outcome because we use standardized dfbeta
    sdVec <- c(31:40/10)
    corMat <- matrix(ncol=length(sdVec), nrow=nData)
    ##A few structures to collect outputs
    pDfbCollect <- matrix(ncol=6, nrow=length(sdVec))
    row.names(pDfbCollect) <- sdVec
    corMeanMedCollect <- matrix(ncol=2, nrow=length(sdVec))
    rownames(corMeanMedCollect) <- sdVec
    colnames(corMeanMedCollect) <- c("Mean_correlation", "Median_correlation")
    ######################################################
    ##Start the analysis
    total <- length(sdVec)
    if(useProgBar==TRUE){
        pb <- txtProgressBar(min=0, max=total, style=3)
    }
    ##Over all the selected sds
    for(k in 1:length(sdVec)){
        if(useProgBar==TRUE){
            setTxtProgressBar(pb, k)
        }
        ##P and T is generated per gene
        for(j in 1:dim(dataPMatrix)[1]){
            ##sd=4 is arbitrary selected but gives a reasonable spread of the data and the setting does not influence the result.
            dataP <- rnorm(n=nSamples, mean=0, sd=4)
            dataT <- c(rep(NA, length(dataP)))
            for(i in 1:length(dataP)){
                dataT[i] <- rnorm(n=1, mean=dataP[i], sd=sdVec[k])
            }
            dataPMatrix[j,] <- dataP
            dataTMatrix[j,] <- dataT
        }
        ##Calculate and store the correlation per gene
        corVec <- c(rep(NA, nData))
        for(l in 1:length(corVec)){
            corVec[l] <- cor(dataPMatrix[l,], dataTMatrix[l,])
        }
        corMat[,k] <- corVec
        ##Store mean and median correlations
        corMeanMedCollect[k,1] <- mean(corVec)
        corMeanMedCollect[k,2] <- median(corVec)
        ##Run the dfb analysis and store the data
        dfbOut <- anota2seqDfbsSummaryOnly(dataT = dataTMatrix, dataP=dataPMatrix, phenoVec = phenoVec, mode=mode)
        ##Store the obtained dfb summary
        pDfbCollect[k,] <- as.vector(unlist(dfbOut$dfbSummary))
        colnames(pDfbCollect) <- colnames(dfbOut$dfbSummary)
    }
    message("\n\n")
    ########################################
    ##Generate output for return
    outputList <- list(
        "Correlation"=corMeanMedCollect,
        "SimulatedDfbs" = pDfbCollect)
    return(outputList)
}
#############################################################################
#############################################################################
anota2seqDfbsSummaryFull <- function(lmDfb, mode, useDfbSim, filename, nDfbSimData, phenoVec, useProgBar=useProgBar){
    dsfSummary <- anota2seqGetSummaryDfb(lmDfb)
    
    ##Perform dfb simulation to get thresholds
    if(useDfbSim==TRUE){
        message("\tPerforming dfbetas simulation\n")
        ##perform simulation
        nT <- length(phenoVec)
        ##mode decides how the simulation should be performed
        dfbSimOut <- anota2seqSimDfbs(nData=nDfbSimData, phenoVec=phenoVec, mode=mode, useProgBar=useProgBar)
        ##Get the output from anota2seqSimDfbs
        simDfbThreshold <- unlist(dfbSimOut$SimulatedDfbs)
        ##create a plot that compares obtined to simulated
        simDfbThresholdMean <- apply(simDfbThreshold, 2, mean)
        dsfSummary <- as.vector(unlist(dsfSummary))
        names(dsfSummary) <- c("dfb>1", "dfb>2", "dfb>3", "dfb>2/sqrt(N)", "dfb>3/sqrt(N)", "dfb>3.5*iqr")
        dfbDifference <- simDfbThresholdMean-dsfSummary
        dfbDifference <- round(dfbDifference, digits=4)
        tmpDfbMat <- cbind(simDfbThresholdMean, dsfSummary)
        tmpDfbMax <- apply(tmpDfbMat, 1, max)
        tmpVec <- c()
        for(r in 1:6){
            tmpVec <- c(tmpVec, dsfSummary[r], simDfbThresholdMean[r])
        }
        spaceVec <- c(0,0,rep(c(1,0),5))
        colVec <- rep(c(1,2),6)
        dfbsNames = c("dfb>1", "dfb>2", "dfb>3", "dfb>2/sqrt(N)", "dfb>3/sqrt(N)", "dfb>3.5*iqr")
        pdf(file=filename, width=6, height=4, pointsize=1/300)
        plot(x=c(0,17), y=c(0, c(max(tmpDfbMax))+0.01), ylab="Proportion of data points", xlab="Cut off method", main="Proportion of outliers in regression assessment using dfbetas", pch="", xaxt="n")
        barplot(tmpVec, space=spaceVec, col=colVec, add=TRUE)
        legend(x=1, y=max(tmpVec-0.002), legend=c("Obtained", "Simulated"), fill=c(1,2))
        text(y=c(tmpDfbMax+0.002), x=c(1, 4, 7, 10, 13, 16), labels=dfbDifference) 
        dev.off()
    }    
    ##Create a plotting output if no simulation was selected
    if(useDfbSim==FALSE){
        message("\tNo dfbetas simulation is performed\n")
        dfbsNames <- c("dfb>1", "dfb>2", "dfb>3", "dfb>2/sqrt(N)", "dfb>3/sqrt(N)", "dfb>3.5*iqr")
        pdf(file=filename, width=6, height=4, pointsize=1/300)
        barplot(as.vector(unlist(dsfSummary)), names.arg=dfbsNames, cex.names=1, ylab="Proportion of data points", xlab="Cut off method", main="Proportion of outliers in regression assessment using dfbetas")
        dev.off()
    }
    return(dsfSummary)
}
#################################################################################
#################################################################################
anota2seqDfbsSummaryOnly <- function(dataT=NULL, dataP=NULL, phenoVec=NULL, mode=mode){
    ##Warnings
    if(is.null(dataT)){
        stop("No data for total samples\n")
    }
    if(is.null(dataP)){
        stop("No data for polysomal samples\n")
    }
    if(is.null(phenoVec)){
        stop("No phenotypes specified\n")
    }
    if(is.null(mode)){
        stop("No mode specified\n")
    }
    if(identical(rownames(dataT), rownames(dataP))==FALSE){
        stop("Polysomal and Total rownames do not follow the same order\n")
    }
    ################################################
    ##calculate n data points and factorise phenoVec
    nData <- dim(dataP)[1]
    phenoVec <- as.factor(phenoVec)
    ##structures to collect outputs
    lmDfb <- matrix(ncol=dim(dataT)[2], nrow=nData)
    ################################################
    ##Start analysis in a per gene loop
    for(i in 1:nData){
        ##lm model with interaction
        if(mode=="int"){
            tmpLm <- lm(dataP[i,]~dataT[i,]*phenoVec)
        }
        if(mode=="add"){
            tmpLm <- lm(dataP[i,]~dataT[i,]+phenoVec)
        }
        ##get dfbs for the slope i.e. in column 2
        tmpDfb <- dfbetas(tmpLm)
        lmDfb[i,] <- tmpDfb[,2]
    }
    #################################################
    ##Evaluate outputs
    dfbSummary <- anota2seqGetSummaryDfb(lmDfb)
    
    ##Create a return object
    dataOut <- list(
        "dfbSummary"=dfbSummary)
    return(dataOut) 
}
#################################################################################
#################################################################################
anota2seqGetSummaryDfb <- function(lmDfb){
    ##Evaluate outputs using different thresholds for dfbetas that have been suggested in the litterature.
    dfb1 <- lmDfb>1
    dfb2 <- lmDfb>2
    dfb3 <- lmDfb>3
    dfb2Sqrt <- lmDfb>(2/sqrt(dim(lmDfb)[2]))
    dfb3Sqrt <- lmDfb>(3/sqrt(dim(lmDfb)[2]))
    ##3.5X IQR. IQR is calculated per gene and each gene is tested separately but summarized across all.
    dfbIqr <- apply(lmDfb, 1, IQR)    
    dfbIqrTh <- 3.5*dfbIqr
    dfbIqrLog <- lmDfb>dfbIqrTh
    ##generate a summary output
    ##compare as a proportion of all tests    
    nTests <- ncol(lmDfb) * nrow(lmDfb)
    dfb1P <- sum(dfb1)/nTests
    dfb2P <- sum(dfb2)/nTests
    dfb3P <- sum(dfb3)/nTests
    dfb2SqrtP <- sum(dfb2Sqrt)/nTests
    dfb3SqrtP <- sum(dfb3Sqrt)/nTests
    dfbIqrLogP <- sum(dfbIqrLog)/nTests
    dsfSummary <- list(
        "Proportion data points with dfb>1"=dfb1P,
        "Proportion data points with dfb>2"=dfb2P,
        "Proportion data points with dfb>3"=dfb3P,
        "Proportion data points with dfb>2/sqrt(N)"=dfb2SqrtP,
        "Proportion data points with dfb>3/sqrt(N)"=dfb3SqrtP,
        "Proportion data points with dfb>3.5*iqr"=dfbIqrLogP)
    return(dsfSummary)
}
#####################################################################################
#####################################################################################
### Wrapper function for pre-processing and data checks...
anota2seqRNAseqPreProcessing <- function(dataP=NULL,dataT=NULL,
                                         transformation ="TMM-log2",
                                         filterZeroGenes=FALSE,
                                         normalize=TRUE){
    
    if(is.null(dataP) | is.null(dataT)){
        stop("please specify both dataP and dataT\n")
    }
    
    if(filterZeroGenes == TRUE){
        tmpNoZero <- anota2seqRemoveZeroSamples(dataP,dataT)
    }
    if(filterZeroGenes == FALSE){
        tmpNoZero <- list(dataP = dataP,dataT=dataT)
    }
    tmpNorm <- anota2seqNormalize(tmpdataP = tmpNoZero$dataP,
                                  tmpdataT = tmpNoZero$dataT,
                                  transformation = transformation,
                                  normalize=normalize)
    return(list(dataP= tmpNorm$dataP,dataT=tmpNorm$dataT))
}
#############################################################################
#############################################################################
#### functions called within the wrapper function.
anota2seqRemoveZeroSamples <- function(tmpdataP=NULL,tmpdataT=NULL){
    #Remove genes with 0 values in at least one sample
    message("Removing zero count genes from data ...\n")
    tmdataP <- cbind(tmpdataP,tmpdataT)
    tmdataPNoZero <- tmdataP[!(apply(tmdataP, 1, function(y) any(y == 0))),]
    message(paste(nrow(tmdataP)-nrow(tmdataPNoZero)," genes were removed due to having 0 counts in at least one sample.\n"))
    dataP <- tmpdataP[rownames(tmdataPNoZero),]
    dataT <- tmpdataT[rownames(tmdataPNoZero),]
    ### end remove zero samples
    return(list(dataP=dataP,dataT=dataT))
}
#################################################################################
#################################################################################
anota2seqNormalize <- function(tmpdataP=NULL,tmpdataT=NULL,transformation=NULL,normalize=NULL){
    ## make sure that if dataP and cdata have matching colnames they stay separated throughout normalisation...
    colnames(tmpdataP) <- paste(colnames(tmpdataP),"_tmpColumnPoly",sep="")
    colnames(tmpdataT) <- paste(colnames(tmpdataT),"_tmpColumnCyto",sep="")
    tmdataP <- cbind(tmpdataP,tmpdataT)
    
    if(normalize == FALSE){
        if(max(range(tmpdataP,na.rm = TRUE))> 100 | max(range(tmpdataT,na.rm = TRUE))> 100){
            colnames(tmpdataT) <- gsub("_tmpColumnCyto","",colnames(tmpdataT))
            colnames(tmpdataP) <- gsub("_tmpColumnPoly","",colnames(tmpdataP))
            
            stop(" Data range indicates a count scale and normalize was set to FALSE.\nanota2seq requires continuous scale for efficient analysis.\nPlease check the normalization/transformation options (more details in the vignette).") 
        }
        if(max(range(tmpdataP,na.rm = TRUE))< 100  & max(range(tmpdataT,na.rm = TRUE))< 100){
            colnames(tmpdataT) <- gsub("_tmpColumnCyto","",colnames(tmpdataT))
            colnames(tmpdataP) <- gsub("_tmpColumnPoly","",colnames(tmpdataP))
            
            returndataP <- tmpdataP
            returndataT <- tmpdataT
        }
    }
    
    if(normalize == TRUE){
        if(max(range(tmpdataP,na.rm = TRUE))> 100  & max(range(tmpdataT,na.rm = TRUE))> 100){
            if(transformation == "TMM-log2"){
                message("Count data will be normalized and transformed according to the TMM-log2 method.\n")
                tmpNorm <- limma::voom(edgeR::calcNormFactors(edgeR::DGEList(tmdataP)))$E
            }
            if(transformation == "rlog"){
                message("Count data will be normalized and transformed according to the rlog method.\n")
                tmpNorm <- DESeq2::rlog(tmdataP)
                rownames(tmpNorm) <- rownames(tmdataP)
            }
            if(!transformation == "TMM-log2" & !transformation == "rlog"){
                stop("Unknown transformation method specified. Check transformation parameter, must be set to either: rlog or TMM-log2\n")
                
            }
            dataTNorm <- tmpNorm[,colnames(tmpdataT)]
            dataPNorm <- tmpNorm[,colnames(tmpdataP)]
            
            colnames(dataTNorm) <- gsub("_tmpColumnCyto","",colnames(tmpdataT))
            colnames(dataPNorm) <- gsub("_tmpColumnPoly","",colnames(tmpdataP))
            
            returndataP <- dataPNorm
            returndataT <- dataTNorm
        }
        if(max(range(tmpdataP,na.rm = TRUE))< 100  & max(range(tmpdataT,na.rm = TRUE))< 100){
            stop("RNAseq data range indicates continuous scale of input data while normalize is set to TRUE.\nCheck your parameter input.\n")
        }
    }
    return(list(dataP=returndataP,dataT=returndataT))
    ## end function normalize
}
##################################################################################
##################################################################################
anota2seqFiltCheckVar <- function(tmpdataP=NULL,tmpdataT=NULL,varCutOff=NULL,phenoVec=NULL){
    #sampleClasses
    levPheno <- levels(as.factor(phenoVec))
    # collect genes to remove...
    removeGenes <- vector("character")
    if(!is.null(varCutOff)){
        # message("Removing low variance genes from data ...\n")
        for(sampleClasses in 1:length(levPheno)){
            #get per sampleClass data
            scdataP <- tmpdataP[,which(levPheno[sampleClasses] == phenoVec)]
            scdataT <- tmpdataT[,which(levPheno[sampleClasses] == phenoVec)]
            #Calc Variance
            dataPVar <- apply(scdataP,1,var)
            dataTVar <- apply(scdataT,1,var)
            
            #Get genes with variance low variance per mRNA type...
            lowVardataP <- rownames(scdataP[dataPVar < varCutOff,])
            lowVardataT <- rownames(scdataT[dataTVar < varCutOff,])
            removeGenes <- c(removeGenes,lowVardataP,lowVardataT)
        }
        removeGenes <- unique(removeGenes)
        
        
        tmpdataP <- tmpdataP[setdiff(rownames(tmpdataP),removeGenes),,drop=FALSE]
        tmpdataT <- tmpdataT[setdiff(rownames(tmpdataT),removeGenes),,drop=FALSE]
        message(paste(length(removeGenes)," genes were removed by variance filtering using a cutoff of: ",varCutOff,"\n",sep=''))
        
        if(dim(tmpdataP)[1] < 1 | dim(tmpdataT)[1] < 1){
            stop("varCutOff parameter is set too high. Too many genes removed, lower the varCutOff value.\n")
        }
    }
    
    ## Check if filtered data has no genes with variance == 0
    for(sampleClasses in 1:length(levPheno)){
        scdataP <- tmpdataP[,which(levPheno[sampleClasses] == phenoVec)]
        scdataT <- tmpdataT[,which(levPheno[sampleClasses] == phenoVec)]
        dataPVar <- apply(scdataP,1,var)
        dataTVar <- apply(scdataT,1,var)
        if(sum(dataPVar <= 0) != 0 | sum(dataTVar <= 0)){
            stop("The analysis cannot be performed for genes with variance (per condition per RNA type) == 0.
                 Please set a variance threshold > 0 (parameter: varCutOff).")
        }
    }
    
    
    return(list(dataP=tmpdataP,dataT=tmpdataT))
}
##################################################################################
##################################################################################
anota2seqPlotIntPvals <- function(intP, intPAdj, intRvmP, intRvmPAdj, useRVM, correctionMethod, fileStem){
    cAx <- 2
    cLab <- 2
    cMain <- 2
    pdf(paste(fileStem, "_interaction_p_distribution.pdf", sep=""), height=10, width=10, pointsize=1/300)
    par(mar=c(5,6,4,2))
    par(mfrow=c(2,2))
    
    plot(density(intP), main="Omnibus interaction p-values", xlab="p-value", cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    
    tmpMain <- paste("Omnibus interaction adjusted p-values", "(", correctionMethod, ")")
    plot(density(intPAdj), main=tmpMain, xlab="Adjusted p-value", cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    hist(intP,  main="Omnibus interaction p-values",breaks=c(0:40)/40, xlab="p-value",
         cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    
    hist(intPAdj,  main=tmpMain, breaks=c(0:40)/40, xlab="Adjusted p-value",
         cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    
    if(useRVM==TRUE){
        plot(density(intRvmP), main="Omnibus interaction RVM p-values",
             xlab="RVM p-value", cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
        
        tmpMain <- paste("Omnibus interaction adjusted RVM p-values", "(", correctionMethod, ")")
        plot(density(intRvmPAdj), main=tmpMain, xlab="Adjusted RVM p-value",
             cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
        
        hist(intRvmP,  main="Omnibus interaction RVM p-values",
             breaks=c(0:40)/40, xlab="RVM p-value", cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
        
        hist(intRvmPAdj,  main=tmpMain, breaks=c(0:40)/40, xlab="Adjusted RVM p-value",
             cex.axis=cAx, cex.lab=cLab, cex.main=cMain)
    }
    dev.off()    
}
##################################################################################
##################################################################################
anota2seqAdjustPvals <- function(pVals, correctionMethod){
    tmpAdj <- mt.rawp2adjp(as.vector(pVals), proc=correctionMethod)
    pAdj <- tmpAdj$adjp[order(tmpAdj$index),2]
    return(pAdj)
}
##################################################################################
##################################################################################
anota2seqAdjustPvalsQ <- function(pVals){
    tmpAdj <- qvalue(as.vector(pVals))
    pAdj <- as.vector(tmpAdj$qvalues)
    return(pAdj)
}
##################################################################################
##################################################################################
anota2seqPerformRVM <- function(MS, Df, residDf, residMS){
    ab <- anota2seqGetab(residMS, residDf)
    names(ab) <- c("a", "b")
    residMSRvm <- ((residDf*residMS)+2/ab[2])/(residDf+(2*ab[1]))
    residDfRvm <- residDf+2*ab[1]
    rvmFval <- MS / residMSRvm
    rvmP <- 1 - pf(rvmFval, Df, residDfRvm)
    return(list(residMSRvm=residMSRvm,
                residDfRvm=residDfRvm,
                rvmFval=rvmFval,
                rvmP=rvmP,
                ab=ab)
    )
}
##################################################################################
##################################################################################
anota2seqSlopeTest <- function(tmpLm, curSlope, analysis){
    tmpLmSum <- summary(tmpLm)
    if(analysis == "translation"){
        if(curSlope < 0 ){
            ##we are doing a one tailed test compared to the 2 tailed test in the output i.e. divide p by 2
            slopeP <- tmpLmSum$coefficients[2,4] / 2
        }
        if(curSlope>1){
            ##compare if slope is sig different to 1
            tmpSlopeEst <- tmpLmSum$coefficients[2,1] - 1
            tmpSlopeErr <- tmpLmSum$coefficients[2,2]
            tmpSlopeT <- tmpSlopeEst/tmpSlopeErr
            ##using resiual dfs to test p-value
            tmpSlopeDf <- tmpLm$df.residual
            slopeP <- 1-pt(tmpSlopeT, tmpSlopeDf)
        }
    }
    if(analysis == "buffering"){
        if(curSlope < -1 ){
            ##we are doing a one tailed test compared to the 2 tailed test in the output i.e. divide p by 2
            slopeP <- tmpLmSum$coefficients[2,4] / 2
        }
        if(curSlope > 0){
            ##compare if slope is sig different to 1
            tmpSlopeEst <- tmpLmSum$coefficients[2,1] - 1
            tmpSlopeErr <- tmpLmSum$coefficients[2,2]
            tmpSlopeT <- tmpSlopeEst/tmpSlopeErr
            ##using resiual dfs to test p-value
            tmpSlopeDf <- tmpLm$df.residual
            slopeP <- 1-pt(tmpSlopeT, tmpSlopeDf)
        }
    }
    return(slopeP)
}
##################################################################################
##################################################################################
anota2seqGetIntercepts <- function(x, y, slope, phenoLev, phenoVecOrg){
    tmpInt <- rep(NA, length(phenoLev))
    for(k in 1:length(phenoLev)){
        tmpX <- mean(x[phenoVecOrg==phenoLev[k]])
        tmpY <- mean(y[phenoVecOrg==phenoLev[k]])
        tmpInt[k] <- tmpY-(slope * tmpX)
    }
    return(tmpInt)
}
##################################################################################
##################################################################################
anota2seqPlotSingleRegression <- function(x, y, geneName, intercepts, slope, phenoVecOrg, phenoLev){
    plot(x=x, y=y, pch="", main=geneName)
    text(x=x, y=y, labels=phenoVecOrg)
    for(k in 1:length(phenoLev)){
        tmpMinX <- min(x, na.rm=TRUE)
        tmpMaxX <- max(x, na.rm=TRUE)
        lines(x=c(tmpMinX, tmpMaxX), y=c((intercepts[k] + slope*tmpMinX), (intercepts[k] + slope*tmpMaxX)),  lty=k)
    }
}
##################################################################################
##################################################################################
##################################################################################
##################################################################################
anota2seqPerformQcWarnings <- function(nPheno, phenoLev, phenoVecOrg, onlyGroup){
    ##require 3 samples per category unless using onlyGroup=TRUE
    for(i in 1:nPheno){
        if(sum(phenoVecOrg==phenoLev[i])<3){
            if(onlyGroup==FALSE){
                stop("If there are less than 3 sample classes anota2seq requires 3 samples per sample class (in phenoVec) to run\n\tIf there is at least two samples per group and more than two groups, the onlyGroup mode can be used to assess omnibus group effects\n")
            }
            if(onlyGroup==TRUE){
                if(sum(phenoVecOrg==phenoLev[i])==1){
                    stop("\tNo analysis can be run whith only one sample despite using the onlyGroup mode\n")
                }
                if(nPheno<3){
                    stop("\tThree sample classes is required to run the onlyGroup analysis when there is <3 samples per group")
                }
            }
        }
    }
}
#################################################################
#################################################################
anota2seqCalculateDeltaData <- function(dataP,dataT,phenoVec,useIds,contrasts,selContr){
    phenoLev <- levels(as.factor(phenoVec))
    deltaT <- matrix(ncol = length(phenoLev), nrow = length(useIds))
    deltaP <- matrix(ncol = length(phenoLev), nrow = length(useIds))
    for (i in 1:length(useIds)) {
        for (j in 1:length(phenoLev)) {
            if(!is.null(dataT)){
                deltaT[i, j] <- mean(dataT[useIds[i], phenoVec ==
                                               phenoLev[j]],drop=FALSE)
            }
            if(!is.null(dataP)){
                deltaP[i, j] <- mean(dataP[useIds[i], phenoVec ==
                                               phenoLev[j]],drop=FALSE)
            }
        }
    }
    deltaPT <- deltaP - deltaT
    deltaTP <- deltaT - deltaP
    rownames(deltaP) <- rownames(deltaT) <- rownames(deltaPT) <- rownames(deltaTP) <- useIds
    colnames(deltaP) <- colnames(deltaT) <- colnames(deltaPT) <- colnames(deltaTP) <- phenoLev
    
    return(cbind(deltaP = anota2seqCalculateFC(deltaP,useIds,contrasts,selContr),
                 deltaT = anota2seqCalculateFC(deltaT,useIds,contrasts,selContr),
                 deltaPT = anota2seqCalculateFC(deltaPT,useIds,contrasts,selContr),
                 deltaTP = anota2seqCalculateFC(deltaTP,useIds,contrasts,selContr)))
}
anota2seqCalculateFC <- function(tmpData,useIds,contrasts,selContr){
    tmp1 <- t(t(tmpData[useIds, ,drop=FALSE]) * contrasts[,
                                                          selContr])
    tmp2 <- apply(tmp1, 1, sum)
    return(tmp2)
}
##################################################################
##################################################################
#################################################################################
#################################################################################
anota2seqResidOutlierPlot <- function(rs=NULL, xs=NULL, env=NULL, geneName=""){  
    matplot(xs, cbind(rs, env), type="pnn", pch=4, mkh=0.06, axes=FALSE, xlab="", ylab="", main=geneName)
    ##gets the limits and calculate lengths for the bars that mark the boundaries of the envelope
    xyul <- par("usr")
    smidge <- min(diff(c(xyul[1], xs, xyul[2])))/2
    segments(xs-smidge, env[,1], xs+smidge, env[,1])
    segments(xs-smidge, env[,2], xs+smidge, env[,2])
    ##get axis with nice ticks
    xul <- trunc(10*xyul[1:2])/10
    axis(side=1, at=seq(xul[1], xul[2], by=0.1), labels=FALSE, tck=0.01)
    xi <- trunc(xyul[1:2])
    axis(side=1, at=seq(xi[1], xi[2], by=0.5), tck=0.02)
    yul <- trunc(5*xyul[3:4])/5
    axis(side=2, at=seq(yul[1], yul[2], by=0.2), labels=FALSE, tck=0.01)
    yi <- trunc(xyul[3:4])
    axis(side=2, at=yi[1]:yi[2], tck=0.02)
    box(bty="l")
    mtext("Quantiles of Standard Normal", side=1, line=2.5, font=3)
    mtext(expression(R[1]), side=2, line=2, at=yul[2])
}
#####################################################################
#####################################################################
anota2seqResidOutlierPlotAll <- function(all=NULL, xsAll=xsAll,  env=env, obtained, expected, obtRelExpected, confInt){
    allMax=max(all)
    allMin=min(all)
    xsMin=min(xsAll)
    xsMax=max(xsAll)
    ##Get number of outliers per rankposition
    allLog <- all<env[,1] | all>env[,2]
    allLogSum <- apply(allLog, 1, sum)
    allLogSumP <- 100*(allLogSum / dim(allLog)[2])
    allLogSumP <- round(allLogSumP, digits=3)
    ##plot
    plot(x=c(xsMin,xsMax), y=c(allMin-0.2, allMax+0.4), pch="", axes=TRUE, xlab="Quantiles of standard normal", ylab="R", main="Summary of all residuals")
    for(i in 1:dim(all)[2]){
        points(x=xsAll[,i], y=all[,i], pch=16, cex=0.2)
    }
    if(is.null(obtained)==FALSE){
        text(x=xsMin+0.4, y=allMax-0.4-(0.05*allMax), labels=paste("Expected outliers: ", confInt*100, "%", sep=""))
        text(x=xsMin+0.4, y=allMax-0.4-(0.1*allMax), labels=paste("Obtained outliers: ", round(obtRelExpected*confInt*100, digits=3), "%", sep=""))
    }
    segments(xsAll[,i]-0.05, env[,1], xsAll[,i]+0.05, env[,1], col=2, lwd=1.5)
    segments(xsAll[,i]-0.05, env[,2], xsAll[,i]+0.05, env[,2], col=2, lwd=1.5)
    text(x=xsAll[,1], y=env[,2]+0.2, labels=paste(allLogSumP, "%", sep=""))
}
#############################################################################
#############################################################################
anota2seqPlotIGFit <- function(useVar, df, title="Empirical", doQQ=TRUE, qqName=NULL) {
    ##Internal fuctions in anotaPlotIGFit
    ## Generate empirical probability density function for data.
    ## data - vector of data points
    ## q - trimming value.  remove 1-q points as outliers from greater tail.
    my.cum <- function(data, q=.9) {
        len <- getLen(data, quan=q)
        maxi <- sort(data)[len]
        x <- seq(min(data), maxi, length=len)
        x <- seq(min(data), maxi, length=len)
        p <- rep(0,len-1)
        lenny <- length(data)
        for(i in 1:len)
            p[i] <- (sum(data<x[i]))/lenny
        return(cbind(x, p))
    }
    #####################################
    flik <- function(p,y){
        ## log liklihood for a*b*x from an F distribution with m and 2*a degrees of freedom
        ## y is a vector containing data and the m values, p contains a and b.
        x<-y[1:(length(y)/2)]
        m<-y[(length(y)/2+1):length(y)]
        p<-abs(p)
        a<-p[1]
        b<-p[2]
        x<-x*(a*b)
        n<-2*a
        out<-base::log(df(x,m,n))+base::log(a*b)
        sum(-out)
    }
    ######################################
    ## Get rid of the very large data points so graphs scale better
    ## quan - trimming quantile
    getLen <- function(data, quan=.90) {
        return(trunc(quan*length(data)))
    }
    ######################################
    degFreedom1 <- df
    vars <- as.vector(useVar)
    ab <- anota2seqGetab(vars, rep(degFreedom1, length(vars)))
    message("\tThe a and b parameters for the inverse gamma distribution are:\n\t", paste(c("a:",ab[1], " b:", ab[2])), "\n")
    adj <- ab[1]*ab[2]
    adjVars <- vars*adj[1]
    scum <- my.cum(adjVars, q=0.9)
    probF <- pf(scum[,1], degFreedom1, 2*ab[1])
    lineWidth <- 2
    if(doQQ) {
        num <- length(adjVars)
        theoF <- rf(num, degFreedom1, df2=2*ab[1])
        qqplot(theoF, adjVars,xlab="Theoretical", ylab="Empirical", main=qqName)
        abline(0,1)
    }
    ##Turn off warnings temporary as there is always an informative warning message
    options(warn=(-1))
    kRes <- ks.test(x=adjVars, y="pf", degFreedom1, 2*ab[1])$p.value
    options(warn=0)
    plot(scum[,1], scum[,2], type="l", lwd=lineWidth, xlab="Var", ylab="cdf",
         main=paste(title, ": KS p-value=", signif(kRes,3), sep=""))
    lines(scum[,1], scum[,2], col=2)
    lines(scum[,1], probF, col=5, lwd=lineWidth)
    legend(x=c(1), y=c(0.3),  legend=c( title, "Theoretical F"), fill=c(2,5))
    options(warn=(1))
}
#############################################################################
#############################################################################
anota2seqGetab <- function(sig,n){
    flik <- function(p,y){
        ## log liklihood for a*b*x from an F distribution with m and 2*a degrees of freedom
        ## y is a vector containing data and the m values, p contains a and b.
        x<-y[1:(length(y)/2)]
        m<-y[(length(y)/2+1):length(y)]
        p<-abs(p)
        a<-p[1]
        b<-p[2]
        x<-x*(a*b)
        n<-2*a
        out<-base::log(df(x,m,n))+base::log(a*b)
        sum(-out)
    }
    set<-(!is.na(sig)&n>0&sig>0)
    sig<-sig[set]
    n<-n[set]
    set<-n>4
    if (sum(set)>0){
        m1<-(n[set]-2)/((n[set])*sig[set])
        m2<-(n[set]-2)*(n[set]-4)/((n[set])*sig[set])^2
        m1<-mean(m1,na.rm=TRUE)
        m2<-mean(m2,na.rm=TRUE)
        b<-m2/m1-m1
        a<-m1^2/(m2-m1^2)
    }
    else{ a<-b<-1}
    strt<-c(a,b)
    ###PATCH
    g <- function(p,yunq) flik(p,yunq)
    ########
    options(warn=(-1))
    a<-nlm(g,strt, yunq=c(sig,n))
    options(warn=0)
    a$estimate<-abs(a$estimate)
}
