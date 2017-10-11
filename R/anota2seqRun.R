anota2seqRun <- function(Anota2seqDataSet,contrasts = NULL, 
                         performQC = TRUE, onlyGroup = FALSE, performROT = TRUE, 
                         generateSingleGenePlots = FALSE, analyzeBuffering =TRUE, 
                         analyzemRNA = TRUE, thresholds = NULL, useRVM = TRUE, correctionMethod ="BH", useProgBar = TRUE){
    
    if(is.null(Anota2seqDataSet)){
        stop("Please provide an Anota2seqDataSet.")
    }
    if(is.null(performQC)){
        stop("Please provide performQC parameter. Must be either TRUE or FALSE.\n")
    }
    if(!performQC%in%c("TRUE","FALSE")){
        stop("performQC parameter must be either TRUE or FALSE.\n")
    }
    if(is.null(performROT)){
        stop("Please provide performROT parameter. Must be either TRUE or FALSE.\n")
    }
    if(!performROT%in%c("TRUE","FALSE")){
        stop("performROT parameter must be either TRUE or FALSE.\n")
    }
    if(is.null(generateSingleGenePlots)){
        stop("Please provide generateSingleGenePlots parameter. Must be either TRUE or FALSE.\n")
    }
    if(!generateSingleGenePlots%in%c("TRUE","FALSE")){
        stop("generateSingleGenePlots parameter must be either TRUE or FALSE.\n")
    }
    if(is.null(analyzeBuffering)){
        stop("Please provide analyzeBuffering parameter. Must be either TRUE or FALSE.\n")
    }
    if(!analyzeBuffering%in%c("TRUE","FALSE")){
        stop("analyzeBuffering parameter must be either TRUE or FALSE.\n")
    }
    if(is.null(analyzemRNA)){
        stop("Please provide analyzemRNA parameter. Must be either TRUE or FALSE.\n")
    }
    if(!analyzemRNA%in%c("TRUE","FALSE")){
        stop("analyzemRNA parameter must be either TRUE or FALSE.\n")
    }
    if(is.null(useRVM)){
        stop("Please provide useRVM parameter. Must be either TRUE or FALSE.\n")
    }
    if(!useRVM%in%c("TRUE","FALSE")){
        stop("useRVM parameter must be either TRUE or FALSE.\n")
    }
    
    
    
    
    
    anota2seqCheckInput(dataP = Anota2seqDataSet@dataP,
                        dataT = Anota2seqDataSet@dataT,
                        phenoVec = Anota2seqDataSet@phenoVec,
                        batchVec = Anota2seqDataSet@batchVec,
                        contrasts = contrasts,
                        correctionMethod=correctionMethod)
    #significance filtering parameters
    parameters <- list(minSlopeTranslation = -1,
                       maxSlopeTranslation = 2,
                       minSlopeBuffering = -2,
                       maxSlopeBuffering = 1,
                       maxPAdj = 0.15,
                       maxP = NULL,
                       minEff = NULL,
                       deltaPT = log2(1.2),
                       deltaTP = log2(1.2),
                       deltaP = NULL,
                       deltaT = NULL)
    
    ## If parameters are specified check them
    if(!is.null(thresholds)){
        for(paramNames in 1:length(thresholds)){
            if(names(thresholds)[paramNames]%in%names(parameters) == FALSE){
                message(paste("ERROR:" ,
                              names(thresholds)[paramNames],
                              " is not a recognized filtering parameter name.\nPlease read the function help for thresholds names.\n"
                              ,sep=""))
                stop()
            }
            if(names(thresholds)[paramNames]%in%names(parameters) == TRUE){
                tmpName<- names(thresholds)[paramNames]
                parameters[tmpName] <- thresholds[tmpName]
            }
        }
    }
    
    ## Perform anotaFunctions
    if(performQC==TRUE){
        if(is.null(onlyGroup)){
            stop("Please provide onlyGroup parameter. Must be either TRUE or FALSE.\n")
        }
        if(!onlyGroup%in%c("TRUE","FALSE")){
            stop("onlyGroup parameter must be either TRUE or FALSE.\n")
        }
        Anota2seqDataSet <- anota2seqPerformQC(Anota2seqDataSet,
                                               useProgBar=useProgBar,
                                               generateSingleGenePlots = generateSingleGenePlots,
                                               onlyGroup = onlyGroup)
    }
    if(performROT == TRUE){
        if(is.null(Anota2seqDataSet@qualityControl)){
            warning("The residual outlier test could not be performed because it requires anota2seqPerformQC to have successfully been run on the Anota2seqDataSet object.\n")
        } else {
            Anota2seqDataSet <- anota2seqResidOutlierTest(Anota2seqDataSet,useProgBar=useProgBar,generateSingleGenePlots = generateSingleGenePlots)
        }
    }
    
    # 4 possible options... 
    ## Perform specified analysis
    if(analyzemRNA == FALSE & analyzeBuffering == FALSE){
        analyzeVec <- "translation"
    }
    if(analyzemRNA == FALSE & analyzeBuffering == TRUE){
        analyzeVec <- c("translation","buffering")
    }
    if(analyzemRNA == TRUE & analyzeBuffering == FALSE){
        analyzeVec <- c("translated mRNA","total mRNA","translation")
    }
    if(analyzemRNA == TRUE & analyzeBuffering == TRUE){
        analyzeVec <- c("translated mRNA","total mRNA","translation","buffering")
    }
    Anota2seqDataSet <- anota2seqAnalyze(Anota2seqDataSet,
                                         contrasts= contrasts,
                                         correctionMethod = correctionMethod,
                                         useProgBar=useProgBar,
                                         analysis = analyzeVec)
    
    
    #extract the usedContrast from the translation object should be there...
    contrasts <- Anota2seqDataSet@contrasts
    message("Start filtering for significant genes ... \n")
    message("Your filtering parameters are:\n")
    message(paste("minSlopeTranslation: ", parameters$minSlopeTranslation,"\n"),sep="")
    message(paste("maxSlopeTranslation: ", parameters$maxSlopeTranslation,"\n"),sep="")
    message(paste("minSlopeBuffering: ", parameters$minSlopeBuffering,"\n"),sep="")
    message(paste("maxSlopeBuffering: ", parameters$maxSlopeBuffering,"\n"),sep="")
    message(paste("maxPAdj: ", parameters$maxPAdj,"\n"),sep="")
    message(paste("maxP: ", parameters$maxP,"\n"),sep="")
    message(paste("minEff: ", parameters$minEff,"\n"),sep="")
    message(paste("deltaPT: ", parameters$deltaPT,"\n"),sep="")
    message(paste("deltaTP: ", parameters$deltaTP,"\n"),sep="")
    message(paste("deltaP: ", parameters$deltaP,"\n"),sep="")
    message(paste("deltaT: ", parameters$deltaT,"\n"),sep="")
    message("\n")
    Anota2seqDataSet <- anota2seqSelSigGenes(Anota2seqDataSet,
                                             selContrast = c(1:dim(contrasts)[2]),
                                             useRVM=useRVM,
                                             minSlopeTranslation = parameters$minSlopeTranslation,
                                             maxSlopeTranslation = parameters$maxSlopeTranslation,
                                             minSlopeBuffering = parameters$minSlopeBuffering,
                                             maxSlopeBuffering = parameters$maxSlopeBuffering,
                                             maxPAdj = parameters$maxPAdj,
                                             maxP = parameters$maxP,
                                             minEff = parameters$minEff,
                                             selDeltaPT = parameters$deltaPT,
                                             selDeltaT = parameters$deltaT,
                                             selDeltaTP = parameters$deltaTP,
                                             selDeltaP=parameters$deltaP,
                                             analysis=analyzeVec)
    noNull <- TRUE
    # ADD test - if anota2seqSelSigGenes for is NULL for any of the analysis stop here
    for( reg in 1:length(analyzeVec)){
        for(contr in 1:dim(Anota2seqDataSet@contrasts)[2]){
            if(is.null(anota2seqGetOutput(object = Anota2seqDataSet,analysis = analyzeVec[reg],output = "selected",selContrast = contr,getRVM=useRVM))){
                noNull <- FALSE
                warning(paste("No significant genes found for analysis of ", analyzeVec[reg], " contrast ",contr,".\n No assessment of regulatory modes possible.\n"))
            }
        }
    }
    
    
    if(analyzemRNA == TRUE & analyzeBuffering == TRUE & noNull == TRUE){
        Anota2seqDataSet <- anota2seqRegModes(Anota2seqDataSet)
    }
    
    if(analyzemRNA == FALSE | analyzeBuffering == FALSE){
        warning("analyzeBuffering and/or analyzemRNA parameter set to FALSE, no assessment of regulatory modes possible.\n For assessment of regulatory modes analyzemRNA and analyzeBuffering must be set to TRUE.\n")   
    }
    
    return(Anota2seqDataSet)
}

