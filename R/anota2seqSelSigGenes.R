
anota2seqSelSigGenes <- function (Anota2seqDataSet, useRVM = TRUE,
                                  analysis = anota2seqGetAvailableAnalyzes(Anota2seqDataSet),
                                  selIds = NULL,
                                  selContrast = seq(along = 1:dim(anota2seqGetContrasts(Anota2seqDataSet))[2]),
                                  minSlopeTranslation = NULL,
                                  maxSlopeTranslation = NULL, minSlopeBuffering = NULL,
                                  maxSlopeBuffering = NULL,  slopeP = NULL, minEff = NULL, maxP = NULL,
                                  maxPAdj = NULL, selDeltaPT = NULL, selDeltaTP = NULL,
                                  selDeltaP = NULL, selDeltaT = NULL, sortBy = c("rvmP", "none", "Eff", "p"))
{
    
    if(is.null(Anota2seqDataSet)){
        stop("Please provide an Anota2seqDataSet.\n")
    }
    if(class(Anota2seqDataSet)!= "Anota2seqDataSet"){
        stop("Please provide an Anota2seqDataSet.\n")
    }
    
    if(is.null(anota2seqGetOutputClass(Anota2seqDataSet,"translated mRNA","full"))&
       is.null(anota2seqGetOutputClass(Anota2seqDataSet,"total mRNA","full"))&
       is.null(anota2seqGetOutputClass(Anota2seqDataSet,"translation","full"))&
       is.null(anota2seqGetOutputClass(Anota2seqDataSet,"buffering","full"))){
        stop("No anota2seqAnalyze output found in the Anota2seqDataSet. Please run anota2seqAnalyze before using anota2seqSelSigGenes.\n")
    }
    
    if(is.null(analysis) == TRUE){
        stop("Please provide the analysis parameter.\nMust be one of the following: translated mRNA, total mRNA, translation, buffering.")
    }
    
    if(is.null(selContrast)){
        stop("Please provide one or more contrasts using the selContrast parameter.\n")
    }
    
    if(max(selContrast) > dim(Anota2seqDataSet@contrasts)[2]){
        stop("One of the selected contrasts does not exist. selContrast must be a numeric vector with 1 or more contrasts.\n The values cannot be greater than the number of columns in the contrast matrix.")
    }
    if(is.null(useRVM)){
        stop("Please provide useRVM parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!useRVM%in%c(TRUE,FALSE)){
        stop("useRVM parameter must be set to TRUE or FALSE.\n")
    }
    if(length(sortBy) <2){
        if(!sortBy %in% c("p","rvmP","none","Eff")){
            stop("sortBy parameter must be one of the following; p, rvmP, Eff, none.\n")
        }
    }
    if(length(sortBy) == 4){
        if(useRVM == FALSE){
            sortBy <- "p"
        }
        if(useRVM == TRUE){
            sortBy <- "rvmP"
        }
    }
    
    if(length(analysis) < 2){
        if(analysis == "translation"){
            if(is.null(selDeltaT) == FALSE){
                stop("selDeltaT is set and alysis is set to translation only. Please check your parameter settings.\n") 
            }
            if(is.null(selDeltaTP) == FALSE){
                stop("selDeltaTP is set and alysis is set to translation only. Please check your parameter settings.\n") 
            }
            if(is.null(minSlopeBuffering) == FALSE){
                stop("minSlopeBuffering is set and alysis is set to translation only. Please check your parameter settings.\n") 
            }
            if(is.null(maxSlopeBuffering) == FALSE){
                stop("maxSlopeBuffering is set and alysis is set to translation only. Please check your parameter settings.\n")  
            }
        }
        if(analysis == "buffering"){
            if(is.null(selDeltaP) == FALSE){
                stop("selDeltaP is set and alysis is set to buffering only. Please check your parameter settings.\n") 
            }
            if(is.null(selDeltaPT) == FALSE){
                stop("selDeltaPT is set and alysis is set to buffering only. Please check your parameter settings.\n") 
            }
            if(is.null(minSlopeTranslation) == FALSE){
                stop("minSlopeTranslation is set and alysis is set to buffering only. Please check your parameter settings.\n") 
            }
            if(is.null(maxSlopeTranslation) == FALSE){
                stop("maxSlopeTranslation is set and alysis is set to buffering only. Please check your parameter settings.\n")  
            } 
        }
        if(analysis =="translated mRNA"){
            if(is.null(selDeltaT) == FALSE){
                stop("selDeltaT is set and alysis is set to translated mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(selDeltaPT) == FALSE){
                stop("selDeltaPT is set and alysis is set to translated mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(selDeltaTP) == FALSE){
                stop("selDeltaTP is set and alysis is set to translated mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(slopeP) == FALSE){
                stop("slopeP is set and alysis is set to translated mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(minSlopeTranslation) == FALSE){
                stop("minSlopeTranslation is set and alysis is set to translated mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(maxSlopeTranslation) == FALSE){
                stop("maxSlopeTranslation is set and alysis is set to translated mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(minSlopeBuffering) == FALSE){
                stop("minSlopeBuffering is set and alysis is set to translated mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(maxSlopeBuffering) == FALSE){
                stop("maxSlopeBuffering is set and alysis is set to translated mRNA only. Please check your parameter settings.\n") 
            }
        }
        if(analysis =="total mRNA"){
            if(is.null(selDeltaP) == FALSE){
                stop("selDeltaP is set and alysis is set to total mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(selDeltaPT) == FALSE){
                stop("selDeltaPT is set and alysis is set to total mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(selDeltaTP) == FALSE){
                stop("selDeltaTP is set and alysis is set to total mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(slopeP) == FALSE){
                stop("slopeP is set and alysis is set to total mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(minSlopeTranslation) == FALSE){
                stop("minSlopeTranslation is set and alysis is set to total mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(maxSlopeTranslation) == FALSE){
                stop("maxSlopeTranslation is set and alysis is set to total mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(minSlopeBuffering) == FALSE){
                stop("minSlopeBuffering is set and alysis is set to total mRNA only. Please check your parameter settings.\n") 
            }
            if(is.null(maxSlopeBuffering) == FALSE){
                stop("maxSlopeBuffering is set and alysis is set to total mRNA only. Please check your parameter settings.\n") 
            }
        }
    }
    
    # check cutoffs 
    cutOffVec <- c(minSlopeTranslation,maxSlopeTranslation,
                   minSlopeBuffering,maxSlopeBuffering,
                   slopeP,minEff,maxP,maxPAdj,selDeltaPT,
                   selDeltaTP,selDeltaP,selDeltaT)
    for(cutOff in 1:length(cutOffVec)){
        if(class(cutOffVec[cutOff])%in%c("NULL","numeric") == FALSE){
            stop("Supplied filter parameters must be numeric.\n")
        }
    }
    
    anota2seqCheckParameter(analysis=analysis,inFunc="analysis")
    
    if(useRVM == TRUE){
        maxRvmPAdj <- maxPAdj
        maxRvmP <- maxP
        maxPAdj <- NULL
        maxP <- NULL
    }
    if(useRVM == FALSE){
        maxRvmP <- NULL
        maxRvmPAdj <- NULL
    }
    
    for(reg in 1:length(analysis)){
        for(cont in 1:length(selContrast)){
            if(is.null(anota2seqGetOutput(object = Anota2seqDataSet,
                                            analysis = analysis[reg],
                                            selContrast = selContrast[cont],
                                            output = "full",
                                            getRVM = useRVM))){
                stop(paste("No anota2seqAnalyze parameter for ",analysis," analysis found.\nPlease change analysis parameter or run anota2seqAnalyze using the specified analysis parameter.\n"))
            }
        }
    }
    # Need to reset the delta filtering option ...
    deltaVec <- list(selDeltaP = selDeltaP,selDeltaT=selDeltaT,selDeltaPT = selDeltaPT,selDeltaTP =selDeltaTP)
    for(reg in 1:length(analysis)){
        message(paste("Starting filtering of: ",analysis[reg],"\n",sep=""))
        if( analysis[reg] =="translated mRNA"){
            maxSlope <- NULL
            minSlope <- NULL
            slopeP <- NULL
            selDeltaPT <- NULL
            selDeltaTP <- NULL
            selDeltaT <- NULL
            selDeltaP <- deltaVec$selDeltaP
        }
        if( analysis[reg] =="total mRNA"){
            maxSlope <- NULL
            minSlope <- NULL
            slopeP <- NULL
            selDeltaPT <- NULL
            selDeltaTP <- NULL
            selDeltaP <- NULL
            selDeltaT <- deltaVec$selDeltaT
        }
        if( analysis[reg] == "translation"){
            minSlope <- minSlopeTranslation
            maxSlope <- maxSlopeTranslation
            selDeltaTP <- NULL
            selDeltaT <- NULL
            selDeltaP <- deltaVec$selDeltaP
            selDeltaPT <- deltaVec$selDeltaPT
        }
        if(analysis[reg] == "buffering"){
            minSlope <- minSlopeBuffering
            maxSlope <- maxSlopeBuffering
            selDeltaPT <- NULL
            selDeltaP <- NULL
            selDeltaTP <- deltaVec$selDeltaTP
            selDeltaT <- deltaVec$selDeltaT
        }
        
        
        for(contr in 1:length(selContrast)){
            message(paste("\tfiltering contrast ",selContrast[contr],".\n",sep=""))
            ### Set anota2seqSigObj according to analyse
            anota2seqSigObj <- anota2seqGetOutputClass(Anota2seqDataSet,analysis = analysis[reg],output = "full")
            if(is.null(anota2seqGetOutputClass(Anota2seqDataSet,analysis[reg],"selected")) == TRUE){
                Anota2seqDataSet <- anota2seqSetOutput(Anota2seqDataSet,
                                               analysis[reg],
                                               "selected",
                                               new("Anota2seqSelectedOutput",
                                                   selectedData = rep(list(NULL),dim(anota2seqSigObj@usedContrasts)[2]),
                                                   selectedRvmData = rep(list(NULL),dim(anota2seqSigObj@usedContrasts)[2]),
                                                   useRVM = useRVM,
                                                   deltaData = rep(list(NULL),dim(anota2seqSigObj@usedContrasts)[2]),
                                                   usedThresholds = rep(list(NULL),dim(anota2seqSigObj@usedContrasts)[2]),
                                                   regModes = FALSE
                                               )
                )
            }
            dataT <- Anota2seqDataSet@dataT
            dataP <- Anota2seqDataSet@dataP
            phenoVec <- Anota2seqDataSet@phenoVec
            
            if (selContrast[contr] > dim(anota2seqSigObj@usedContrasts)[2]) {
                
                stop("Specified contrast does not exist")
            }
            tmpData <- anota2seqSigObj@apvStats[[selContrast[contr]]]
            tmpDataRvm <- NULL
            if (useRVM == TRUE) {
                tmpDataRvm <- anota2seqSigObj@apvStatsRvm[[selContrast[contr]]]
            }
            
            deltaP <- as.matrix(Anota2seqDataSet@deltaData[[selContrast[contr]]][,"deltaP"])
            deltaT <- as.matrix(Anota2seqDataSet@deltaData[[selContrast[contr]]][,"deltaT"])
            deltaPT <- as.matrix(Anota2seqDataSet@deltaData[[selContrast[contr]]][,"deltaPT"])
            deltaTP <- as.matrix(Anota2seqDataSet@deltaData[[selContrast[contr]]][,"deltaTP"])
            colnames(deltaP) <- "deltaP"
            colnames(deltaT) <- "deltaT"
            colnames(deltaPT) <- "deltaPT"
            colnames(deltaTP) <- "deltaTP"
            rownames(deltaP) <- rownames(dataP)
            rownames(deltaT) <- rownames(dataP)
            rownames(deltaPT) <- rownames(dataP)
            rownames(deltaTP) <- rownames(dataP)
            if (is.null(selIds) == FALSE) {
                useIds <- selIds
                tmpData <- tmpData[useIds, ,drop=FALSE]
                if (useRVM == TRUE) {
                    tmpDataRvm <- tmpDataRvm[useIds, ,drop=FALSE]
                }
            }
            if (is.null(selIds) == TRUE) {
                if (is.null(minSlope) == FALSE) {
                    tmpData <- tmpData[tmpData[, "apvSlope"] > minSlope,
                                       ,drop=FALSE]
                    if (useRVM == TRUE) {
                        tmpDataRvm <- tmpDataRvm[tmpDataRvm[, "apvSlope"] >
                                                     minSlope, ,drop=FALSE]
                    }
                }
                if (is.null(maxSlope) == FALSE) {
                    tmpData <- tmpData[tmpData[, "apvSlope"] < maxSlope,
                                       ,drop=FALSE]
                    if (useRVM == TRUE) {
                        tmpDataRvm <- tmpDataRvm[tmpDataRvm[, "apvSlope"] <
                                                     maxSlope, ,drop=FALSE]
                    }
                }
                if (is.null(slopeP) == FALSE) {
                    tmpData <- tmpData[tmpData[, "apvSlopeP"] > slopeP,
                                       ,drop=FALSE]
                    if (useRVM == TRUE) {
                        tmpDataRvm <- tmpDataRvm[tmpDataRvm[, "apvSlopeP"] >
                                                     slopeP, ,drop=FALSE]
                    }
                }
                if (is.null(minEff) == FALSE) {
                    tmpData <- tmpData[abs(tmpData[, "apvEff"]) > minEff,
                                       ,drop=FALSE]
                    if (useRVM == TRUE) {
                        tmpDataRvm <- tmpDataRvm[abs(tmpDataRvm[, "apvEff"]) >
                                                     minEff, ,drop=FALSE]
                    }
                }
                if (is.null(maxP) == FALSE) {
                    tmpNames <- rownames(tmpData[tmpData[, "apvP"] <
                                                     maxP, ,drop=FALSE])
                    tmpData <- tmpData[tmpNames, ,drop=FALSE]
                    if (useRVM == TRUE) {
                        tmpDataRvm <- tmpDataRvm[tmpNames, ,drop=FALSE]
                    }
                }
                if (is.null(maxPAdj) == FALSE) {
                    tmpNames <- rownames(tmpData[tmpData[, "apvPAdj"] <
                                                     maxPAdj, ,drop=FALSE])
                    tmpData <- tmpData[tmpNames, ,drop=FALSE]
                    if (useRVM == TRUE) {
                        tmpDataRvm <- tmpDataRvm[tmpNames, ,drop=FALSE]
                    }
                }
                if (is.null(maxRvmP) == FALSE & useRVM == TRUE) {
                    tmpNames <- rownames(tmpDataRvm[tmpDataRvm[, "apvRvmP"] <
                                                        maxRvmP, ,drop=FALSE])
                    tmpData <- tmpData[tmpNames, ,drop=FALSE]
                    tmpDataRvm <- tmpDataRvm[tmpNames,,drop=FALSE]
                }
                if (is.null(maxRvmPAdj) == FALSE & useRVM == TRUE) {
                    tmpNames <- rownames(tmpDataRvm[tmpDataRvm[, "apvRvmPAdj"] <
                                                        maxRvmPAdj, ,drop=FALSE])
                    tmpData <- tmpData[tmpNames, ,drop=FALSE]
                    tmpDataRvm <- tmpDataRvm[tmpNames, ,drop=FALSE]
                }
                useIds <- rownames(tmpData)
            }
            if (is.null(selDeltaP) == FALSE & analysis[reg] %in% c("translation", "translated mRNA")) {
                tmpDelta <- (deltaP[useIds,] > selDeltaP & anota2seqSigObj@apvStats[[selContrast[contr]]][useIds,
                                                                                                          "apvEff", drop = FALSE] > 0) | (deltaP[useIds,] < (-selDeltaP) &
                                                                                                                                              anota2seqSigObj@apvStats[[selContrast[contr]]][useIds, "apvEff",
                                                                                                                                                                                             drop = FALSE] < 0)
                
                useIds <- useIds[tmpDelta == TRUE]
                
                tmpData <- tmpData[useIds, , drop = FALSE]
                if (useRVM == TRUE) {
                    tmpDataRvm <- tmpDataRvm[useIds, , drop = FALSE]
                }
            }
            
            
            if (is.null(selDeltaT) == FALSE & analysis[reg] %in% c("buffering", "total mRNA")) {
                tmpDelta <- (deltaT[useIds,] < -(selDeltaT) & anota2seqSigObj@apvStats[[selContrast[contr]]][useIds,
                                                                                                             "apvEff", drop = FALSE] < 0) | (deltaT[useIds,] > (selDeltaT) &
                                                                                                                                                 anota2seqSigObj@apvStats[[selContrast[contr]]][useIds, "apvEff",
                                                                                                                                                                                                drop = FALSE] > 0)
                useIds <- useIds[tmpDelta == TRUE]
                tmpData <- tmpData[useIds, , drop = FALSE]
                if (useRVM == TRUE) {
                    tmpDataRvm <- tmpDataRvm[useIds, , drop = FALSE]
                }
            }
            
            if (is.null(selDeltaPT) == FALSE & analysis[reg] == "translation") {
                
                tmpDelta <- (deltaPT[useIds,] > selDeltaPT & anota2seqSigObj@apvStats[[selContrast[contr]]][useIds,
                                                                                                            "apvEff",drop=FALSE] > 0) | (deltaPT[useIds,] < (-selDeltaPT) & anota2seqSigObj@apvStats[[selContrast[contr]]][useIds,
                                                                                                                                                                                                                           "apvEff",drop=FALSE] < 0)
                useIds <- useIds[tmpDelta == TRUE]
                tmpData <- tmpData[useIds, ,drop=FALSE]
                if (useRVM == TRUE) {
                    tmpDataRvm <- tmpDataRvm[useIds, ,drop=FALSE]
                }
                
            }
            
            
            if(is.null(selDeltaTP) == FALSE & analysis[reg] == "buffering"){
                tmpDelta <- (deltaTP[useIds,] < -(selDeltaTP) & anota2seqSigObj@apvStats[[selContrast[contr]]][useIds,
                                                                                                               "apvEff", drop = FALSE] < 0) | (deltaTP[useIds,] > (selDeltaTP) &
                                                                                                                                                   anota2seqSigObj@apvStats[[selContrast[contr]]][useIds, "apvEff",
                                                                                                                                                                                                  drop = FALSE] > 0)
                useIds <- useIds[tmpDelta == TRUE]
                tmpData <- tmpData[useIds, , drop = FALSE]
                if (useRVM == TRUE) {
                    tmpDataRvm <- tmpDataRvm[useIds, , drop = FALSE]
                }
            }
            
            if (dim(tmpData)[1] < 1) {
                warning(paste("No genes pass selected thresholds for analysis of ", analysis[reg],".\n",sep=""))
            }
            
            if (is.null(sortBy) == FALSE & !sortBy %in% "none") {
                if (sortBy == "Eff") {
                    tmpSorter <- tmpData[useIds, ,drop = FALSE][order(tmpData[useIds, "apvEff", drop = FALSE]), , drop = FALSE]
                    useIds <- rownames(tmpSorter)
                }
                if (sortBy == "p") {
                    tmpSorter <- tmpData[useIds, , drop = FALSE][order(tmpData[useIds, "apvP", drop = FALSE]), , drop = FALSE]
                    useIds <- rownames(tmpSorter)
                }
                if (sortBy == "p" & (useRVM == TRUE)) {
                    stop("You have selected non-RVM based sorting but useRVM parameter is set to TRUE\n To filter based on p, set the useRVM paramter to FALSE\n")
                }
                if (sortBy == "rvmP" & (useRVM == TRUE)) {
                    tmpSorter <- tmpDataRvm[useIds, ,drop = FALSE][order(tmpDataRvm[useIds, "apvRvmP",drop = FALSE]), , drop = FALSE]
                    useIds <- rownames(tmpSorter)
                }
                if (sortBy == "rvmP" & (useRVM == FALSE)) {
                    stop("You have selected RVM based sorting but useRVM parameter is set to FALSE.\n To filter based on rvmP, set the useRVM paramter to TRUE.\n")
                }
            }
            
            tmpDataOut <- tmpData[useIds, ,drop=FALSE]
            tmpDataRvmOut <- NULL
            if (useRVM == TRUE) {
                tmpDataRvmOut <- tmpDataRvm[useIds, ,drop=FALSE]
            }
            
            if(analysis[reg] == "translated mRNA"){
                deltaT <- NULL   
                deltaTP <- NULL
                deltaPT <- NULL
            }
            if(analysis[reg] == "total mRNA"){
                deltaP <- NULL
                deltaTP <- NULL
                deltaPT <- NULL
            }
            if(analysis[reg] == "translation"){
                deltaTP <- NULL
                deltaT <- NULL
            }
            if(analysis[reg] == "buffering"){
                deltaPT <- NULL
                deltaP <- NULL
            }
            Anota2seqDataSet <- anota2seqSetSelectedOutput(Anota2seqDataSet,
                                                    analysis[reg],
                                                    selContrast[contr],
                                                    list(selectedData = as.data.frame(tmpDataOut),
                                                         selectedRvmData = as.data.frame(tmpDataRvmOut),
                                                         useRVM = useRVM,
                                                         deltaData = cbind(deltaP = deltaP[useIds, ,drop=FALSE],
                                                                           deltaT = deltaT[useIds, ,drop=FALSE],
                                                                           deltaPT = deltaPT[useIds, ,drop=FALSE],
                                                                           deltaTP = deltaTP[useIds, , drop=FALSE]),
                                                         usedThresholds = list(selContrast = selContrast[contr], minSlope = minSlope,
                                                                               maxSlope = maxSlope, slopeP = slopeP, minEff = minEff,
                                                                               maxP = maxP, maxPAdj = maxPAdj, maxRvmP = maxRvmP,
                                                                               maxRvmPAdj = maxRvmPAdj, selDeltaPT = selDeltaPT,selDeltaTP = selDeltaTP,
                                                                               selDeltaP = selDeltaP, selDeltaT = selDeltaT),
                                                         regModes = FALSE))
            message(paste("\tContrast ",selContrast[contr]," done.\n",sep=""))
        }
        message("Filtering for analysis ",analysis[reg]," done.\n")
    }
    return(Anota2seqDataSet)
}