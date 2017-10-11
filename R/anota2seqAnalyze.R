anota2seqAnalyze <- function (Anota2seqDataSet, contrasts=NULL, correctionMethod = "BH",
                              useProgBar = TRUE, fileStem = "ANOTA2SEQ",
                              analysis = c("translation","buffering","translated mRNA","total mRNA"))
{
    if(is.null(Anota2seqDataSet)){
        stop("Please provide an Anota2seqDataSet.\n")
    }
    if(class(Anota2seqDataSet)!= "Anota2seqDataSet"){
        stop("Please provide an Anota2seqDataSet.\n")
    }
    anota2seqCheckInput(dataP = Anota2seqDataSet@dataP,
                        dataT = Anota2seqDataSet@dataT,
                        phenoVec = Anota2seqDataSet@phenoVec,
                        batchVec = Anota2seqDataSet@batchVec,
                        contrasts = contrasts,
                        correctionMethod=correctionMethod)
    if(is.null(contrasts) == FALSE){
        Anota2seqDataSet@contrasts <- contrasts
    }
    
    #########################################################
    anota2seqCheckParameter(analysis = analysis,
                            inFunc = "analysis")
    phenoVec <- Anota2seqDataSet@phenoVec
    batchVec <- Anota2seqDataSet@batchVec
    contrasts <- Anota2seqDataSet@contrasts
    
    
    for(reg in 1:length(analysis)){
        
        if(analysis[reg]=="translated mRNA"){
            message("Starting analysis of translated mRNA (RPFs). \n")
            dataP <- Anota2seqDataSet@dataP
            dataT <- NULL
        }
        if(analysis[reg]=="total mRNA"){
            message("Starting analysis of total mRNA. \n")
            dataP <- Anota2seqDataSet@dataT
            dataT <- NULL
        }
        if(analysis[reg] =="translation"){
            dataP <- Anota2seqDataSet@dataP
            dataT <- Anota2seqDataSet@dataT
            message("Starting analysis of translation.\n")
        }
        if(analysis[reg] == "buffering"){
            dataP <- Anota2seqDataSet@dataT
            dataT <- Anota2seqDataSet@dataP
            message("Starting analysis of translational buffering.\n")
        }
        
        phenoVecOrg <- phenoVec
        phenoVec <- as.factor(phenoVec)
        if (is.null(batchVec) == FALSE) {
            batchVec <- as.factor(batchVec)
        }
        nPheno <- length(levels(phenoVec))
        phenoLev <- levels(phenoVec)
        nData <- dim(dataP)[1]
        
        
        #message("Running anotaGetSigGenes\n")
        if (is.null(contrasts)) {
            message("\tUsing default \"treatment\" contrasts between (custom contrasts can be set):\n")
            message("\t", paste(levels(phenoVec),sep="  "))
            message("\n")
        }
        
        
        residRvmDf <- apvP <- apvRvmP <- apvPAdj <- apvRvmPAdj <- apvEff <- apvMSerror <- apvRvmMSerror <- apvT <- apvF <- apvFP <- apvRvmF <- apvN <- apvRvmDf <- apvDf <- matrix(nrow = nData,
                                                                                                                                                                                   ncol = nPheno - 1)
        residDf <- unadjustedResidError <- apvSlope <- apvSlopeP <- c(rep(NA,
                                                                          nData))
        rownames(apvP) <- rownames(apvRvmP) <- rownames(apvRvmDf) <- rownames(apvDf) <- rownames(apvEff) <- rownames(apvRvmF) <- rownames(apvRvmF) <- rownames(apvRvmMSerror) <- names(residDf) <- names(apvSlope) <- names(apvSlopeP) <- rownames(residRvmDf) <- rownames(dataP)
        total <- nData
        for (i in 1:nData) {
            tmpList <- list(PolyRNA = dataP[i, ], TotalRNA = dataT[i,], phenoType = phenoVec, batchVec = batchVec)
            tmpApvSum <- c()
            tmpDim <- c()
            tmpApv <- c()
            if (is.null(contrasts) == FALSE) {
                contrasts(tmpList$phenoType) <- contrasts
            }
            if (is.null(contrasts) == TRUE) {
                tmpCont <- contrasts(tmpList$phenoType)
                tmpCont[1, ] <- (-1)
                contrasts(tmpList$phenoType) <- tmpCont
                Anota2seqDataSet@contrasts <- tmpCont
            }
            if (i == 1) {
                tmpOut <- contrasts(tmpList$phenoType)
                tmpColSum <- apply(tmpOut, 2, sum)
                if (sum(tmpColSum == 0) != length(tmpColSum)) {
                    message()
                    stop("\nIncorrect contrasts. Each column of the contrast matrix should sum to 0\n\n")
                }
                colnames(tmpOut) <- paste("contrast", c(1:dim(tmpOut)[2]))
                message("\tThese contrasts will be evaluated:\n\n")
                print(tmpOut)
                message("\n")
                if (useProgBar == TRUE) {
                    pb <- txtProgressBar(min = 0, max = total, style = 3)
                }
            }
            ## attach(tmpList)
            ## rewritten to avoid use of attach; should remove the error msg
            modelAlt <- NA
            if (is.null(batchVec) == TRUE & is.null(dataT) == FALSE) {
                tmpApv <- lm(dataP[i,] ~ dataT[i,] + phenoVec)
                modelAlt <- 1
            }
            if (is.null(batchVec) == FALSE & is.null(dataT) == FALSE) {
                tmpApv <- lm(dataP[i,] ~ dataT[i,] + phenoVec + batchVec)
                modelAlt <- 1
            }
            if (is.null(batchVec) == TRUE & is.null(dataT) == TRUE) {
                tmpApv <- lm(dataP[i,] ~ phenoVec)
                modelAlt <- 2
            }
            if (is.null(batchVec) == FALSE & is.null(dataT) == TRUE) {
                tmpApv <- lm(dataP[i,] ~ phenoVec + batchVec)
                modelAlt <- 2
            }
            ##detach(tmpList)
            tmpApvSum <- summary(tmpApv)
            if (modelAlt == 1) {
                tmpError <- tmpApvSum$sigma
                unadjustedResidError[i] <- tmpError
                tmpSlope <- tmpApvSum$coefficients[2, 1]
                apvSlope[i] <- tmpSlope
                apvSlopeP[i] <- 1
                
                if(analysis[reg] =="translation"){
                    if (apvSlope[i] < 0 | apvSlope[i] > 1) {
                        apvSlopeP[i] <- anota2seqSlopeTest(tmpLm = tmpApv,
                                                           curSlope = apvSlope[i],analysis ="translation")
                    }
                }
                if(analysis[reg] == "buffering"){
                    if (apvSlope[i] < -1 | apvSlope[i] > 0) {
                        apvSlopeP[i] <- anota2seqSlopeTest(tmpLm = tmpApv,
                                                           curSlope = apvSlope[i],analysis = "buffering")
                    }
                }
                
                contMat <- contrasts(tmpList$phenoType)
                for (j in 1:dim(contMat)[2]) {
                    contMatRed <- contMat[contMat[, j] != 0, j]
                    tmpContMat <- matrix(nrow = length(contMatRed),
                                         ncol = 4)
                    rownames(tmpContMat) <- names(contMatRed)
                    colnames(tmpContMat) <- c("contrast", "estimate",
                                              "groupN", "T_mean")
                    for (k in 1:dim(tmpContMat)[1]) {
                        tmpGroup <- rownames(tmpContMat)[k]
                        tmpEst <- mean(dataP[i, phenoVecOrg == tmpGroup]) -
                            (tmpSlope * mean(dataT[i, phenoVecOrg ==
                                                       tmpGroup]))
                        tmpCovMean <- mean(dataT[i, phenoVecOrg ==
                                                     tmpGroup])
                        tmpContMat[k, "contrast"] <- contMatRed[k]
                        tmpContMat[k, "estimate"] <- tmpEst
                        tmpContMat[k, "groupN"] <- sum(phenoVecOrg ==
                                                           tmpGroup)
                        tmpContMat[k, "T_mean"] <- tmpCovMean
                    }
                    tmpContMatCopy <- tmpContMat
                    tmpContMatCopy[tmpContMatCopy[, "contrast"] >
                                       0, "contrast"] <- tmpContMatCopy[tmpContMatCopy[,
                                                                                       "contrast"] > 0, "contrast"]/sum(abs(tmpContMatCopy[tmpContMatCopy[,
                                                                                                                                                          "contrast"] > 0, "contrast"]))
                    tmpContMatCopy[tmpContMatCopy[, "contrast"] <
                                       0, "contrast"] <- tmpContMatCopy[tmpContMatCopy[,
                                                                                       "contrast"] < 0, "contrast"]/sum(abs(tmpContMatCopy[tmpContMatCopy[,
                                                                                                                                                          "contrast"] < 0, "contrast"]))
                    tmpDiffEff <- sum(tmpContMatCopy[, "contrast"] *
                                          tmpContMatCopy[, "estimate"])
                    lmT <- lm(dataT[i, ] ~ as.factor(phenoVecOrg))
                    lmTAov <- anova(lmT)
                    tmpCovSS <- lmTAov[2, 2]
                    tmpCovMean <- mean(dataT[i, ])
                    tmpCov <- sum(tmpContMatCopy[, "contrast"] *
                                      (tmpContMat[, "T_mean"] - tmpCovMean))
                    tmpN <- sum(tmpContMatCopy[, "contrast"] * tmpContMatCopy[,
                                                                              "contrast"]/tmpContMat[, "groupN"])
                    tmpErrorAdjusted <- tmpError * tmpError * c(tmpN +
                                                                    (tmpCov * tmpCov/tmpCovSS))
                    apvMSerror[i, j] <- tmpErrorAdjusted
                    apvEff[i, j] <- tmpDiffEff
                    apvDf[i, j] <- 1
                }
                residDf[i] <- tmpApvSum$df[2]
                if (useProgBar == TRUE) {
                    setTxtProgressBar(pb, i)
                }
            }
            if (modelAlt == 2) {
                tmpError <- tmpApvSum$sigma
                unadjustedResidError[i] <- tmpError
                tmpSlope <- NA
                apvSlope[i] <- tmpSlope
                apvSlopeP[i] <- NA
                contMat <- contrasts(tmpList$phenoType)
                for (j in 1:dim(contMat)[2]) {
                    contMatRed <- contMat[contMat[, j] != 0, j]
                    tmpContMat <- matrix(nrow = length(contMatRed),
                                         ncol = 4)
                    rownames(tmpContMat) <- names(contMatRed)
                    colnames(tmpContMat) <- c("contrast", "estimate",
                                              "groupN", "T_mean")
                    for (k in 1:dim(tmpContMat)[1]) {
                        tmpGroup <- rownames(tmpContMat)[k]
                        tmpEst <- mean(dataP[i, phenoVecOrg == tmpGroup])
                        tmpContMat[k, "contrast"] <- contMatRed[k]
                        tmpContMat[k, "estimate"] <- tmpEst
                        tmpContMat[k, "groupN"] <- sum(phenoVecOrg ==
                                                           tmpGroup)
                        tmpContMat[k, "T_mean"] <- NA
                    }
                    tmpContMatCopy <- tmpContMat
                    tmpContMatCopy[tmpContMatCopy[, "contrast"] >
                                       0, "contrast"] <- tmpContMatCopy[tmpContMatCopy[,
                                                                                       "contrast"] > 0, "contrast"]/sum(abs(tmpContMatCopy[tmpContMatCopy[,
                                                                                                                                                          "contrast"] > 0, "contrast"]))
                    tmpContMatCopy[tmpContMatCopy[, "contrast"] <
                                       0, "contrast"] <- tmpContMatCopy[tmpContMatCopy[,
                                                                                       "contrast"] < 0, "contrast"]/sum(abs(tmpContMatCopy[tmpContMatCopy[,
                                                                                                                                                          "contrast"] < 0, "contrast"]))
                    tmpDiffEff <- sum(tmpContMatCopy[, "contrast"] *
                                          tmpContMatCopy[, "estimate"])
                    tmpN <- sum(tmpContMatCopy[, "contrast"] * tmpContMatCopy[,
                                                                              "contrast"]/tmpContMat[, "groupN"])
                    tmpErrorAdjusted <- tmpError * tmpError * c(tmpN)
                    apvMSerror[i, j] <- tmpErrorAdjusted
                    apvEff[i, j] <- tmpDiffEff
                    apvDf[i, j] <- 1
                }
                residDf[i] <- tmpApvSum$df[2]
                if (useProgBar == TRUE) {
                    setTxtProgressBar(pb, i)
                }
            }
        }
        message("\n\n")
        
        
        apvF <- apvEff * apvEff/apvMSerror
        for (j in 1:(nPheno - 1)) {
            apvFP[, j] <- 1 - pf(apvF[, j], apvDf[, j], residDf)
            if (correctionMethod != "qvalue") {
                apvPAdj[, j] <- anota2seqAdjustPvals(pVals = apvFP[,
                                                                   j], correctionMethod = correctionMethod)
            }
            if (correctionMethod == "qvalue") {
                apvPAdj[, j] <- anota2seqAdjustPvalsQ(pVals = apvFP[,
                                                                    j])
            }
        }
        
        statsList <- list()
        for (j in 1:dim(contMat)[2]) {
            tmpMat <- cbind(apvSlope, apvSlopeP, unadjustedResidError,
                            apvEff[, j], apvMSerror[, j], apvF[, j], residDf,
                            apvFP[, j], apvPAdj[, j])
            colnames(tmpMat) <- c("apvSlope", "apvSlopeP", "unadjustedResidError",
                                  "apvEff", "apvMSerror", "apvF", "residDf", "apvP",
                                  "apvPAdj")
            statsList[[j]] <- tmpMat
        }
        contrastMat <- contrasts(tmpList$phenoType)
        contrastNames <- colnames(contMat)
        abList <- list()
        statsListRvm <- list()
        
        message("\tCalculating RVM statistics\n")
        jpeg(paste(fileStem, "_rvm_fit_for_all_contrasts_group.jpg",
                   sep = ""), width = 800, height = c((nPheno - 1) *
                                                          400), quality = 100)
        par(mfrow = c(c(nPheno - 1), 2))
        for (j in 1:c(nPheno - 1)) {
            message("\tAssessing contrast", j, "\n")
            
            anota2seqPlotIGFit(apvMSerror[, j], residDf[1], qqName = paste("Fit for contrast",
                                                                           j))
            
            tmpRVM <- anota2seqPerformRVM(MS = apvEff[, j] * apvEff[,
                                                                    j], Df = apvDf[, j], residDf = residDf, residMS = apvMSerror[,
                                                                                                                                 j])
            apvRvmMSerror[, j] <- tmpRVM$residMSRvm
            residRvmDf[, j] <- tmpRVM$residDfRvm
            apvRvmF[, j] <- tmpRVM$rvmFval
            apvRvmP[, j] <- tmpRVM$rvmP
            abList[[j]] <- tmpRVM$ab
            if (correctionMethod != "qvalue") {
                apvRvmPAdj[, j] <- anota2seqAdjustPvals(pVals = apvRvmP[,
                                                                        j], correctionMethod = correctionMethod)
            }
            
            if (correctionMethod == "qvalue") {
                apvRvmPAdj[, j] <- anota2seqAdjustPvalsQ(pVals = apvRvmP[,
                                                                         j])
            }
            
        }
        dev.off()
        for (j in 1:dim(contMat)[2]) {
            tmpMat <- cbind(apvSlope, apvSlopeP, apvEff[, j],
                            apvRvmMSerror[, j], apvRvmF[, j], residRvmDf[,
                                                                         j], apvRvmP[, j], apvRvmPAdj[, j])
            tmpColnames <- c("apvSlope", "apvSlopeP", "apvEff",
                             "apvRvmMSerror", "apvRvmF", "residRvmDf", "apvRvmP",
                             "apvRvmPAdj")
            colnames(tmpMat) <- tmpColnames
            statsListRvm[[j]] <- tmpMat
        }
        
        if(analysis[reg]=="total mRNA"){
            dataT <- dataP
            dataP <- NULL
        }
        
        if(analysis[reg] == "buffering"){
            tmpDataT <- dataT
            tmpDataP <- dataP
            
            dataP <- tmpDataT
            dataT <- tmpDataP
            
        }
        
        if(is.null(Anota2seqDataSet@deltaData)){
            Anota2seqDataSet@deltaData <- rep(list(NULL),dim(Anota2seqDataSet@contrasts)[2])
            for(cont in 1:dim(Anota2seqDataSet@contrasts)[2]){
                Anota2seqDataSet@deltaData[[cont]] <- anota2seqCalculateDeltaData(dataP=Anota2seqDataSet@dataP,
                                                                                  dataT=Anota2seqDataSet@dataT,
                                                                                  phenoVec=Anota2seqDataSet@phenoVec,
                                                                                  useIds = rownames(Anota2seqDataSet@dataP),
                                                                                  contrasts=Anota2seqDataSet@contrasts,
                                                                                  selContr=cont)
            }
        }
        outputList <- new("Anota2seqOutput", 
                          apvStats = statsList, 
                          apvStatsRvm = statsListRvm,
                          correctionMethod = correctionMethod, 
                          usedContrasts = contrastMat,
                          abList = abList)
        
        if(analysis[reg] == "translated mRNA"){
            Anota2seqDataSet@translatedmRNA <- outputList
            message("End of translated mRNA analysis.\n")
        }
        if(analysis[reg] == "total mRNA"){
            Anota2seqDataSet@totalmRNA <- outputList
            message("End of total mRNA analysis.\n")
        }
        
        if(analysis[reg]== "translation"){
            Anota2seqDataSet@translation <- outputList
            message("End of analysis of translation.\n")
        }
        
        if(analysis[reg] == "buffering"){
            Anota2seqDataSet@buffering <- outputList
            message("End of analysis of translational buffering.\n")
        }
    }
    return(Anota2seqDataSet)
}

