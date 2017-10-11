anota2seqPerformQC <- function(Anota2seqDataSet, 
                               generateSingleGenePlots=FALSE, fileName="ANOTA2SEQ_translation_vs_mRNA_individual_regressions.pdf",
                               nReg=200, correctionMethod="BH", useDfb=TRUE, useDfbSim=TRUE, nDfbSimData=2000, 
                               useRVM=TRUE, onlyGroup=FALSE, useProgBar=TRUE, fileStem="ANOTA2SEQ"){
    
    if(is.null(Anota2seqDataSet)){
        stop("Please provide an Anota2seqDataSet.\n")
    }
    if(class(Anota2seqDataSet)!= "Anota2seqDataSet"){
        stop("Please provide an Anota2seqDataSet.\n")
    }
    if(is.null(generateSingleGenePlots)){
        stop("Please provide generateSingleGenePlots parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!generateSingleGenePlots%in%c(TRUE,FALSE)){
        stop("generateSingleGenePlots parameter must be set to TRUE or FALSE.\n")
    }
    if(is.null(useDfb)){
        stop("Please provide useDfb parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!useDfb%in%c(TRUE,FALSE)){
        stop("useDfb parameter must be set to TRUE or FALSE.\n")
    }
    if(is.null(useDfbSim)){
        stop("Please provide useDfbSim parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!useDfbSim%in%c(TRUE,FALSE)){
        stop("useDfbSim parameter must be set to TRUE or FALSE.\n")
    }
    if(is.null(useRVM)){
        stop("Please provide useRVM parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!useDfbSim%in%c(TRUE,FALSE)){
        stop("useRVM parameter must be set to TRUE or FALSE.\n")
    }
    if(is.null(onlyGroup)){
        stop("Please provide onlyGroup parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!onlyGroup%in%c(TRUE,FALSE)){
        stop("onlyGroup parameter must be set to TRUE or FALSE.\n")
    }
    if(is.null(useProgBar)){
        stop("Please provide useProgBar parameter. Must be set to TRUE or FALSE.\n")
    }
    if(!useProgBar%in%c(TRUE,FALSE)){
        stop("useProgBar parameter must be set to TRUE or FALSE.\n")
    }
    if(is.null(nReg)){
        stop("Please provide nReg parameter. Must be a numeric value.\n")
    }
    if(!class(nReg)%in%c("numeric")){
        stop("nReg parameter must be a numeric value.\n")
    }
    if(is.null(nDfbSimData)){
        stop("Please provide nDfbSimData parameter. Must be a numeric value of 0 or greater.\n")
    }
    if(!class(nDfbSimData)%in%c("numeric")){
        stop("nDfbSimData parameter must be a numeric value of 0 or greater.\n")
    }
    if(is.null(useRVM) | !useRVM%in%c(TRUE,FALSE)){
        stop("useRVM parameter must be set to either TRUE or FALSE.\n")
    }
    
    
    
    anota2seqCheckInput(dataP = Anota2seqDataSet@dataP,
                        dataT = Anota2seqDataSet@dataT,
                        phenoVec = Anota2seqDataSet@phenoVec,
                        batchVec = Anota2seqDataSet@batchVec,
                        contrasts = NULL,
                        correctionMethod=correctionMethod)
    
    dataP <- Anota2seqDataSet@dataP
    dataT <- Anota2seqDataSet@dataT
    phenoVec <- Anota2seqDataSet@phenoVec
    nData <- dim(dataP)[1]
    phenoVecOrg <- phenoVec
    phenoVec <- as.factor(phenoVec)
    phenoLev <- levels(phenoVec)
    nPheno <- length(phenoLev)
    ##Warnings
    ##is there sufficient replication?
    anota2seqPerformQcWarnings(nPheno=nPheno, phenoLev=phenoLev, phenoVecOrg=phenoVecOrg, onlyGroup=onlyGroup)
    ##############################
    ##initiation of objects for qc (dfbetas, interactions, slopes) and analysis (omnibus group, intercepts and rvm)
    ##dfbetas structures
    lmFittedValsAdd <- lmResidAdd <- lmDfbAdd <- matrix(ncol=dim(dataT)[2], nrow=nData)
    colnames(lmFittedValsAdd) <- colnames(lmDfbAdd) <- colnames(lmResidAdd) <- colnames(dataP)
    rownames(lmFittedValsAdd) <- rownames(lmDfbAdd) <- rownames(lmResidAdd) <- rownames(dataP)
    
    ##interactions structures
    intP <- intPAdj<- intMS <- intRvmFval <- intDf <- intResidMS <- intResidMSRvm <- intRvmP <- intRvmPAdj <- intResidDf <- intResidDfRvm <- c(rep(NA, nData))
    names(intP) <- names(intPAdj) <- names(intMS) <- names(intRvmFval) <- names(intDf) <- names(intResidMS) <- names(intResidMSRvm) <- names(intRvmP) <- names(intRvmPAdj) <- names(intResidDf) <- names(intResidDfRvm) <- rownames(dataP)
    
    ##slope structures
    groupSlope <- groupSlopeP <- c(rep(NA, nData))
    names(groupSlope) <- names(groupSlopeP) <- rownames(dataP)
    
    ##group structures
    groupP <- groupPAdj <- groupMS <- groupRvmFval <-  groupDf <-groupResidMS <-  groupResidMSRvm  <- groupRvmP <-groupRvmPAdj<- groupResidDf <- groupResidDfRvm <- c(rep(NA, nData))
    names(groupP) <- names(groupPAdj) <- names(groupMS) <- names(groupRvmFval) <-  names(groupDf) <- names(groupResidMS) <-  names(groupResidMSRvm)  <- names(groupRvmP) <- names(groupRvmPAdj) <- names(groupResidDf) <- names(groupResidDfRvm) <- rownames(dataP)
    
    ##intercept structure
    groupIntercepts <- matrix(nrow=nData, ncol=nPheno)
    colnames(groupIntercepts) <- phenoLev
    rownames(groupIntercepts) <- rownames(dataP)
    
    ##rvm structures
    names <- rownames(dataP)
    abInt <- rvmSummary <- rvmSummaryGroup <- abGroup <- dsfSummaryAdd <- NULL
    ##############################
    ##initiate optional regression plot
    if(generateSingleGenePlots==1){
        geneNames <- rownames(dataP)
        pdf(fileName, width=8, height=11, pointsize=1/600)
        par(mfrow=c(4,2))
    }
    message("Running anota2seqPerformQc quality control\n")
    ##############################
    ##Start analysis in a per gene loop
    message("\tCalculating omnibus interactions & effects and dfbetas\n")
    total <- nData
    if(useProgBar==TRUE){
        pb <- txtProgressBar(min=0, max=total, style=3)
    }
    for(i in 1:nData){
        if(useProgBar==TRUE){
            setTxtProgressBar(pb, i)
        }
        tmpList <- list("PolyRNA"=dataP[i,], "TotalRNA"=dataT[i,], "phenoType"=phenoVec)
        #attach(tmpList)
        if(onlyGroup==FALSE){
            ##perform regression with interactions
            tmpLm <- lm(dataP[i,]~dataT[i,]*phenoVec)
            ##Get omnibus interactions and stats
            tmpLmAov <- anova(tmpLm)
            intDf[i] <- tmpLmAov[3,1]
            intMS[i] <- tmpLmAov[3,3]
            intP[i] <- tmpLmAov[3,5]
            intResidDf[i] <- tmpLmAov[4,1]
            intResidMS[i] <- tmpLmAov[4,3]
        }
        ##get omnibus group effects without interaction
        tmpLm <- lm(dataP[i,]~dataT[i,]+phenoVec)
        tmpLmAov <- anova(tmpLm)
        #detach(tmpList)
        groupDf[i] <- tmpLmAov[2,1]
        groupMS[i] <- tmpLmAov[2,3]
        groupP[i] <- tmpLmAov[2,5]
        groupResidDf[i] <- tmpLmAov[3,1]
        groupResidMS[i] <- tmpLmAov[3,3]
        
        ##if slope is <0 or >1 test if there is significant difference else set to 1.
        groupSlope[i] <- tmpLm$coefficients[2]
        groupSlopeP[i] <- 1
        if(groupSlope[i]>1 | groupSlope[i]<0){
            groupSlopeP[i] <- anota2seqSlopeTest(tmpLm=tmpLm, curSlope=groupSlope[i],"translation")
        }
        
        ##get dfbetas for the slope i.e. in column 2 from the no interaction model
        if(useDfb==TRUE){
            tmpDfb <- dfbetas(tmpLm)
            lmDfbAdd[i,] <- tmpDfb[,2]
        }
        
        ##Collect residuals and fittedvalues
        lmResidAdd[i,] <- tmpLm$residuals
        lmFittedValsAdd[i,] <- tmpLm$fitted.values
        
        ##get class intercepts
        groupIntercepts[i,] <- anota2seqGetIntercepts(x=dataT[i,], y=dataP[i,], slope=groupSlope[i], phenoVecOrg=phenoVecOrg, phenoLev=phenoLev)
        
        ##Plot single gene regression
        if(generateSingleGenePlots==1 & i<=nReg){
            anota2seqPlotSingleRegression(x=dataT[i,], y=dataP[i,], geneName=geneNames[i], intercepts=groupIntercepts[i,], slope=groupSlope[i], phenoVecOrg=phenoVecOrg, phenoLev=phenoLev)
        }
    }
    message("\n\n")
    ##End plotting
    if(generateSingleGenePlots==1){
        dev.off()
    }
    ##done per gene analysis
    #############################################
    ##Dfb analysis with or without simulation
    if(useDfb==TRUE){
        message("\tAssessing dfbetas for model without interaction\n")
        dsfSummaryAdd <- anota2seqDfbsSummaryFull(lmDfb=lmDfbAdd, mode="add", filename=paste(fileStem, "_simulated_vs_obt_dfbetas_without_interaction.pdf", sep=""), useDfbSim=useDfbSim, nDfbSimData, phenoVec=phenoVecOrg, useProgBar=useProgBar)
    }
    
    #############################################
    ##RVM analysis
    ##Interactions
    if(useRVM==TRUE & onlyGroup==FALSE){
        message("\tUsing RVM for omnibus interaction statistics\n")
        jpeg(paste(fileStem, "_rvm_fit_for_interactions.jpg", sep=""), width=800, height=400, quality=100)
        par(mfrow=c(1,2))
        anota2seqPlotIGFit(intResidMS, intResidDf[1], qqName="Fit for interactions")
        dev.off()
        ##adjust the intResidMSs based on the ab parameters
        tmpRVM <- anota2seqPerformRVM(MS=intMS, Df=intDf, residMS=intResidMS, residDf=intResidDf)
        intResidMSRvm <- tmpRVM$residMSRvm
        intResidDfRvm <- tmpRVM$residDfRvm
        intRvmFval <- tmpRVM$rvmFval
        intRvmP <- tmpRVM$rvmP
        abInt <- tmpRVM$ab
    }
    rvmSummary <- cbind(intMS, intDf, intResidMS, intResidDf, intResidMSRvm, intResidDfRvm,  intRvmFval, intP, intRvmP)
    
    ##Group effect
    if(useRVM==TRUE){
        message("\tUsing RVM for omnibus group statistics\n")
        jpeg(paste(fileStem, "_rvm_fit_for_omnibus_group.jpg", sep=""), width=800, height=400, quality=100)
        par(mfrow=c(1,2))
        anota2seqPlotIGFit(groupResidMS, groupResidDf[1], qqName="Fit for omnibus group")
        dev.off()
        ##adjust the groupResidMSs based on the ab parameters
        tmpRVM <- anota2seqPerformRVM(MS=groupMS, Df=groupDf, residMS=groupResidMS, residDf=groupResidDf)
        groupResidMSRvm <- tmpRVM$residMSRvm
        groupResidDfRvm <- tmpRVM$residDfRvm
        groupRvmFval <- tmpRVM$rvmFval
        groupRvmP <- tmpRVM$rvmP
        abGroup <- tmpRVM$ab
    }
    rvmSummaryGroup <- cbind(groupSlope,groupSlopeP, groupMS, groupDf, groupResidMS, groupResidDf, groupResidMSRvm, groupResidDfRvm,  groupRvmFval, groupP, groupRvmP)
    
    ########################################
    ##Multiple testing adjustments
    message("\tAdjusting p-values for multiple testing\n\n")    
    ##one set for when multtest is used
    if(correctionMethod!="qvalue"){
        if(onlyGroup==FALSE){
            intPAdj <- anota2seqAdjustPvals(pVals=intP, correctionMethod=correctionMethod)
        }
        groupPAdj <- anota2seqAdjustPvals(pVals=groupP, correctionMethod=correctionMethod)
        if(useRVM==TRUE){
            if(onlyGroup==FALSE){
                intRvmPAdj <- anota2seqAdjustPvals(pVals=intRvmP, correctionMethod=correctionMethod)
            }
            groupRvmPAdj <- anota2seqAdjustPvals(pVals=groupRvmP, correctionMethod=correctionMethod)
        }
        rvmSummary <- cbind(rvmSummary, intPAdj, intRvmPAdj)
        rvmSummaryGroup <- cbind(rvmSummaryGroup, groupPAdj, groupRvmPAdj) 
    }
    
    ##another set for when it is storey qvalue
    if(correctionMethod=="qvalue"){
        if(onlyGroup==FALSE){
            intPAdj <- anota2seqAdjustPvalsQ(intP)
        }
        groupPAdj <- anota2seqAdjustPvalsQ(groupP)
        if(useRVM==TRUE){
            if(onlyGroup==FALSE){
                intRvmPAdj <- anota2seqAdjustPvalsQ(intRvmP)
            }
            groupRvmPAdj <- anota2seqAdjustPvalsQ(groupRvmP)
        }
        rvmSummary <- cbind(rvmSummary, intPAdj, intRvmPAdj)
        rvmSummaryGroup <- cbind(rvmSummaryGroup, groupPAdj, groupRvmPAdj) 
    }
    
    ###############################
    ##Plot for interaction p-values
    if(onlyGroup==FALSE){
        anota2seqPlotIntPvals(intP=intP, intPAdj=intPAdj, intRvmP=intRvmP, intRvmPAdj=intRvmPAdj, useRVM=useRVM, correctionMethod=correctionMethod, fileStem=fileStem)
    }
    ################################
    ##Create a return object
    dataOut <- new("Anota2seqQualityControl",
                   omniIntStats = rvmSummary,
                   omniGroupStats = rvmSummaryGroup,
                   groupIntercepts = groupIntercepts,
                   correctionMethod = correctionMethod,
                   dsfSummary = dsfSummaryAdd,
                   dfbetas = lmDfbAdd,
                   residuals = lmResidAdd,
                   fittedValues = lmFittedValsAdd,
                   phenoClasses = levels(phenoVec),
                   sampleNames = colnames(dataP),
                   abParametersInt = abInt,
                   abParametersGroup = abGroup)
    
    Anota2seqDataSet@qualityControl <- dataOut
    
    return(Anota2seqDataSet)
}

