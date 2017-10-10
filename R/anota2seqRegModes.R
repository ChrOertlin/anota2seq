anota2seqRegModes <- function(Anota2seqDataSet = NULL, mRNASelect = c(TRUE,TRUE)){
    if(is.null(Anota2seqDataSet)){
        stop("Please provide an Anota2seqDataSet")
    }
    if(is.null(anota2seqGetOutputClass(Anota2seqDataSet,"translated mRNA","selected")) == TRUE &
       is.null(anota2seqGetOutputClass(Anota2seqDataSet,"total mRNA","selected")) == TRUE &
       is.null(anota2seqGetOutputClass(Anota2seqDataSet,"translation","selected")) == TRUE &
       is.null(anota2seqGetOutputClass(Anota2seqDataSet,"buffering","selected")) == TRUE){
        stop("No output of anota2seqSelSigGenes detected for any regulation. Please run anota2seqSelSigGenes function for all analysis on the Anota2seqDataSet before running the anota2seqRegModes function ... \n") 
    }
    if(is.null(mRNASelect)){
        stop("Please provide the mRNASelect parameter. Must be one of the following: c(TRUE,TRUE), c(TRUE, FALSE) or c(FALSE,TRUE).\n.")
    }
    if(mRNASelect[1] == FALSE & mRNASelect[2] == FALSE){
        stop("mRNASelect parameter must be one of the following: c(TRUE,TRUE), c(TRUE, FALSE) or c(FALSE,TRUE).\n")
    }
    
    assessmRNA <- TRUE
    # Check for availability of selectedOutputs
    if(is.null(anota2seqGetOutputClass(Anota2seqDataSet,"translated mRNA","selected")) == TRUE |
       is.null(anota2seqGetOutputClass(Anota2seqDataSet,"total mRNA","selected")) == TRUE){
        stop("No anota2seqSelSigGenes output for translated mRNA or total mRNA detected.\n No assessment of mRNA Abundance regulation possible.\n Please run anota2seqSelSigGenes for both transalted mRNA and total mRNA.\n")
        assessmRNA <- FALSE
    }
    if(is.null(anota2seqGetOutputClass(Anota2seqDataSet,"translation","selected"))){
        stop("No selected output for differential translation detected. Please run anota2seqSelSigGenes for analysis of translation.\n")
    }
    if(is.null(anota2seqGetOutputClass(Anota2seqDataSet,"buffering","selected"))){
        stop("No selected output for buffering detected. Please run anota2seqSelSigGenes for analysis of buffering.\n")
    }
    
    usedRVM <- Anota2seqDataSet@selectedTranslation@useRVM
    tmpContrast <- Anota2seqDataSet@translation@usedContrasts
    
    mRNAAbundance <- new("Anota2seqMRNAabundanceOutput",
                         totalmRNA = rep(list(NULL),dim(tmpContrast)[2]),
                         translatedmRNA = rep(list(NULL),dim(tmpContrast)[2]),
                         useRVM = usedRVM,
                         mRNASelect = mRNASelect,
                         regModes = TRUE)
    
    for (cont in 1:dim(tmpContrast)[2]){
        if(assessmRNA == TRUE){
            if(mRNASelect[1] == TRUE & mRNASelect[2] == FALSE){
                ### Get mRNA abundance changes... i.e. differential expression in both polysome-associated mRNA as well as total mRNA (based on FDR cutoff)
                tmpPoly <- anota2seqGetOutput(Anota2seqDataSet,"translated mRNA","selected",cont,usedRVM)
                tmpTotal <- as.data.frame(anota2seqGetOutput(Anota2seqDataSet,"total mRNA","full",cont,usedRVM))
                tmpOverlap <- intersect(rownames(tmpPoly),rownames(tmpTotal))
                tmpPoly <- tmpPoly[tmpOverlap,]
                tmpTotal <- tmpTotal[tmpOverlap,]
                tmpPolyAbund <- tmpPoly[
                    (sign(tmpPoly[,"apvEff"]) == sign(tmpTotal[,"apvEff"])),
                    ,drop=FALSE]
                tmpTotalAbund <- tmpTotal[rownames(tmpPolyAbund),,drop=FALSE]
                mRNAAbundance@translatedmRNA[[cont]] <- tmpPolyAbund
                mRNAAbundance@totalmRNA[[cont]] <- tmpTotalAbund
            }
            
            if(mRNASelect[1] == FALSE & mRNASelect[2] == TRUE){
                ### Get mRNA abundance changes... i.e. differential expression in both polysome-associated mRNA as well as total mRNA (based on FDR cutoff)
                tmpTotal <- anota2seqGetOutput(Anota2seqDataSet,"total mRNA","selected",cont,usedRVM)
                tmpPoly <- as.data.frame(anota2seqGetOutput(Anota2seqDataSet,"translated mRNA","full",cont,usedRVM))
                tmpOverlap <- intersect(rownames(tmpPoly),rownames(tmpTotal))
                tmpPoly <- tmpPoly[tmpOverlap,]
                tmpTotal <- tmpTotal[tmpOverlap,]
                tmpTotalAbund <- tmpTotal[
                    (sign(tmpPoly[,"apvEff"]) == sign(tmpTotal[,"apvEff"])),
                    ,drop=FALSE]
                tmpPolyAbund <- tmpPoly[rownames(tmpTotalAbund),,drop=FALSE]
                mRNAAbundance@translatedmRNA[[cont]] <- tmpPolyAbund
                mRNAAbundance@totalmRNA[[cont]] <- tmpTotalAbund
            }
            if(mRNASelect[1] == TRUE & mRNASelect[2] == TRUE){
                ### Get mRNA abundance changes... i.e. differential expression in both polysome-associated mRNA as well as total mRNA (based on FDR cutoff)
                tmpTotal <- anota2seqGetOutput(Anota2seqDataSet,"total mRNA","selected",cont,usedRVM)
                tmpPoly <- anota2seqGetOutput(Anota2seqDataSet,"translated mRNA","selected",cont,usedRVM)
                tmpOverlap <- intersect(rownames(tmpPoly),rownames(tmpTotal))
                tmpPoly <- tmpPoly[tmpOverlap,]
                tmpTotal <- tmpTotal[tmpOverlap,]
                tmpPolyAbund <- tmpPoly[
                    (sign(tmpPoly[,"apvEff"]) == sign(tmpTotal[,"apvEff"])),
                    ,drop=FALSE]
                
                tmpTotalAbund <- tmpTotal[rownames(tmpPolyAbund),,drop=FALSE]
                mRNAAbundance@translatedmRNA[[cont]] <- tmpPolyAbund
                mRNAAbundance@totalmRNA[[cont]] <- tmpTotalAbund
            }
            ###
        } 
        # Add a column indicating the regulatory mode
        
        ## Translation
        if(usedRVM == TRUE){
            if(!"singleRegMode" %in% colnames(Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]])){
                if(nrow(Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]]) > 0){
                    Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]] <- cbind(Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]],singleRegMode = "translation")
                    Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]]$singleRegMode)
                }
                if(nrow(Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]]) < 1){
                    Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]] <- cbind(Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]],singleRegMode = character(0))
                }
            }
            
            Anota2seqDataSet@selectedTranslation@regModes <- TRUE
        }
        if(usedRVM == FALSE){
            if(nrow(Anota2seqDataSet@selectedTranslation@selectedData[[cont]]) > 0){
                Anota2seqDataSet@selectedTranslation@selectedData[[cont]] <- cbind(Anota2seqDataSet@selectedTranslation@selectedData[[cont]],singleRegMode = "translation")
                Anota2seqDataSet@selectedTranslation@selectedData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTranslation@selectedData[[cont]]$singleRegMode)
            }
            if(nrow(Anota2seqDataSet@selectedTranslation@selectedData[[cont]]) < 1){
                Anota2seqDataSet@selectedTranslation@selectedData[[cont]] <- cbind(Anota2seqDataSet@selectedTranslation@selectedData[[cont]],singleRegMode = character(0))
            }
            
            Anota2seqDataSet@selectedTranslation@regModes <- TRUE
        }
        
        
        # mRNA abundance
        if(is.null(anota2seqGetOutput(Anota2seqDataSet,"translation","selected",cont,usedRVM)) == FALSE & 
           is.null(mRNAAbundance@totalmRNA[[cont]])==FALSE){
            #Add singeRegMode column
            if(!"singleRegMode" %in% colnames(mRNAAbundance@totalmRNA[[cont]])){
                if(nrow(mRNAAbundance@totalmRNA[[cont]]) > 0){
                    mRNAAbundance@totalmRNA[[cont]] <- cbind(mRNAAbundance@totalmRNA[[cont]],singleRegMode = "abundance")
                    mRNAAbundance@totalmRNA[[cont]]$singleRegMode <- as.character(mRNAAbundance@totalmRNA[[cont]]$singleRegMode)
                }
                if(nrow(mRNAAbundance@totalmRNA[[cont]]) < 1){
                    mRNAAbundance@totalmRNA[[cont]] <- cbind(mRNAAbundance@totalmRNA[[cont]],singleRegMode = character(0))
                }
                
            }
            if(nrow(mRNAAbundance@totalmRNA[[cont]]) > 0){
                # Mark genes also found in selectedTranslation as translation as this mode is prioritized
                mRNAAbundance@totalmRNA[[cont]][intersect(rownames(mRNAAbundance@totalmRNA[[cont]]),
                                                          rownames(anota2seqGetOutput(Anota2seqDataSet,"translation","selected",cont,usedRVM))),"singleRegMode"] <- "translation"
            }
            
            #Add singeRegMode column
            if(!"singleRegMode" %in% colnames(mRNAAbundance@translatedmRNA[[cont]])){
                if(nrow(mRNAAbundance@translatedmRNA[[cont]]) > 0){
                    mRNAAbundance@translatedmRNA[[cont]] <- cbind(mRNAAbundance@translatedmRNA[[cont]],singleRegMode = "abundance")
                    mRNAAbundance@translatedmRNA[[cont]]$singleRegMode <- as.character(mRNAAbundance@translatedmRNA[[cont]]$singleRegMode)
                }
                if(nrow(mRNAAbundance@translatedmRNA[[cont]]) < 1){
                    mRNAAbundance@translatedmRNA[[cont]] <- cbind(mRNAAbundance@translatedmRNA[[cont]],singleRegMode = character(0))
                }
            }
            if(nrow(mRNAAbundance@translatedmRNA[[cont]]) > 0){
                # Mark genes also found in selectedTranslation as translation as this mode is prioritized
                mRNAAbundance@translatedmRNA[[cont]][intersect(rownames(mRNAAbundance@translatedmRNA[[cont]]),
                                                               rownames(anota2seqGetOutput(Anota2seqDataSet,"translation","selected",cont,usedRVM))),"singleRegMode"] <- "translation"
            }
        }
        
        
        if(is.null(anota2seqGetOutput(Anota2seqDataSet,"buffering","selected",cont,usedRVM)) == FALSE &
           is.null(anota2seqGetOutput(Anota2seqDataSet,"translation","selected",cont,usedRVM)) == FALSE & 
           is.null(mRNAAbundance@totalmRNA[[cont]]) == FALSE){
            
            
            if(usedRVM == TRUE){
                #Add singleRegMode column
                if(!"singleRegMode" %in% colnames(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]])){
                    if(nrow(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]])>0){
                        Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]] <- cbind(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]],singleRegMode = "buffering")
                        Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]]$singleRegMode)
                    }
                    if(nrow(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]])<1){
                        Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]] <- cbind(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]],singleRegMode = character(0)) 
                    }
                }
                
                Anota2seqDataSet@selectedBuffering@regModes <- TRUE
                if(nrow(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]])>0){
                    # Mark any genes found in translation as translation
                    Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]][
                        intersect(rownames(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]]),
                                  rownames(anota2seqGetOutput(Anota2seqDataSet,"translation","selected",cont,usedRVM))),"singleRegMode"] <- "translation"
                    # mark any genes found in abundance but not in translation as abundance
                    Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]][
                        intersect(rownames(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]]),
                                  rownames(mRNAAbundance@totalmRNA[[cont]])[which(mRNAAbundance@totalmRNA[[cont]][,"singleRegMode"] == "abundance")]),"singleRegMode"] <- "abundance"
                }
            }
            if(usedRVM == FALSE){
                #Add singleRegMode column
                if(!"singleRegMode" %in% colnames(Anota2seqDataSet@selectedBuffering@selectedData[[cont]])){
                    if(nrow(Anota2seqDataSet@selectedBuffering@selectedData[[cont]]) > 0){
                        Anota2seqDataSet@selectedBuffering@selectedData[[cont]] <- cbind(Anota2seqDataSet@selectedBuffering@selectedData[[cont]],singleRegMode = "buffering")
                        Anota2seqDataSet@selectedBuffering@selectedData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedBuffering@selectedData[[cont]]$singleRegMode)
                    }
                    if(nrow(Anota2seqDataSet@selectedBuffering@selectedData[[cont]]) < 1){
                        Anota2seqDataSet@selectedBuffering@selectedData[[cont]] <- cbind(Anota2seqDataSet@selectedBuffering@selectedData[[cont]],singleRegMode = character(0))
                    }
                }
                
                Anota2seqDataSet@selectedBuffering@regModes <- TRUE
                if(nrow(Anota2seqDataSet@selectedBuffering@selectedData[[cont]]) > 0){
                    # Mark any genes found in translation as translation
                    Anota2seqDataSet@selectedBuffering@selectedData[[cont]][
                        intersect(rownames(Anota2seqDataSet@selectedBuffering@selectedData[[cont]]),
                                  rownames(anota2seqGetOutput(Anota2seqDataSet,"translation","selected",cont,usedRVM))),"singleRegMode"] <- "translation"
                    # mark any genes found in abundance but not in translation as abundance
                    Anota2seqDataSet@selectedBuffering@selectedData[[cont]][
                        intersect(rownames(Anota2seqDataSet@selectedBuffering@selectedData[[cont]]),
                                  rownames(mRNAAbundance@totalmRNA[[cont]])[which(mRNAAbundance@totalmRNA[[cont]][,"singleRegMode"] == "abundance")]),"singleRegMode"] <- "abundance"
                }
            }
        }
        ## add regulations to total and poly slots 
        if(usedRVM == TRUE){
            if(!"singleRegMode" %in% colnames(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]])){
                if(nrow(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]]) > 0){
                    Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]] <- cbind(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]],"singleRegMode"= "none")
                    Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]]$singleRegMode)
                }
                if(nrow(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]]) < 1){
                    Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]] <- cbind(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]],"singleRegMode"= character(0))
                    Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]]$singleRegMode)
                }
            }
            if(!"singleRegMode" %in% colnames(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]])){
                if(nrow(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]]) > 0){
                    Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]] <- cbind(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]],"singleRegMode"= "none")
                    Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]]$singleRegMode)
                }
                if(nrow(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]]) < 1){
                    Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]] <- cbind(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]],"singleRegMode"= character(0))
                    Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]]$singleRegMode)
                }
            }
            
            if(nrow(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]]) > 0){
                Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]])
                    ,rownames(
                        Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]])[
                            which(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]][
                                ,"singleRegMode"] == "buffering")])
                    ,"singleRegMode"] <- "buffering"
                
                Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]])
                    ,rownames(
                        mRNAAbundance@totalmRNA[[cont]])[
                            which(mRNAAbundance@totalmRNA[[cont]][
                                ,"singleRegMode"] == "abundance")])
                    ,"singleRegMode"] <- "abundance"
                
                
                Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTotalmRNA@selectedRvmData[[cont]])
                    ,rownames(
                        Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]])[
                            which(Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]][
                                ,"singleRegMode"] == "translation")])
                    ,"singleRegMode"] <- "translation"
            }
            if(nrow(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]]) > 0){
                Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]])
                    ,rownames(
                        Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]])[
                            which(Anota2seqDataSet@selectedBuffering@selectedRvmData[[cont]][
                                ,"singleRegMode"] == "buffering")])
                    ,"singleRegMode"] <- "buffering"
                
                Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]])
                    ,rownames(
                        mRNAAbundance@translatedmRNA[[cont]])[
                            which(mRNAAbundance@translatedmRNA[[cont]][
                                ,"singleRegMode"] == "abundance")])
                    ,"singleRegMode"] <- "abundance"
                
                Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTranslatedmRNA@selectedRvmData[[cont]])
                    ,rownames(
                        Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]])[
                            which(Anota2seqDataSet@selectedTranslation@selectedRvmData[[cont]][
                                ,"singleRegMode"] == "translation")])
                    ,"singleRegMode"] <- "translation"
            }
        }
        
        ## add regulations to total and poly slots 
        if(usedRVM == FALSE){
            if(!"singleRegMode" %in% colnames(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]])){
                if(nrow(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]]) > 0){
                    Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]] <- cbind(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]],"singleRegMode"= "none")
                    Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]]$singleRegMode)
                }
                if(nrow(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]]) < 1){
                    Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]] <- cbind(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]],"singleRegMode"= character(0))
                    Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]]$singleRegMode)
                }
            }
            if(!"singleRegMode" %in% colnames(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]])){
                if(nrow(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]]) > 0){
                    Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]] <- cbind(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]],"singleRegMode"= "none")
                    Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]]$singleRegMode)
                }
                if(nrow(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]]) < 1){
                    Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]] <- cbind(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]],"singleRegMode"= character(0))
                    Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]]$singleRegMode <- as.character(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]]$singleRegMode)
                }
            }
            
            if(nrow(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]]) > 0){
                Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]])
                    ,rownames(
                        Anota2seqDataSet@selectedBuffering@selectedData[[cont]])[
                            which(Anota2seqDataSet@selectedBuffering@selectedData[[cont]][
                                ,"singleRegMode"] == "buffering")])
                    ,"singleRegMode"] <- "buffering"
                
                Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]])
                    ,rownames(
                        mRNAAbundance@totalmRNA[[cont]])[
                            which(mRNAAbundance@totalmRNA[[cont]][
                                ,"singleRegMode"] == "abundance")])
                    ,"singleRegMode"] <- "abundance"
                
                
                Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTotalmRNA@selectedData[[cont]])
                    ,rownames(
                        Anota2seqDataSet@selectedTranslation@selectedData[[cont]])[
                            which(Anota2seqDataSet@selectedTranslation@selectedData[[cont]][
                                ,"singleRegMode"] == "translation")])
                    ,"singleRegMode"] <- "translation"
            }
            if(nrow(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]]) > 0){
                Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]])
                    ,rownames(
                        Anota2seqDataSet@selectedBuffering@selectedData[[cont]])[
                            which(Anota2seqDataSet@selectedBuffering@selectedData[[cont]][
                                ,"singleRegMode"] == "buffering")])
                    ,"singleRegMode"] <- "buffering"
                
                Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]])
                    ,rownames(
                        mRNAAbundance@translatedmRNA[[cont]])[
                            which(mRNAAbundance@translatedmRNA[[cont]][
                                ,"singleRegMode"] == "abundance")])
                    ,"singleRegMode"] <- "abundance"
                
                Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]][intersect(
                    rownames(Anota2seqDataSet@selectedTranslatedmRNA@selectedData[[cont]])
                    ,rownames(
                        Anota2seqDataSet@selectedTranslation@selectedData[[cont]])[
                            which(Anota2seqDataSet@selectedTranslation@selectedData[[cont]][
                                ,"singleRegMode"] == "translation")])
                    ,"singleRegMode"] <- "translation"
            }
            
        }
        
    }
    
    
    if(assessmRNA == TRUE){
        Anota2seqDataSet@mRNAAbundance <- mRNAAbundance
    }
    
    return(Anota2seqDataSet)
    
}
