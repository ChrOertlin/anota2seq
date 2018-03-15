setMethod("anota2seqGetOutput","Anota2seqDataSet",
          function(object, analysis, output,selContrast,getRVM= TRUE) {
              s4MethodChecks(object=object,analysis=analysis,output=output,selContrast=selContrast,getRVM=getRVM,inFunc = "output")
              if(output == "full"){
                  if(getRVM != TRUE & getRVM != FALSE){
                      stop("Please specify RVM status of output with TRUE or FALSE.")
                  }
                  if(analysis == "translated mRNA"){
                      if(is.null(object@translatedmRNA) == FALSE & output == "full") {
                          if(getRVM == FALSE){
                              return(object@translatedmRNA@apvStats[[selContrast]])
                          }
                          if(getRVM == TRUE){
                              return(object@translatedmRNA@apvStatsRvm[[selContrast]])
                          }
                          else{
                              stop("No full translated mRNA output found...")
                          }
                      }
                  }
                  if(analysis == "total mRNA"){
                      if(is.null(object@totalmRNA) == FALSE & output == "full"){
                          if(getRVM == FALSE){
                              return(object@totalmRNA@apvStats[[selContrast]])
                          }
                          if(getRVM == TRUE){
                              return(object@totalmRNA@apvStatsRvm[[selContrast]])
                          }
                          
                          
                          else {
                              stop("No full total mRNA output found ...")
                          }
                      }
                  }
                  
                  if(analysis == "translation"){
                      if(is.null(object@translation) == FALSE & output == "full"){
                          if(getRVM == FALSE){
                              return(object@translation@apvStats[[selContrast]])
                          }
                          if(getRVM == TRUE){
                              return(object@translation@apvStatsRvm[[selContrast]])
                          }
                          else{
                              stop("No full translation output found ...")
                          }
                      }
                  }
                  
                  if(analysis == "buffering"){
                      if(is.null(object@buffering) == FALSE & output == "full"){
                          if(getRVM == FALSE){
                              return(object@buffering@apvStats[[selContrast]])
                          }
                          if(getRVM == TRUE){
                              return(object@buffering@apvStatsRvm[[selContrast]])
                          }
                          else{
                              stop("No full buffering output found ...")
                          }
                      }
                  }
              }# if output
              if(output =="selected"){
                  
                  if(analysis =="translated mRNA"){
                      if(is.null(object@selectedTranslatedmRNA) == FALSE){
                          if(getRVM != object@selectedTranslatedmRNA@useRVM){
                              stop(paste("Can only retrieve selected data with RVM status used in anota2seqSelSigGenes.\n",
                                         "getRVM parameter is set to ",
                                         getRVM,
                                         " while anota2seqSelSigGenes was run with parameter useRVM ",
                                         object@selectedTranslatedmRNA@useRVM,".\nPlease provide correct RVM parameter.\n"))
                          }
                          getRVM <- object@selectedTranslatedmRNA@useRVM
                          if(getRVM == FALSE){
                              return(object@selectedTranslatedmRNA@selectedData[[selContrast]])
                          }
                          if(getRVM == TRUE){
                              return(object@selectedTranslatedmRNA@selectedRvmData[[selContrast]])
                          }
                          
                          else {
                              stop("No selected translated mRNA output found...")
                          }
                      }
                  }
                  if(analysis == "total mRNA"){
                      if(is.null(object@selectedTotalmRNA) == FALSE){
                          if(getRVM != object@selectedTotalmRNA@useRVM){
                              stop(paste("Can only retrieve selected data with RVM status used in anota2seqSelSigGenes.\n",
                                         "getRVM parameter is set to ",
                                         getRVM,
                                         " while anota2seqSelSigGenes was run with parameter useRVM ",
                                         object@totalmRNA@useRVM,".\nPlease provide correct RVM parameter.\n"))
                          }
                          getRVM <- object@selectedTotalmRNA@useRVM
                          if(getRVM == FALSE){
                              return(object@selectedTotalmRNA@selectedData[[selContrast]])
                          } 
                          if(getRVM == TRUE){
                              return(object@selectedTotalmRNA@selectedRvmData[[selContrast]])
                          }
                          
                          else{
                              stop("No selected total mRNA output found ...")
                          }
                      }
                  }
                  if(analysis == "translation"){
                      if(is.null(object@selectedTranslation) == FALSE){
                          if(getRVM != object@selectedTranslation@useRVM){
                              stop(paste("Can only retrieve selected data with RVM status used in anota2seqSelSigGenes.\n",
                                         "getRVM parameter is set to ",
                                         getRVM,
                                         " while anota2seqSelSigGenes was run with parameter useRVM ",
                                         object@selectedTranslation@useRVM,".\nPlease provide correct RVM parameter.\n"))
                          }
                          getRVM <- object@selectedTranslation@useRVM
                          if(getRVM == FALSE){
                              return(object@selectedTranslation@selectedData[[selContrast]])
                          }
                          if(getRVM == TRUE){
                              return(object@selectedTranslation@selectedRvmData[[selContrast]])
                          }
                          else {
                              stop("No selected translation output found ...")
                          }
                      }
                  }
                  
                  if(analysis == "buffering"){
                      if(is.null(object@selectedBuffering) == FALSE){
                          if(getRVM != object@selectedBuffering@useRVM){
                              stop(paste("Can only retrieve selected data with RVM status used in anota2seqSelSigGenes.\n",
                                         "getRVM parameter is set to ",
                                         getRVM,
                                         " while anota2seqSelSigGenes was run with parameter useRVM ",
                                         object@selectedBuffering@useRVM,".\nPlease provide correct RVM parameter.\n"))
                          }
                          getRVM <- object@selectedBuffering@useRVM
                          if(getRVM == FALSE){
                              return(object@selectedBuffering@selectedData[[selContrast]])
                          }
                          if(getRVM == TRUE){
                              return(object@selectedBuffering@selectedRvmData[[selContrast]])
                          }
                          else {
                              stop("No selected buffering output found ...")
                          }
                      }
                  }
                  
                  if(analysis == "mRNA abundance"){
                      if(is.null(object@mRNAAbundance) == FALSE){
                          if(getRVM != object@mRNAAbundance@useRVM){
                              stop(paste("Can only retrieve selected data with RVM status used in anota2seqSelSigGenes.\n",
                                         "getRVM parameter is set to ",
                                         getRVM,
                                         " while anota2seqSelSigGenes was run with parameter useRVM ",
                                         object@mRNAAbundance@useRVM,".\nPlease provide correct RVM parameter.\n"))
                          }
                          mRNASelect <- object@mRNAAbundance@mRNASelect
                          if(mRNASelect[1] == TRUE & mRNASelect[2] == TRUE){
                              return(list("abundance total mRNA" =object@mRNAAbundance@totalmRNA[[selContrast]],
                                          "abundance translated mRNA"= object@mRNAAbundance@translatedmRNA[[selContrast]])
                              )
                          }
                          if(mRNASelect[1] == TRUE & mRNASelect[2] == FALSE){
                              return(object@mRNAAbundance@translatedmRNA[[selContrast]])
                          }
                          if(mRNASelect[1] == FALSE & mRNASelect[2] == TRUE){
                              return(object@mRNAAbundance@totalmRNA[[selContrast]])
                          }
                      }
                      if(is.null(object@mRNAAbundance) == TRUE){
                          stop("No mRNA abundance genes output found ...")
                      }
                  }
              }#if selected
              if(output =="regModes"){
                  if(analysis %in% c("translated mRNA","total mRNA")){
                      stop("When selecting regModes output, analysis parameter can only be set to: translation, buffering or mRNA abundance.\n")
                  }
                  if(analysis == "translation"){
                      if(is.null(object@selectedTranslation) == FALSE){
                          if(object@selectedTranslation@regModes == FALSE){
                              stop("No anota2seqRegModes output found for translation analysis.\nPlease run the anota2seqRegModes function on the Anota2seqDataSet.\n")
                          }
                          if(getRVM != object@selectedTranslation@useRVM){
                              stop(paste("Can only retrieve selected data with RVM status used in anota2seqSelSigGenes.\n",
                                         "getRVM parameter is set to ",
                                         getRVM,
                                         " while anota2seqSelSigGenes was run with parameter useRVM ",
                                         object@selectedTranslation@useRVM,".\nPlease provide correct RVM parameter.\n"))
                          }
                          getRVM <- object@selectedTranslation@useRVM
                          if(getRVM == FALSE){
                              return(
                                  object@selectedTranslation@selectedData[[selContrast]][which(object@selectedTranslation@selectedData[[selContrast]][,"singleRegMode"] == "translation"),]
                              )
                          }
                          if(getRVM == TRUE){
                              return(
                                  object@selectedTranslation@selectedRvmData[[selContrast]][which(object@selectedTranslation@selectedRvmData[[selContrast]][,"singleRegMode"] == "translation"),]
                              )
                          }
                      }
                  }
                  
                  if(analysis == "buffering"){
                      if(is.null(object@selectedBuffering) == FALSE){
                          if(object@selectedBuffering@regModes == FALSE){
                              stop("No anota2seqRegModes output found for buffering analysis.\nPlease run the anota2seqRegModes function on the Anota2seqDataSet.\n")
                          }
                          if(getRVM != object@selectedBuffering@useRVM){
                              stop(paste("Can only retrieve selected data with RVM status used in anota2seqSelSigGenes.\n",
                                         "getRVM parameter is set to ",
                                         getRVM,
                                         " while anota2seqSelSigGenes was run with parameter useRVM ",
                                         object@selectedBuffering@useRVM,".\nPlease provide correct RVM parameter.\n"))
                          }
                          getRVM <- object@selectedBuffering@useRVM
                          if(getRVM == FALSE){
                              return(
                                  object@selectedBuffering@selectedData[[selContrast]][which(object@selectedBuffering@selectedData[[selContrast]][,"singleRegMode"] == "buffering"),]
                              )
                          }
                          if(getRVM == TRUE){
                              return(
                                  object@selectedBuffering@selectedRvmData[[selContrast]][which(object@selectedBuffering@selectedRvmData[[selContrast]][,"singleRegMode"] == "buffering"),]
                              )
                          }
                      }
                  }
                  
                  if(analysis == "mRNA abundance"){
                      if(is.null(object@mRNAAbundance) == FALSE){
                          if(getRVM != object@mRNAAbundance@useRVM){
                              stop(paste("Can only retrieve selected data with RVM status used in anota2seqSelSigGenes.\n",
                                         "getRVM parameter is set to ",
                                         getRVM,
                                         " while anota2seqSelSigGenes was run with parameter useRVM ",
                                         object@mRNAAbundance@useRVM,".\nPlease provide correct RVM parameter.\n"))
                          }
                          if(object@mRNAAbundance@regModes == FALSE){
                              stop("No anota2seqRegModes output found for mRNA abundance analysis.\nPlease run the anota2seqRegModes function on the Anota2seqDataSet.\n")
                          }
                          mRNASelect <- object@mRNAAbundance@mRNASelect
                          if(mRNASelect[1] == TRUE & mRNASelect[2] == TRUE){
                              return(list("abundance total mRNA" =object@mRNAAbundance@totalmRNA[[selContrast]][which(object@mRNAAbundance@totalmRNA[[selContrast]][,"singleRegMode"] == "abundance"),],
                                          "abundance translated mRNA"= object@mRNAAbundance@translatedmRNA[[selContrast]][which(object@mRNAAbundance@translatedmRNA[[selContrast]][,"singleRegMode"] == "abundance"),])
                              )
                          }
                          if(mRNASelect[1] == TRUE & mRNASelect[2] == FALSE){
                              return(object@mRNAAbundance@translatedmRNA[[selContrast]][which(object@mRNAAbundance@translatedmRNA[[selContrast]][,"singleRegMode"] == "abundance"),])
                          }
                          if(mRNASelect[1] == FALSE & mRNASelect[2] == TRUE){
                              return(object@mRNAAbundance@totalmRNA[[selContrast]][which(object@mRNAAbundance@totalmRNA[[selContrast]][,"singleRegMode"] == "abundance"),])
                          }
                      }
                  }
              }#if regModes
          })
setMethod("anota2seqGetQualityControl","Anota2seqDataSet",
          function(object){
              s4MethodChecks(object=object,inFunc = "NA")
              if(is.null(object@qualityControl) == FALSE){
                  return(list(
                      omniIntStats = object@qualityControl@omniIntStats,
                      omniGroupStats = object@qualityControl@omniGroupStats,
                      groupIntercepts = object@qualityControl@groupIntercepts,
                      correctionMethod = object@qualityControl@correctionMethod,
                      dsfSummary = object@qualityControl@dsfSummary,
                      dfbetas = object@qualityControl@dfbetas,
                      residuals = object@qualityControl@residuals,
                      fittedValues = object@qualityControl@fittedValues,
                      phenoClasses = object@qualityControl@phenoClasses,
                      sampleNames = object@qualityControl@sampleNames,
                      abParametersInt = object@qualityControl@abParametersInt,
                      abParametersGroup = object@qualityControl@abParametersGroup  
                  ))
              }
              
              if(is.null(object@qualityControl) == TRUE){
                  stop("No quality control detected. Please run the anota2seqPerformQC function on the Anota2seqDataSet.\n")
              }
              
          })
setMethod("anota2seqGetResidOutlierTest","Anota2seqDataSet",
          function(object){
              s4MethodChecks(object=object,inFunc = "NA")
              
              if(is.null(object@residOutlierTest) == FALSE){
                  return(list(              
                      confInt = object@residOutlierTest@confInt,
                      rnormIter = object@residOutlierTest@rnormIter,
                      outlierMatrixLog = object@residOutlierTest@outlierMatrixLog,
                      meanOutlierPerIteration = object@residOutlierTest@meanOutlierPerIteration,
                      obtainedComparedToExpected = object@residOutlierTest@obtainedComparedToExpected,
                      nExpected = object@residOutlierTest@nExpected,
                      nObtained= object@residOutlierTest@nObtained
                  ))
              }
              if(is.null(object@residOutlierTest) == TRUE){
                  stop("No residOutlierTest detected. Please run the anota2seqResidOutlierTest function on the Anota2seqDataSet. \n")
              }
          })
setMethod("anota2seqGetDeltaData","Anota2seqDataSet",
          function(object, output, analysis, selContrast){
              s4MethodChecks(object=object,output=output,analysis=analysis,selContrast=selContrast,inFunc = "delta")
              if(output == "full"){
                  if(is.null(object@deltaData)){
                      stop("No deltaData found. Please run the anota2seqAnalyze function on the Anota2seqDataSet.\n")
                  }
                  if(analysis == "translation"){
                      return(object@deltaData[[selContrast]][,c("deltaP","deltaPT"),drop=FALSE])
                  }
                  if(analysis == "buffering"){
                      return(object@deltaData[[selContrast]][,c("deltaT","deltaTP"),drop=FALSE])
                  }
                  if(analysis == "translated mRNA"){
                      return(as.matrix(object@deltaData[[selContrast]][,"deltaP",drop=FALSE]))
                  }
                  
                  if(analysis == "total mRNA"){
                      return(as.matrix(object@deltaData[[selContrast]][,"deltaT",drop=FALSE]))
                  }
              }
              
              if(output == "selected"){
                  if(is.null(anota2seqGetOutputClass(object,output="selected",analysis=analysis)) == TRUE){
                      stop("selected deltaData (i.e. output set to selected) can only be retrieved if anota2seqSelSigGenes has been run for specified analysis.\nPlease run anota2seqSelSigGenes on the Anota2seqDataSet.\n")
                  }
                  if(is.null(anota2seqGetOutputClass(object,output="selected",analysis=analysis)) == FALSE){
                      tmpDat <- anota2seqGetOutputClass(object,output="selected",analysis=analysis)
                      return(tmpDat@deltaData[[selContrast]])
                  }
              }
          })
setMethod("anota2seqGetThresholds","Anota2seqDataSet",
          function(object = NULL, analysis, selContrast){
              s4MethodChecks(object=object,analysis=analysis,selContrast=selContrast,inFunc = "tresholds")
              if(is.null(anota2seqGetOutputClass(object,output="selected",analysis=analysis)) == TRUE){
                  stop("thresholds can only be retrieved if anota2seqSelSigGenes has been run for specified analysis.\nPlease run anota2seqSelSigGenes on the Anota2seqDataSet.\n")
              }
              if(is.null(anota2seqGetOutputClass(object,output="selected",analysis=analysis)) == FALSE){
                  tmpDat <- anota2seqGetOutputClass(object,output="selected",analysis=analysis)
                  return(tmpDat@usedThresholds[[selContrast]])
              }
          })
setMethod("anota2seqGetNormalizedData","Anota2seqDataSet",
          definition = function(object){ 
              s4MethodChecks(object=object,inFunc = "NA")
              if(is.null(object@dataP) | is.null(object@dataT)){
                  stop("No dataP or dataT found in Anota2seqDataSet.\n")
              }
              return(list(dataP = object@dataP,dataT = object@dataT))
          })
setMethod("anota2seqGetCovariates","Anota2seqDataSet",
          function(object) {
              s4MethodChecks(object=object,inFunc = "NA")
              if(is.null(object@phenoVec)){
                  stop("No phenoVec found in Anota2seqDataSet.\n")
              }
              return(list(phenoVec = object@phenoVec,batchVec=object@batchVec))
          })
setMethod("anota2seqGetContrasts","Anota2seqDataSet",
          function(object){
              s4MethodChecks(object=object,inFunc = "NA")
              if(is.null(object@contrasts)){
                  warning("No contrasts found. Contrasts are set with the anota2seqAnalyze or anota2seqRun function.\n")
              }
              return(object@contrasts)})
setMethod("anota2seqPlotFC","Anota2seqDataSet",
          function(object,visualizeRegModes="all",selContrast, contrastName = NULL,fileStem = "ANOTA2SEQ_FoldchangePlot",plotToFile = TRUE, myYlim = NULL, myXlim = NULL,...){    
              message("Creating Fold-change plots.\n")
              if(is.null(object@buffering)&is.null(object@translation)&is.null(object@translatedmRNA)&is.null(object@totalmRNA)){
                  stop("No anota2seqAnalyze output detected. Please run the anota2seqAnalyze function before using the anota2seqPlotFC function.")
              }
              if(visualizeRegModes == "translation"){
                  if(is.null(anota2seqGetOutputClass(object,analysis = "translation","selected"))){
                      stop("No anota2seqAnalyse output for translation detected. Please run  anota2seqSelSigGenes with analysis parameter set to translation before proceeding.\n")
                  }
              }
              if(visualizeRegModes == "buffering"){
                  if(is.null(anota2seqGetOutputClass(object,analysis = "buffering","selected"))){
                      stop("No anota2seqAnalyse output for buffering detected. Please run anota2seqSelSigGenes with analysis parameter set to buffering before proceeding.\n")
                  }
              }
              if(visualizeRegModes == "all"){
                  if(anota2seqGetOutputClass(object,analysis = "translation","selected")@regModes == FALSE){
                      stop("No regModes found. Please run the anota2seqRegModes function on the object before generating fold change plots.\n")
                  }
              }
              if(!is.null(myXlim) | !is.null(myXlim)){
                  if(is.null(myXlim)){
                      stop("myYlim is set but not myXlim. Please specify both parameters or set both to NULL.\n")
                  }
                  if(is.null(myYlim)){
                      stop("myXlim is set but not myYlim. Please specify both parameters or set both to NULL.\n")
                  }
              }
              
              if(!is.null(contrastName)){
                  if(length(contrastName) != length(selContrast)){
                      stop("Number of contrast names do not match the number of selected contrasts.\nPlease supply a contrast name for each selected contrast.\n")
                  }
              }
              if(is.null(contrastName)){
                  contrastName <- paste("contrast ",1:length(selContrast))
              }
              
              s4MethodChecks(object=object,
                             selContrast=selContrast,
                             visualizeRegModes=visualizeRegModes,
                             plotToFile=plotToFile,
                             inFunc = "anota2seqPlotFC")
              
              cols <- c(RColorBrewer::brewer.pal(8,"Reds")[c(4,8)],
                        RColorBrewer::brewer.pal(8,"Blues")[c(4,8)],
                        RColorBrewer::brewer.pal(8,"Greens")[c(4,8)])
              
              names(cols) <- c("translation up","translation down","buffering down","buffering up","mRNA abundance up","mRNA abundance down")
              tmpContrasts <- object@contrasts
              
              graphArgs <- list(...)
              if(length(graphArgs) < 1){
                  graphArgs <- list(mar=c(5,5,3,2)+0.1)
              }
              
              if(length(graphArgs) > 0){
                  if(!"mar" %in% names(graphArgs)){
                      graphArgs[["mar"]] <- c(5,5,3,2) + 0.1
                  }
              }
              
              par(graphArgs)
              for(i in 1:length(selContrast)){
                  deltaP <- object@deltaData[[selContrast[i]]][,"deltaP"]
                  deltaT <- object@deltaData[[selContrast[i]]][,"deltaT"]
                  
                  if(plotToFile ==TRUE){
                      pdf(paste(fileStem, gsub(" ","",contrastName[i]) ,".pdf",sep=""))
                  }
                  # plot foldchanges
                  maxVal <- max(abs(cbind(deltaT,deltaP)),na.rm = TRUE)
                  
                  if(is.null(myXlim) & is.null(myYlim)){              
                      myXlim <- c(maxVal*-1,maxVal)
                      myYlim <- c(maxVal*-1,maxVal)
                  }
                  
                  if(!is.null(myXlim) & !is.null(myYlim)){  
                      myXlim <- myXlim
                      myYlim <- myYlim
                  }
                  
                  
                  plot(x=deltaT,y=deltaP, pch=19, cex=.8,
                       xlab=paste("total mRNA log2FC\n(",contrastName[i],")", sep = ""),
                       ylab = paste("translated mRNA log2FC\n(",contrastName[i],")", sep = ""),
                       ylim = myYlim,
                       xlim = myXlim,
                       col="grey")
                  abline(h=0,lty=2)
                  abline(v=0,lty=2)
                  abline(a = 0, b = 1, lty = 2)
                  
                  if(visualizeRegModes == "all"){
                      
                      useRVM <- object@selectedTranslation@useRVM
                      if(!is.null(object@mRNAAbundance@totalmRNA[[selContrast[i]]])){
                          # mRNA abundance up
                          points(x= deltaT[rownames(object@mRNAAbundance@totalmRNA[[selContrast[i]]])[
                              which(object@mRNAAbundance@totalmRNA[[selContrast[i]]][,"apvEff"] > 0 & object@mRNAAbundance@totalmRNA[[selContrast[i]]][,"singleRegMode"] == "abundance")
                              ]],
                              y= deltaP[rownames(object@mRNAAbundance@totalmRNA[[selContrast[i]]])[
                                  which(object@mRNAAbundance@totalmRNA[[selContrast[i]]][,"apvEff"] > 0& object@mRNAAbundance@totalmRNA[[selContrast[i]]][,"singleRegMode"] == "abundance")
                                  ]],
                              pch=19,col=cols["mRNA abundance up"],cex=.5)
                          #mRNA abundance down
                          points(x= deltaT[rownames(object@mRNAAbundance@totalmRNA[[selContrast[i]]])[
                              which(object@mRNAAbundance@totalmRNA[[selContrast[i]]][,"apvEff"] < 0& object@mRNAAbundance@totalmRNA[[selContrast[i]]][,"singleRegMode"] == "abundance")
                              ]],
                              y= deltaP[rownames(object@mRNAAbundance@totalmRNA[[selContrast[i]]])[
                                  which(object@mRNAAbundance@totalmRNA[[selContrast[i]]][,"apvEff"] < 0& object@mRNAAbundance@totalmRNA[[selContrast[i]]][,"singleRegMode"] == "abundance")
                                  ]],
                              pch=19,col=cols["mRNA abundance down"],cex=.5)
                          
                      }
                      if(!is.null(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM))& !is.null(anota2seqGetOutput(object,"translation","full",selContrast[i],useRVM))){
                          
                          #up regulated differential translation
                          points(x= deltaT[rownames(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM))[
                              which(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"apvEff"] >= 0 & anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"singleRegMode"] == "translation"  )]],
                              y= deltaP[rownames(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM))[
                                  which(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"apvEff"] >= 0& anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"singleRegMode"] == "translation"   )]],
                              pch=19,col=cols["translation up"],cex=.5)
                          # down regulated differential translation
                          points(x= deltaT[rownames(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM))[
                              which(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"apvEff"] < 0 & anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"singleRegMode"] == "translation"  )]],
                              y= deltaP[rownames(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM))[
                                  which(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"apvEff"] < 0& anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"singleRegMode"] == "translation"   )]],
                              pch=19,col=cols["translation down"],cex=.5)
                      }
                      
                      if(!is.null(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM))& !is.null(anota2seqGetOutput(object,"buffering","full",selContrast[i],useRVM))){
                          
                          #up regulated differential buffering
                          points(x= deltaT[rownames(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM))[
                              which(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"apvEff"] < 0 & anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"singleRegMode"] == "buffering")]],
                              y= deltaP[rownames(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM))[
                                  which(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"apvEff"]< 0& anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"singleRegMode"] == "buffering")]],
                              pch=19,col=cols["buffering up"],cex=.5)
                          # down regulated differential buffering
                          points(x= deltaT[rownames(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM))[
                              which(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"apvEff"] > 0 & anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"singleRegMode"] == "buffering")]],
                              y= deltaP[rownames(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM))[
                                  which(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"apvEff"] > 0 & anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"singleRegMode"] == "buffering")]],
                              pch=19,col=cols["buffering down"],cex=.5)
                      }
                      
                      transl <- anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[which(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"singleRegMode"] == "translation"),]
                      buff <- anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[which(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"singleRegMode"] == "buffering"),]
                      abund <- rownames(object@mRNAAbundance@totalmRNA[[selContrast[i]]])[which(object@mRNAAbundance@totalmRNA[[selContrast[i]]][,"singleRegMode"] == "abundance")]
                  }
                  if(visualizeRegModes == "translation"){
                      
                      useRVM <- object@selectedTranslation@useRVM
                      
                      #up regulated differential translation
                      points(x= deltaT[rownames(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM))[
                          which(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"apvEff"] >= 0)]],
                          y= deltaP[rownames(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM))[
                              which(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"apvEff"] >= 0)]],
                          pch=19,col=cols["translation up"],cex=.5)
                      # down regulated differential translation
                      points(x= deltaT[rownames(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM))[
                          which(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"apvEff"] < 0 )]],
                          y= deltaP[rownames(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM))[
                              which(anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)[,"apvEff"] < 0)]],
                          pch=19,col=cols["translation down"],cex=.5)
                      buff <- NULL
                      abund <- NULL
                      transl <- anota2seqGetOutput(object,"translation","selected",selContrast[i],useRVM)
                  }
                  if(visualizeRegModes == "buffering"){
                      
                      useRVM <- object@selectedBuffering@useRVM
                      
                      #bufferin down
                      points(x= deltaT[rownames(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM))[
                          which(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"apvEff"] > 0)]],
                          y= deltaP[rownames(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM))[
                              which(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"apvEff"] > 0)]],
                          pch=19,col=cols["buffering down"],cex=.5)
                      #buffering up
                      points(x= deltaT[rownames(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM))[
                          which(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"apvEff"] < 0)]],
                          y= deltaP[rownames(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM))[
                              which(anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)[,"apvEff"] < 0)]],
                          pch=19,col=cols["buffering up"],cex=.5)
                      transl <- NULL
                      abund <- NULL
                      buff <- anota2seqGetOutput(object,"buffering","selected",selContrast[i],useRVM)
                  }
                  
                  if(visualizeRegModes%in%c("all","buffering","translation")){
                      legendVec <- c(paste("Translation up (", sum(transl[, "apvEff"] > 0), ")", sep = ""),
                                     paste("Translation down (", sum(transl[, "apvEff"] < 0), ")", sep = ""),
                                     paste("Buffered (mRNA up) (", sum(buff[, "apvEff"] > 0), ")", sep = ""),
                                     paste("Buffered (mRNA down) (", sum(buff[, "apvEff"] < 0), ")", sep = ""),
                                     paste("mRNA abundance up (", sum(deltaP[abund] > 0), ")", sep = ""),
                                     paste("mRNA abundance down (", sum(deltaP[abund] < 0), ")", sep = ""))
                      
                      anyGeneVec <- c(rep(!is.null(transl), 2), rep(!is.null(buff), 2),
                                      rep(!is.null(abund),2))
                      
                      
                      if(sum(anyGeneVec) != 0){
                          legend("topleft", legend = legendVec[anyGeneVec], col = cols[anyGeneVec], pch = 19)
                      }
                      
                      
                      
                      thresholdsP <- NULL
                      thresholdsT <- NULL
                      if(visualizeRegModes %in% c("all", "translation")){
                          if(is.null( anota2seqGetThresholds(object,analysis = "translation",selContrast = selContrast[i]) ) == FALSE){
                              thresholdsP <- anota2seqGetThresholds(object,analysis = "translation",selContrast = selContrast[i]) 
                          }
                      }
                      
                      if(visualizeRegModes %in% c("all", "buffering")){
                          if(is.null( anota2seqGetThresholds(object,analysis = "buffering",selContrast = selContrast[i]) ) == FALSE){
                              thresholdsT <- anota2seqGetThresholds(object,analysis = "buffering",selContrast = selContrast[i]) 
                          }
                      }
                      if(is.null(thresholdsT$selDeltaT) == FALSE){
                          vert <- c(-thresholdsT$selDeltaT, thresholdsT$selDeltaT)
                          abline(v = vert, lty=1, lwd=2)
                      }
                      if(is.null(thresholdsP$selDeltaP) == FALSE){
                          horiz <- c(-thresholdsP$selDeltaP, thresholdsP$selDeltaP)
                          abline(h = horiz, lty=1, lwd=2)
                      }
                      if(is.null(thresholdsP$selDeltaPT) == FALSE){
                          
                          diag <- c(-thresholdsP$selDeltaPT, thresholdsP$selDeltaPT)
                          abline(a = diag[1], b = 1, lty = 1, lwd = 2)
                          abline(a = diag[2], b = 1, lty = 1, lwd = 2)
                      }
                  }
                  
                  if(plotToFile==TRUE){
                      dev.off()
                  }
              }
              
          })

setMethod("anota2seqPlotPvalues","Anota2seqDataSet",
          function(object,useRVM = TRUE,selContrast,contrastName = NULL,myBw = 0.05,plotToFile=TRUE, fileStem = "ANOTA2SEQ_pvalue_density", ...){ 
              message("Creating pvalue and FDR density plots.\n") 
              if(is.null(object@buffering)&is.null(object@translation)&is.null(object@translatedmRNA)&is.null(object@totalmRNA)){
                  stop("No anota2seqAnalyze output detected. Please run the anota2seqAnalyze function before using the anota2seqPlotFC function.")
              }
              s4MethodChecks(object=object,selContrast=selContrast,useRVM=useRVM,myBw=myBw,plotToFile=plotToFile,inFunc = "anota2seqPlotPvalues")
              #Check whether to plot RVM values or non-RVM values.
              tmpCols <- c("apvP","apvPAdj")
              if(useRVM == TRUE){
                  tmpCols <- c("apvRvmP","apvRvmPAdj")
              }
              tmpContrasts <- object@contrasts
              # make a list of lists for all created contrasts...
              plotList <- rep(list(NULL),length(selContrast))
              names(plotList) <- paste("contrast",1:length(selContrast),sep="")
              regulations <- c("translated mRNA", "total mRNA","translation","buffering")
              nameList <- rep(list(NULL),length(selContrast))
              graphArgs <- list(...)
              for( contr in 1:length(selContrast)){
                  for(regs in 1:length(regulations)){
                      if(is.null(anota2seqGetOutput(object,regulations[regs],output = "full",selContrast  = selContrast[contr],getRVM = useRVM)) == FALSE){
                          plotList[[contr]][[regs]] <- anota2seqGetOutput(object,regulations[regs],output = "full",selContrast = selContrast[contr],getRVM = useRVM)
                          nameList[[contr]] <- c(nameList[[contr]],regulations[regs])
                      }
                  }
                  names(plotList[[contr]]) <- regulations
              }
              
              tmpColours <- c(RColorBrewer::brewer.pal(8,"Reds")[8],
                              RColorBrewer::brewer.pal(8,"Blues")[8],
                              RColorBrewer::brewer.pal(8,"Set1")[c(4,5)])
              names(tmpColours) <- c("translation","buffering","total mRNA","translated mRNA")
              
              # #Create a densityPlot per analysis and contrast
              # #First split the plots by contrast
              # # TO DO: FIX the ylim issue... (i.e. get some max value to set ylim...) 
              # # First get densities then plot then - retrieve max value....
              fdrDens <- rep(list(NULL),length(selContrast))
              pvalDens <-rep(list(NULL),length(selContrast))
              fdrMax <- vector("numeric")
              pvalMax <- vector("numeric")
              
              if(!is.null(contrastName)){
                  if(length(contrastName) != length(selContrast)){
                      stop("More contrasts selected (selContrast) than contrast names supplied (contrastName).\nPlease supply a contrast name for each selected contrast.\n")
                  }
              }
              if(is.null(contrastName)){
                  contrastName <- paste("contrast ",1:length(selContrast),sep="")
              }
              
              for(cont in 1:length(selContrast)){
                  for(names in 1:length(plotList[[cont]])){
                      
                      if(is.null(plotList[[cont]][[nameList[[cont]][names]]])==FALSE){
                          
                          pvalDens[[cont]][[names]] <- density(as.numeric(plotList[[cont]][[
                              nameList[[cont]][
                                  names]]][
                                      ,tmpCols[1]]),bw=myBw)
                          
                          fdrDens[[cont]][[names]] <- density(as.numeric(plotList[[cont]][[
                              nameList[[cont]][
                                  names]]][
                                      ,tmpCols[2]]),bw=myBw)
                      }
                  }
                  names(pvalDens[[cont]]) <- nameList[[cont]]
                  names(fdrDens[[cont]]) <- nameList[[cont]]
              }
              par(graphArgs)
              for(cont in 1:length(selContrast)){
                  maxFDR <-  max(unlist(lapply(fdrDens[[cont]],function(x) max(x$y))))
                  maxPval <- max(unlist(lapply(pvalDens[[cont]],function(x) max(x$y))))
                  if(plotToFile == TRUE){
                      pdf(paste(fileStem, "_",gsub(" ","",contrastName[cont]),".pdf",sep=""))
                  }
                  
                  for(reg in 1:length(pvalDens[[cont]])){
                      if(reg == 1){
                          plot(pvalDens[[cont]][[reg]],
                               main = contrastName[cont],lwd=2,ylim=c(0,maxPval+1),xlim=c(-0.2,1.2),
                               xlab = "P-value",col =tmpColours[names(pvalDens[[cont]])[reg]])
                      }
                      else{
                          lines(pvalDens[[cont]][[reg]],lwd=2,col =tmpColours[names(pvalDens[[cont]])[reg]])
                      }
                      
                  }
                  legend("top",col=tmpColours[names(pvalDens[[cont]])],lty=1,lwd=2,legend = names(pvalDens[[cont]]))
                  for(reg in 1:length(fdrDens[[cont]])){
                      if(reg == 1){
                          plot(fdrDens[[cont]][[reg]],
                               main = contrastName[cont],lwd=2,ylim=c(0,maxFDR+1),xlim=c(-0.2,1.2),
                               xlab = "FDR",col =tmpColours[names(fdrDens[[cont]])[reg]])
                      }
                      else{
                          lines(fdrDens[[cont]][[reg]],lwd=2,col =tmpColours[names(fdrDens[[cont]])[reg]])
                      }
                      
                  }
                  legend("top",col=tmpColours[names(pvalDens[[cont]])],
                         lty=1,lwd=2,legend = names(pvalDens[[cont]]))
                  
                  if(plotToFile == TRUE){
                      dev.off()
                  }
              }
          })

setMethod("anota2seqPlotGenes","Anota2seqDataSet",
          function(object,selContrast,analysis,geneNames = NULL,plotToFile = TRUE,fileStem="ANOTA2SEQ_significantGenes_plot"){
              s4MethodChecks(object=object,selContrast=selContrast,analysis = analysis,plotToFile=plotToFile,inFunc = "anota2seqPlotGenes")
              if(is.null(anota2seqGetOutputClass(object,analysis,"selected"))){
                  stop("No anota2seqSelSigGenes output in Anota2seqDataSet found.\n Please run the anota2seqSelSigGenes function on the Anota2seqDataSet.\n")
              }
              
              phenoVec <- object@phenoVec
              anota2seqSigObj <- anota2seqGetOutputClass(object,analysis = analysis,output = "full")
              if (analysis == "translation"){
                  useRVM <- object@selectedTranslation@useRVM
                  useIds <- rownames(anota2seqGetOutput(object,analysis,"selected",selContrast))
                  dataX <- object@dataT
                  dataY <- object@dataP
                  labx = "total mRNA"
                  laby = "translated mRNA"
              } else if (analysis == "buffering"){
                  useRVM <- object@selectedTranslation@useRVM
                  useIds <- rownames(anota2seqGetOutput(object,analysis,"selected",selContrast))
                  dataX <- object@dataP
                  dataY <- object@dataT
                  labx = "translated mRNA"
                  laby = "total mRNA"
              }
              
              useGeneNames <- FALSE
              if (is.null(geneNames) == FALSE) {
                  useGeneNames <- TRUE
                  tmp <- rownames(geneNames)
                  geneNames <- as.vector(geneNames)
                  names(geneNames) <- tmp
              }
              
              if(plotToFile == TRUE){
                  pdf(paste(fileStem, ".pdf", sep = ""), width = 12,height = 12)
              }
              par(mfrow = c(3, 3))
              for (i in 1:length(useIds)) {
                  tmpSlope <- anota2seqSigObj@apvStats[[1]][useIds[i],
                                                            "apvSlope", drop = FALSE]
                  
                  tmpXmin <- min(dataX[useIds, , drop = FALSE])
                  tmpXmax <- max(dataX[useIds, , drop = FALSE])
                  tmpYmin <- min(dataY[useIds, , drop = FALSE])
                  tmpYmax <- max(dataY[useIds, , drop = FALSE])
                  mainTitle <- paste(useIds[i], "Slope:", round(tmpSlope,
                                                                digits = 2))
                  if (useGeneNames == TRUE) {
                      mainTitle <- paste(useIds[i], geneNames[useIds[i]],
                                         "Slope:", round(tmpSlope, digits = 2))
                  }
                  plot(x = c(tmpXmin, tmpXmax), y = c(tmpYmin,
                                                      tmpYmax), pch = "", main = mainTitle, xlab = labx,
                       ylab = laby)
                  phenoLev <- levels(as.factor(phenoVec))
                  for (j in 1:length(phenoLev)) {
                      text(x = dataX[useIds[i], phenoVec == phenoLev[j],
                                     drop = FALSE], y = dataY[useIds[i], phenoVec ==
                                                                  phenoLev[j], drop = FALSE], labels = phenoVec[phenoVec ==
                                                                                                                    phenoLev[j]], col = j)
                      tmpX <- mean(dataX[useIds[i], phenoVec == phenoLev[j],
                                         drop = FALSE])
                      tmpY <- mean(dataY[useIds[i], phenoVec == phenoLev[j],
                                         drop = FALSE])
                      tmpInt <- tmpY - (tmpSlope * tmpX)
                      lines(x = c(tmpXmin, tmpXmax), y = c((tmpInt +
                                                                tmpSlope * tmpXmin), (tmpInt + tmpSlope * tmpXmax)),
                            col = j)
                  }
                  deltaLine <- 1
                  deltaLine2 <- 0.5
                  xShift <- 4
                  plot(y = c(0, 10), x = c(0, 10), main = paste(useIds[i],
                                                                "APV statistics without RVM"), pch = "", xaxt = "n",
                       yaxt = "n", xlab = "", ylab = "")
                  nCont <- dim(anota2seqSigObj@usedContrasts)[2]
                  lineCount <- 10
                  xPos <- 2
                  
                  for (j in 1:nCont) {
                      sampCl1 <- rownames(anota2seqSigObj@usedContrasts)[anota2seqSigObj@usedContrasts[,
                                                                                                       j] < 0]
                      sampCl2 <- rownames(anota2seqSigObj@usedContrasts)[anota2seqSigObj@usedContrasts[,
                                                                                                       j] > 0]
                      tmpEff <- anota2seqSigObj@apvStats[[j]][useIds[i],
                                                              "apvEff",drop=FALSE]
                      tmpP <- anota2seqSigObj@apvStats[[j]][useIds[i],
                                                            "apvP",drop=FALSE]
                      tmpPadj <- anota2seqSigObj@apvStats[[j]][useIds[i],
                                                               "apvPAdj",drop=FALSE]
                      text(y = lineCount, x = xPos, labels = paste("Contrast:",
                                                                   j), font = 2, cex = 1.2)
                      lineCount <- lineCount - deltaLine2
                      text(y = lineCount, x = xPos, labels = paste("Effect:",
                                                                   round(tmpEff, digits = 2)))
                      lineCount <- lineCount - deltaLine2
                      text(y = lineCount, x = xPos, labels = paste("p-value:",
                                                                   round(tmpP, digits = 4)))
                      lineCount <- lineCount - deltaLine2
                      text(y = lineCount, x = xPos, labels = paste("adjusted p-value:",
                                                                   round(tmpPadj, digits = 3)))
                      lineCount <- lineCount - deltaLine
                      if (lineCount < 3) {
                          lineCount = 10
                          xPos = xPos + xShift
                      }
                  }
                  plot(y = c(0, 10), x = c(0, 10), main = paste(useIds[i],
                                                                "APV statistics with RVM"), pch = "", xaxt = "n",
                       yaxt = "n", xlab = "", ylab = "")
                  if (useRVM == TRUE) {
                      nCont <- dim(anota2seqSigObj@usedContrasts)[2]
                      lineCount <- 10
                      xPos <- 2
                      for (j in 1:nCont) {
                          sampCl1 <- rownames(anota2seqSigObj@usedContrasts)[anota2seqSigObj@usedContrasts[,
                                                                                                           j] < 0,drop=FALSE]
                          sampCl2 <- rownames(anota2seqSigObj@usedContrasts)[anota2seqSigObj@usedContrasts[,
                                                                                                           j] > 0,drop=FALSE]
                          tmpEff <- anota2seqSigObj@apvStatsRvm[[j]][useIds[i],
                                                                     "apvEff",drop=FALSE]
                          tmpP <- anota2seqSigObj@apvStatsRvm[[j]][useIds[i],
                                                                   "apvRvmP",drop=FALSE]
                          tmpPadj <- anota2seqSigObj@apvStatsRvm[[j]][useIds[i],
                                                                      "apvRvmPAdj",drop=FALSE]
                          text(y = lineCount, x = xPos, labels = paste("Contrast:",
                                                                       j), font = 2, cex = 1.2)
                          lineCount <- lineCount - deltaLine2
                          text(y = lineCount, x = xPos, labels = paste("Effect:",
                                                                       round(tmpEff, digits = 2)))
                          lineCount <- lineCount - deltaLine2
                          text(y = lineCount, x = xPos, labels = paste("p-value:",
                                                                       round(tmpP, digits = 4)))
                          lineCount <- lineCount - deltaLine2
                          text(y = lineCount, x = xPos, labels = paste("adjusted p-value:",
                                                                       round(tmpPadj, digits = 3)))
                          lineCount <- lineCount - deltaLine
                          if (lineCount < 3) {
                              lineCount = 10
                              xPos = xPos + xShift
                          }
                      }
                  }
              }
              if(plotToFile == TRUE){
                  dev.off()
              }
              
          })

setMethod("anota2seqGetOutputClass","Anota2seqDataSet",
          function(object , analysis, output) {
              
              if(!analysis %in% c("translated mRNA","total mRNA","translation","buffering","mRNA abundance")){
                  stop("analysis parameter wrong ... must be one of the following\n translated mRNA, total mRNA, translation or buffering ... ")
              }
              if(!output %in% c("full","selected")){
                  stop("output parameter wrong ... must be either full or selected")
              }
              
              if(analysis == "translated mRNA"){
                  if(is.null(object@translatedmRNA) == FALSE & output == "full") {
                      return(object@translatedmRNA)
                  } 
                  if(is.null(object@translatedmRNA) == TRUE & output == "full") {
                      return(NULL)
                  }
                  
                  if(is.null(object@selectedTranslatedmRNA) == FALSE & output == "selected"){
                      
                      return(object@selectedTranslatedmRNA)
                  } 
                  if(is.null(object@selectedTranslatedmRNA) == TRUE & output == "selected")
                  {
                      return(NULL)
                  }
              }
              
              if(analysis == "total mRNA"){
                  if(is.null(object@totalmRNA) == FALSE & output == "full"){
                      return(object@totalmRNA)
                  } 
                  if(is.null(object@totalmRNA) == TRUE & output == "full") {
                      return(NULL)
                  }
                  if(is.null(object@selectedTotalmRNA) == FALSE & output == "selected"){
                      return(object@selectedTotalmRNA)
                  } 
                  if(is.null(object@selectedTotalmRNA) == TRUE & output == "selected"){
                      return(NULL)
                  }
              }
              
              if(analysis == "translation"){
                  if(is.null(object@translation) == FALSE & output == "full"){
                      return(object@translation)
                  } 
                  if(is.null(object@translation) == TRUE & output == "full"){
                      return(NULL)
                  }
                  if(is.null(object@selectedTranslation) == FALSE & output == "selected"){
                      return(object@selectedTranslation)
                  } 
                  if(is.null(object@selectedTranslation) == TRUE & output == "selected") {
                      return(NULL)
                  }
              }
              
              if(analysis == "buffering"){
                  if(is.null(object@buffering) == FALSE & output == "full"){
                      return(object@buffering)
                  } 
                  if(is.null(object@buffering) == TRUE & output == "full"){
                      return(NULL)
                  }
                  if(is.null(object@selectedBuffering) == FALSE & output == "selected"){
                      return(object@selectedBuffering)
                  } 
                  if(is.null(object@selectedBuffering) == TRUE & output == "selected") {
                      return(NULL)
                  }
              }
              if(analysis == "mRNA abundance"){
                  if(is.null(object@mRNAAbundance) == FALSE){
                      return(object@mRNAAbundance)
                  }
                  if(is.null(object@mRNAAbundance) == TRUE){
                      return(NULL)
                  }
                  
              }
          })
setMethod("anota2seqSetOutput", "Anota2seqDataSet",
          function(object,analysis,output,input){
              
              if(!analysis %in% c("translated mRNA","total mRNA","translation","buffering")){
                  stop("analysis parameter wrong ... must be one of the following\n translated mRNA, total mRNA, translation or buffering ... ")
              }
              if(!output %in% c("full","selected")){
                  stop("output parameter wrong ... must be either full or selected")
              }
              
              if(analysis == "translated mRNA"){
                  if(output == "full") {
                      object@translatedmRNA <- input
                  } 
                  if(output == "selected"){
                      
                      object@selectedTranslatedmRNA <- input
                  }
              }
              
              if(analysis == "total mRNA"){
                  if(output == "full"){
                      object@totalmRNA <- input
                  }
                  if(output == "selected"){
                      object@selectedTotalmRNA <- input
                  }
              }
              
              if(analysis == "translation"){
                  if(output == "full"){
                      object@translation <- input
                  }
                  if(output == "selected"){
                      object@selectedTranslation <- input
                  }
              }
              
              if(analysis == "buffering"){
                  if(output == "full"){
                      object@buffering <- input
                  }
                  if(output == "selected"){
                      object@selectedBuffering <- input
                  }
              }
              return(object)
          })
setMethod("anota2seqSetSelectedOutput","Anota2seqDataSet",
          function(object,analysis,selContrast, input){
              if(!analysis %in% c("translated mRNA","total mRNA","translation","buffering")){
                  stop("analysis parameter wrong ... must be one of the following\n translated mRNA, total mRNA, translation or buffering ... ")
              }
              if(analysis == "translated mRNA"){
                  
                  object@selectedTranslatedmRNA@selectedData[[selContrast]] <-input[["selectedData"]]
                  object@selectedTranslatedmRNA@selectedRvmData[[selContrast]] <-input[["selectedRvmData"]]
                  #object@selectedPolysomeassociatedmRNA@groupIntercepts[[contrast]] <-input[["groupIntercepts"]]
                  object@selectedTranslatedmRNA@deltaData[[selContrast]] <-input[["deltaData"]]
                  object@selectedTranslatedmRNA@usedThresholds[[selContrast]] <-input[["usedThresholds"]]
                  object@selectedTranslatedmRNA@regModes[[selContrast]] <-input[["regModes"]]
                  
                  
              }
              if(analysis == "total mRNA"){
                  object@selectedTotalmRNA@selectedData[[selContrast]] <-input[["selectedData"]]
                  object@selectedTotalmRNA@selectedRvmData[[selContrast]] <-input[["selectedRvmData"]]
                  #object@selectedTotalmRNA@groupIntercepts[[contrast]] <-input[["groupIntercepts"]]
                  object@selectedTotalmRNA@deltaData[[selContrast]] <-input[["deltaData"]]
                  object@selectedTotalmRNA@usedThresholds[[selContrast]] <-input[["usedThresholds"]]
                  object@selectedTotalmRNA@regModes[[selContrast]] <-input[["regModes"]]
              }
              if(analysis == "translation"){
                  
                  object@selectedTranslation@selectedData[[selContrast]] <-input[["selectedData"]]
                  object@selectedTranslation@selectedRvmData[[selContrast]] <-input[["selectedRvmData"]]
                  #object@selectedTranslation@groupIntercepts[[contrast]] <-input[["groupIntercepts"]]
                  object@selectedTranslation@deltaData[[selContrast]] <-input[["deltaData"]]
                  object@selectedTranslation@usedThresholds[[selContrast]] <-input[["usedThresholds"]]
                  object@selectedTranslation@regModes[[selContrast]] <-input[["regModes"]]
                  
              }
              if(analysis == "buffering"){
                  object@selectedBuffering@selectedData[[selContrast]] <-input[["selectedData"]]
                  object@selectedBuffering@selectedRvmData[[selContrast]] <-input[["selectedRvmData"]]
                  #object@selectedBuffering@groupIntercepts[[contrast]] <-input[["groupIntercepts"]]
                  object@selectedBuffering@deltaData[[selContrast]] <-input[["deltaData"]]
                  object@selectedBuffering@usedThresholds[[selContrast]] <-input[["usedThresholds"]]
                  object@selectedBuffering@regModes[[selContrast]] <-input[["regModes"]]
                  
              }
              return(object)
          })
setMethod("anota2seqGetAvailableAnalyzes","Anota2seqDataSet",
          function(object){
              availableAnalyzes <- c("translated mRNA", "total mRNA", "translation", "buffering")[
                  c(!is.null(object@translatedmRNA), !is.null(object@totalmRNA),
                    !is.null(object@translation), !is.null(object@buffering))]
              if(length(availableAnalyzes) <= 0){
                  availableAnalyzes <- NULL
              }
              return(availableAnalyzes)
              
          })