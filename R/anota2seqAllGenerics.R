# accessor normalizedData
setGeneric("anota2seqGetNormalizedData",
           function(object=NULL) standardGeneric("anota2seqGetNormalizedData"))
# accessor phenoVec
setGeneric("anota2seqGetCovariates",
           function(object=NULL) standardGeneric("anota2seqGetCovariates"))
# accessor selContrast
setGeneric("anota2seqGetContrasts",
           function(object=NULL) standardGeneric("anota2seqGetContrasts"))
setGeneric("anota2seqGetOutput",
           function(object=NULL, analysis=NULL, output=NULL,selContrast=NULL,getRVM =NULL) standardGeneric("anota2seqGetOutput"))
setGeneric("anota2seqGetQualityControl",
           function(object=NULL) standardGeneric("anota2seqGetQualityControl"))
setGeneric("anota2seqGetResidOutlierTest",
           function(object=NULL) standardGeneric("anota2seqGetResidOutlierTest"))
setGeneric("anota2seqGetDeltaData",
           function(object=NULL,output=NULL,analysis=NULL,selContrast=NULL) standardGeneric("anota2seqGetDeltaData"))
setGeneric("anota2seqGetThresholds",
           function(object=NULL,analysis=NULL,selContrast=NULL) standardGeneric("anota2seqGetThresholds"))
setGeneric("anota2seqGetAvailableAnalyzes",
           function(object = NULL) standardGeneric("anota2seqGetAvailableAnalyzes"))

setGeneric("anota2seqPlotFC",
           function(object=NULL,visualizeRegModes="all",selContrast= NULL,fileName= NULL,plotToFile = TRUE, ...) standardGeneric("anota2seqPlotFC"))
setGeneric("anota2seqPlotPvalues",
           function(object=NULL,useRVM = TRUE,selContrast = NULL,myBw = 0.05,plotToFile=TRUE, fileName= NULL, ...) standardGeneric("anota2seqPlotPvalues"))

setGeneric("anota2seqPlotGenes",
           function(object= NULL,selContrast=NULL,analysis=NULL,geneNames = NULL,plotToFile = TRUE,fileName=NULL) standardGeneric("anota2seqPlotGenes"))

setGeneric("anota2seqSetSelectedOutput", 
           function(object=NULL,analysis=NULL,selContrast=NULL,input=NULL) standardGeneric("anota2seqSetSelectedOutput"))

setGeneric("anota2seqGetOutputClass",
           function(object=NULL, analysis=NULL, output=NULL) standardGeneric("anota2seqGetOutputClass"))
setGeneric("anota2seqSetOutput",
           function(object=NULL,analysis=NULL,output=NULL,input=NULL) standardGeneric("anota2seqSetOutput"))



