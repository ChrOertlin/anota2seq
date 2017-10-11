# accessor normalizedData
setGeneric("anota2seqGetNormalizedData",
           function(object) standardGeneric("anota2seqGetNormalizedData"))
# accessor phenoVec
setGeneric("anota2seqGetCovariates",
           function(object) standardGeneric("anota2seqGetCovariates"))
# accessor selContrast
setGeneric("anota2seqGetContrasts",
           function(object) standardGeneric("anota2seqGetContrasts"))
setGeneric("anota2seqGetOutput",
           function(object, analysis, output,selContrast,getRVM=TRUE) standardGeneric("anota2seqGetOutput"))
setGeneric("anota2seqGetQualityControl",
           function(object) standardGeneric("anota2seqGetQualityControl"))
setGeneric("anota2seqGetResidOutlierTest",
           function(object) standardGeneric("anota2seqGetResidOutlierTest"))
setGeneric("anota2seqGetDeltaData",
           function(object,output,analysis,selContrast) standardGeneric("anota2seqGetDeltaData"))
setGeneric("anota2seqGetThresholds",
           function(object,analysis,selContrast) standardGeneric("anota2seqGetThresholds"))
setGeneric("anota2seqGetAvailableAnalyzes",
           function(object) standardGeneric("anota2seqGetAvailableAnalyzes"))

setGeneric("anota2seqPlotFC",
           function(object,visualizeRegModes="all",selContrast,fileStem = "ANOTA2SEQ_FoldchangePlot",plotToFile = TRUE, ...) standardGeneric("anota2seqPlotFC"))
setGeneric("anota2seqPlotPvalues",
           function(object,useRVM = TRUE,selContrast,myBw = 0.05,plotToFile=TRUE, fileStem = "ANOTA2SEQ_pvalue_density", ...) standardGeneric("anota2seqPlotPvalues"))

setGeneric("anota2seqPlotGenes",
           function(object,selContrast,analysis,geneNames = NULL,plotToFile = TRUE,fileStem = "ANOTA2SEQ_significantGenes_plot") standardGeneric("anota2seqPlotGenes"))

setGeneric("anota2seqSetSelectedOutput", 
           function(object=NULL,analysis=NULL,selContrast=NULL,input=NULL) standardGeneric("anota2seqSetSelectedOutput"))

setGeneric("anota2seqGetOutputClass",
           function(object=NULL, analysis=NULL, output=NULL) standardGeneric("anota2seqGetOutputClass"))
setGeneric("anota2seqSetOutput",
           function(object=NULL,analysis=NULL,output=NULL,input=NULL) standardGeneric("anota2seqSetOutput"))



