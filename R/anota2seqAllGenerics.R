# accessor normalizedData
setGeneric("anota2seq.get.normalizedData",
           function(object=NULL) standardGeneric("anota2seq.get.normalizedData"))
# accessor phenoVec
setGeneric("anota2seq.get.covariates",
           function(object=NULL) standardGeneric("anota2seq.get.covariates"))
# accessor selContrast
setGeneric("anota2seq.get.contrasts",
           function(object=NULL) standardGeneric("anota2seq.get.contrasts"))
setGeneric("anota2seq.get.output",
           function(object=NULL, analysis=NULL, output=NULL,selContrast=NULL,getRVM =NULL) standardGeneric("anota2seq.get.output"))
setGeneric("anota2seq.get.qualityControl",
           function(object=NULL) standardGeneric("anota2seq.get.qualityControl"))
setGeneric("anota2seq.get.residOutlierTest",
           function(object=NULL) standardGeneric("anota2seq.get.residOutlierTest"))
setGeneric("anota2seq.get.deltaData",
           function(object=NULL,output=NULL,analysis=NULL,selContrast=NULL) standardGeneric("anota2seq.get.deltaData"))
setGeneric("anota2seq.get.thresholds",
           function(object=NULL,analysis=NULL,selContrast=NULL) standardGeneric("anota2seq.get.thresholds"))
setGeneric("anota2seq.get.availableAnalyzes",
           function(object = NULL) standardGeneric("anota2seq.get.availableAnalyzes"))

setGeneric("anota2seqPlotFC",
           function(object=NULL,visualizeRegModes="all",selContrast= NULL,fileName= NULL,plotToFile = TRUE, ...) standardGeneric("anota2seqPlotFC"))
setGeneric("anota2seqPlotPvalues",
           function(object=NULL,useRVM = TRUE,selContrast = NULL,myBw = 0.05,plotToFile=TRUE, fileName= NULL, ...) standardGeneric("anota2seqPlotPvalues"))

setGeneric("anota2seqPlotGenes",
           function(object= NULL,selContrast=NULL,analysis=NULL,geneNames = NULL,plotToFile = TRUE,fileName=NULL) standardGeneric("anota2seqPlotGenes"))

setGeneric("set.selected.output", 
           function(object=NULL,analysis=NULL,selContrast=NULL,input=NULL) standardGeneric("set.selected.output"))

setGeneric("anota2seq.get.output.class",
           function(object=NULL, analysis=NULL, output=NULL) standardGeneric("anota2seq.get.output.class"))
setGeneric("set.output",
           function(object=NULL,analysis=NULL,output=NULL,input=NULL) standardGeneric("set.output"))



