# GO over this again later to make sure it is done the right way .... 
anota2seqDataSetFromSE<- function(
    se = NULL,
    assayNum = NULL,
    dataType = NULL,
    normalize = FALSE,
    transformation = "TMM-log2",
    filterZeroGenes = FALSE,
    varCutOff = NULL)
{
    if(is.null(se) == TRUE){
        stop("Please provide a SummarizedExperiment ...")
    }
    anota2seqCheckParameter(normalize,dataType,transformation,filterZeroGenes,varCutOff,inFunc="dataset")
    #if summarizedExperiment is supplied dataP,dataT must be NULL ...
    if(is.null(se) == FALSE){
        if(class(se) != "SummarizedExperiment"){
            stop("Parameter se must be an object of class SummerizedExperiment.")
        }
    }
    if(is.null(assayNum) & length(assays(se)) > 1){
        stop("More than 1 assay detected in se but assayNum parameter not specified.\nPlease specify which assay to use.")
    }
    if(length(assays(se)) == 1 & is.null(assayNum)){
        assayNum <- 1
        warning("assayNum parameter not specified and length of assays(se) == 1. Will take the first assay data as input.\n")
    }
    
    anot<- colData(se)
    
    if("RNA" %in% colnames(anot) == FALSE){
        stop("colData in the provided SummarizedExperoment has no RNA column. Please read the anota2seqDataSetFromSE function help about the RNA column in colData of your SummarizedExperiment. \n")
    }
    if(length(levels(as.factor(anot[,"RNA"]))) > 2 | length(levels(as.factor(anot[,"RNA"]))) < 2 ){
        stop("RNA column in colData must describe exactly two RNA sources. Please read the anota2seqDataSet function help about the RNA column in colData of your SummarizedExperiment.")
    }
    if(length(levels(as.factor(anot[,"RNA"]))) == 2){
        msg <- NULL
        if("P"%in%anot[,"RNA"] == FALSE){
            msg <- c(msg,"No P found in RNA column in colData.\n")
        }
        if("T"%in%anot[,"RNA"] == FALSE){
            msg <- c(msg,"No T found in RNA column in colData.\n")
        }
        if(is.null(msg) == FALSE){
            stop(paste(msg,"Please read the anota2seqDataSet function help about the RNA column in colData of your SummarizedExperiment.\n",sep=""))
        }
    }
    if("samplePairs" %in% colnames(anot) == FALSE){
        stop("colData provided has no samplePairs column.  Please read the anota2seqDataSet function about the samplePairs column in colData of your SummarizedExperiment. \n")
    }
    if(length(unique(anot[,"samplePairs"])) < ncol(assays(se)[[assayNum]])/2){
        stop("Too few sample pairs described in colData. Please read the anota2seqDataSet function about the samplePairs column in colData of your SummarizedExperiment.\n")
    }
    if(length(unique(anot[,"samplePairs"])) > ncol(assays(se)[[assayNum]])/2){
        stop("Too many sample pairs described in colData. Please read the anota2seqDataSet function about the samplePairs column in colData of your SummarizedExperiment.\n")
    }
    
    if("treatment" %in% colnames(anot) == FALSE){
        stop("colData provided has no treatment column. Please read the anota2seqDataSet function help about the treatment column in colData of your SummarizedExperiment. \n")
    }
    
    anotP <- anot[anot[,"RNA"] == "P",]
    anotT <- anot[anot[,"RNA"] == "T",]
    
    anotP <- anotP[order(anotP[,"samplePairs"]),]
    anotT <- anotT[order(anotT[,"samplePairs"]),]
    
    # get DataP and dataT
    dataP <- assays(se)[[assayNum]][,rownames(anotP)]
    dataT <- assays(se)[[assayNum]][,rownames(anotT)]
    # order dataP and dataT by samplePairs 
    dataP <- dataP[,rownames(anotP)]
    dataT <- dataT[,rownames(anotT)]
    # Get phenoVec - Treatment column of anotP/anotT should correspond to this
    phenoVec <- as.vector(anotP[,"treatment"])
    # get batchVec if present ...
    
    if("batch"%in%colnames(anotP) == TRUE){
        batchVec <- as.vector(anotP[,"batch"])
    }
    if("batch"%in%colnames(anotP) == FALSE){
        batchVec <- NULL
    }
    anota2seqCheckInput(dataP,
               dataT,
               phenoVec,
               batchVec,
               NULL,
               "BH",
               inFunc="fromSE")
    

    
    if(dataType == "RNAseq"){
        preProcess <- anota2seqRNAseqPreProcessing(dataP=dataP,
                                                   dataT=dataT,
                                                   transformation =transformation,
                                                   filterZeroGenes=filterZeroGenes,
                                                   normalize=normalize)
        
        dataT <- preProcess$dataT
        dataP <- preProcess$dataP
    }
    

    
    # Check for genes that show low variance in the dataset ...
    varCheck <- anota2seqFiltCheckVar(tmpdataP = dataP,
                                  tmpdataT = dataT,
                                  varCutOff = varCutOff,
                                  phenoVec = phenoVec)
    
    dataP <- varCheck$dataP
    dataT <- varCheck$dataT
    if(max(range(dataP)) > 100 | max(range(dataT)) > 100){
        message()
        stop("Input data range indicates a non continuous scale.\nMake sure the input data is normalized and if coming from RNAsequencing transformed to a continuous scale.\n")
    }
    # initialize the anota2seqDataSet class first so that checks on phenoVec and contrast get performed.
    anota2seqClass <- new("anota2seqDataSet",
                          dataP = dataP,
                          dataT = dataT,
                          phenoVec = phenoVec,
                          batchVec = batchVec,
                          contrasts = NULL,
                          qualityControl = NULL,
                          residOutlierTest = NULL,
                          translatedmRNA = NULL,
                          totalmRNA = NULL,
                          translation = NULL,
                          buffering = NULL,
                          selectedTranslatedmRNA = NULL,
                          selectedTotalmRNA = NULL,
                          selectedTranslation = NULL,
                          selectedBuffering = NULL,
                          mRNAAbundance = NULL,
                          deltaData = NULL)
    message("All input checkpoints passed.\n")
    return(anota2seqClass)
    
}
