# GO over this again later to make sure it is done the right way .... 
anota2seqDataSetFromMatrix <- function(
    dataP = NULL,
    dataT = NULL,
    phenoVec = NULL,
    batchVec = NULL,
    dataType = NULL,
    normalize = FALSE,
    transformation = "TMM-log2",
    filterZeroGenes = FALSE,
    varCutOff = NULL)
{
    
    if(is.null(batchVec) == FALSE){
        batchVec <- as.character(batchVec)
    }
    anota2seqCheckParameter(normalize,dataType,transformation,filterZeroGenes,varCutOff,inFunc="dataset")
    anota2seqCheckInput(dataP,
                        dataT,
                        phenoVec,
                        batchVec,
                        NULL,
                        "BH",
                        inFunc="fromMatrix")
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
        stop("Input data range indicates a non continuous scale. \n
             Make sure the input data is normalized and if coming from RNAsequencing transformed to a continuous scale.\n")
    }
    # initialize the anota2seqDataSet class first so that checks on phenoVec and contrast get performed.
    anota2seqClass <- new("anota2seqDataSet",
                          dataP = dataP,
                          dataT = dataT,
                          phenoVec = as.character(phenoVec),
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




