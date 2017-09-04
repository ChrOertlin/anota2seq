## anota2seq S4 class implementation
setClassUnion("listOrNULL",members=c("list", "NULL"))
setClassUnion("numOrNULL",members =c("numeric","NULL"))
setClassUnion("batchVec",members=c("numeric","character","NULL"))
setClassUnion("matrixOrNULL",members = c("matrix","NULL"))
setClass("anota2seqQualityControl",
         slots = c(
             omniIntStats = "matrix",
             omniGroupStats = "matrix",
             groupIntercepts = "matrix",
             correctionMethod = "character",
             dsfSummary = "numeric",
             dfbetas = "matrix",
             residuals = "matrix",
             fittedValues = "matrix",
             phenoClasses = "character",
             sampleNames = "character",
             abParametersInt = "numOrNULL",
             abParametersGroup = "numOrNULL"
         ))
setClass("anota2seqResidOutlierTest",
         slots = c(
             confInt = "numeric",
             inputResiduals = "matrix",
             rnormIter = "numeric",
             outlierMatrixLog = "matrix",
             meanOutlierPerIteration = "numeric",
             obtainedComparedToExpected = "numeric",
             nExpected = "numeric",
             nObtained= "numeric"
         ))
## Class that holds the standard output of the anota2seqAnalyse function (slightly reduced to avoid redundancy)
setClass("anota2seqOutput",
         slots=c(
             apvStats = "listOrNULL",
             apvStatsRvm = "listOrNULL",
             correctionMethod = "character",
             usedContrasts = "matrix",  
             abList = "list"
         ))
## Class that holds the output of anota2seqSelectSignificantGenes 
setClass("anota2seqSelectedOutput",
         slots= c(
             selectedData = "listOrNULL",
             selectedRvmData = "listOrNULL",
             useRVM = "logical",
             deltaData = "listOrNULL",
             usedThresholds = "listOrNULL",
             regModes = "logical"
         ))
setClass("mRNAabundanceOutput",
         slots = c(
             totalmRNA = "listOrNULL",
             translatedmRNA = "listOrNULL",
             useRVM = "logical",
             mRNASelect = "logical",
             regModes = "logical"
         ))
setClassUnion("qcOrNULL",members = c("anota2seqQualityControl", "NULL"))
setClassUnion("residOrNULL", members = c("anota2seqResidOutlierTest", "NULL"))
setClassUnion("outOrNULL", members = c("anota2seqOutput", "NULL"))
setClassUnion("selOutOrNULL", members = c("anota2seqSelectedOutput", "NULL"))
setClassUnion("abundOrNULL", members =  c("mRNAabundanceOutput", "NULL"))
setClass("anota2seqDataSet",
                             slots = c(
                                 dataT = "matrix",
                                 dataP = "matrix",
                                 phenoVec = "character",
                                 batchVec = "batchVec",
                                 contrasts = "matrixOrNULL",
                                 qualityControl = "qcOrNULL",
                                 residOutlierTest = "residOrNULL",
                                 translatedmRNA = "outOrNULL",
                                 totalmRNA = "outOrNULL",
                                 translation = "outOrNULL",
                                 buffering = "outOrNULL",
                                 selectedTranslatedmRNA="selOutOrNULL",
                                 selectedTotalmRNA = "selOutOrNULL",
                                 selectedTranslation = "selOutOrNULL",
                                 selectedBuffering = "selOutOrNULL",
                                 mRNAAbundance = "abundOrNULL",
                                 deltaData = "listOrNULL"
                             ))

