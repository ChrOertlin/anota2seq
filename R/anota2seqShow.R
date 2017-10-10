## show method for Anota2seqDataSet
## This method should give you a general overview of your data

setMethod("show","Anota2seqDataSet",
          definition = function(object=NULL){
              cat("Information on input data.\n")
              cat(paste("dataP and dataT contain data for",nrow(object@dataP),"mRNA and",ncol(object@dataP),"samples.","\n",sep=" "))
              cat("\n")
              cat("Your sample classes are:\n","\t",levels(as.factor(object@phenoVec)),"\n",sep=" ")
              cat("\n")
              if(is.null(object@contrasts) == FALSE){
                  cat("contrasts used in the analysis are:\n")
                  cat("\t\t\t",c(paste(colnames(object@contrasts),sep="\t"),"\n"))
                  for(rows in 1:dim(object@contrasts)[1]){
                      cat(c("\t",rownames(object@contrasts)[rows],"\t"))
                      for(cols in 1:dim(object@contrasts)[2]){
                          cat(c("\t",object@contrasts[rows,cols],"\t"))
                      }
                      cat("\n")
                  }
                  
                  cat("\n")
              }
              if(is.null(object@contrasts) == TRUE){
                  cat("\nNo custom contrasts set, default will be used.\n")
              }
              if(is.null(object@translatedmRNA) == TRUE){
                  cat("No differential expression analysis for translated mRNA has been detected.\n")
              }
              if(is.null(object@totalmRNA) == TRUE){
                  cat("No differential expression analysis for total mRNA has been detected.\n")
              }
              if(is.null(object@translation) == TRUE){
                  cat("No analysis of translation has been detected.\n")
              }
              if(is.null(object@buffering) == TRUE){
                  cat("No analysis of buffering has been detected.\n")
              }
              if(is.null(object@translatedmRNA) == FALSE){
                  cat("\nOutput for differential expression analysis of translated mRNA detected.\n")
              }
              if(is.null(object@totalmRNA) == FALSE){
                  cat("\nOutput for differential expression analysis of total mRNA detected.\n")
              }
              if(is.null(object@translation) == FALSE){
                  cat("\nOutput for analysis of translation detected.\n")
              }
              if(is.null(object@buffering) == FALSE){
                  cat("\nOutput for analysis of buffering detected.\n")
              }
              
              if(is.null(object@selectedTranslatedmRNA) == TRUE){
                  cat("No selected data of differential expression analysis for translated mRNA has been detected.\n")
              }
              if(is.null(object@selectedTotalmRNA) == TRUE){
                  cat("No selected data of differential expression analysis for total mRNA has been detected.\n")
              }
              if(is.null(object@selectedTranslation) == TRUE){
                  cat("No selected data of analysis of translation has been detected.\n")
              }
              if(is.null(object@selectedBuffering) == TRUE){
                  cat("No selected data of analysis of buffering has been detected.\n")
              }
              
              if(is.null(object@selectedTranslatedmRNA) == FALSE){
                  cat("\nOutput for differential expression analysis of translated mRNA detected:\n")
                  showSelectedOutput(object,"translated mRNA")
              }
              if(is.null(object@selectedTotalmRNA) == FALSE){
                  cat("\nOutput for differential expression analysis of total mRNA detected:\n")
                  showSelectedOutput(object,"total mRNA")
              }
              if(is.null(object@selectedTranslation) == FALSE){
                  cat("\nOutput for analysis of translation detected:\n")
                  showSelectedOutput(object,"translation")
              }
              if(is.null(object@selectedBuffering) == FALSE){
                  cat("\nOutput for analysis of buffering detected:\n")
                  showSelectedOutput(object,"buffering")
              }
              if(is.null(object@mRNAAbundance) == TRUE){
                  cat("\nNo output for mRNA abundance genes.\n")
              }
              if(is.null(object@mRNAAbundance) == FALSE){
                  cat("\nOutput for mRNA abundance genes detected.\n")
                  showSelectedOutput(object,"mRNA abundance")
              }
              if(is.null(object@mRNAAbundance) == FALSE){
                  cat("\nRegulatory modes selection:\n")
                  showRegModeOutput(object,"translation","translation")
                  showRegModeOutput(object,"abundance","mRNA abundance")
                  showRegModeOutput(object,"buffering","buffering")
                  
              }   
              
          })



