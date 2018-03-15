\docType{methods}
\name{anota2seq-plots}
\alias{anota2seq-plots}
\alias{anota2seqPlotFC}
\alias{anota2seqPlotPvalues}
\alias{anota2seqPlotGenes}
\alias{anota2seqPlotFC,Anota2seqDataSet-method}
\alias{anota2seqPlotPvalues,Anota2seqDataSet-method}
\alias{anota2seqPlotGenes,Anota2seqDataSet-method}
\title{Visualization of anota2seq results}
\usage{
anota2seqPlotFC(object, visualizeRegModes = "all",
  selContrast, contrastName = NULL, fileStem = "ANOTA2SEQ_FoldchangePlot", plotToFile = TRUE, myYlim = NULL, myXlim = NULL, ...)

anota2seqPlotPvalues(object, useRVM = TRUE, selContrast, contrastName = NULL,
  myBw = 0.05, plotToFile = TRUE, fileStem = "ANOTA2SEQ_pvalue_density", ...)

anota2seqPlotGenes(object, selContrast, analysis, geneNames = NULL, 
plotToFile = TRUE, fileStem = "ANOTA2SEQ_significantGenes_plot")

\S4method{anota2seqPlotFC}{Anota2seqDataSet}(object,
  visualizeRegModes = "all", selContrast, contrastName = NULL, fileStem = "ANOTA2SEQ_FoldchangePlot",
  plotToFile = TRUE, myYlim = NULL, myXlim = NULL, ...)

\S4method{anota2seqPlotPvalues}{Anota2seqDataSet}(object, 
  useRVM = TRUE, selContrast, contrastName = NULL, myBw = 0.05, plotToFile = TRUE,
  fileStem = "ANOTA2SEQ_pvalue_density", ...)

\S4method{anota2seqPlotGenes}{Anota2seqDataSet}(object,
  selContrast, analysis, geneNames = NULL,
  plotToFile = TRUE, fileStem = "ANOTA2SEQ_significantGenes_plot")
}
\arguments{
\item{object}{An object of class Anota2seqDataSet. Should contain output of 
the anota2seqAnalyze function. To visualize significant identifiers in the 
fold change plot (\code{anota2seqPlotFC}) and for regression plots 
(\code{anota2seqPlotGenes}), the object should also contain the output of 
anota2seqSelSigGenes and/or anota2seqRegModes (see Details).}

\item{visualizeRegModes}{Can be set to "all" (default, genes which were 
selected by \code{\link{anota2seqSelSigGenes}} are colored according to 
their allocated regulatory mode (change in mRNA abundance or change in 
translational efficiency leading to altered protein levels or buffering), 
"none" (no colors), "translation" (only identifiers regulated by change in 
translational efficiency leading to altered protein levels are colored), 
"buffering" (only identifiers regulated by change in translational 
efficiency leading to buffering are colored). The object will be required 
to contain different outputs depending on the value given to this 
parameter (see Details).}

\item{selContrast}{Which contrast(s) should be considered? Descriptions of
the contrasts can be found in the output from the anota2seqAnalyze object
in the usedContrasts slot. Indicate the contrast(s) by a numeric vector of
the column number(s).}

\item{contrastName}{Custom name annotation for the selected contrast(s). Provide a character name for each selected contrast, this name will be used as plot title in the anota2seqPlotPvalues plots or as xlab and ylab annotation in the anota2seqPlotFC plots.}

\item{fileStem}{if plotToFile is TRUE, this stem will be added in front of the
output filename(s).}

\item{plotToFile}{Boolean. If set to TRUE (default) the function will output 
the plot as a PDF file. If set to FALSE plots will be plotted to the R 
default graphics device.}

\item{myXlim}{Specify the x-axis limits used for the scatterplots(s) in the anota2seqPlotFC function. Default is NULL.}

\item{myYlim}{Specify the y-axis limits used for the scattersplot(s) in the anota2seqPlotFC function. Default is NULL}

\item{...}{Graphical parameters that are passed to the par() function.}

\item{useRVM}{Should the density plot be performed on RVM p-values/FDR 
(default) or no-RVM p-values/FDR?}

\item{myBw}{Smoothing bandwidth parameter for used in the stats::density 
function. Within one plot, the same bandwidth will be used for all density 
curves.}

\item{analysis}{For which analysis should anota2seqPlotGenes be performed. Can be set 
to "translation" or "buffering".}

\item{geneNames}{When anota2seqPlotGenes performs the individual gene plots they will 
be named by the original row names supplied to the anota2seqAnalyze 
function. geneNames allows the user to add additional names when plotting 
to e.g. include gene symbols. Input is a matrix with one column where the 
original row names match the row names of the input matrix and the desired 
new names are given in column 1. Default is NULL i.e. no additional names.}
}
\value{
No value is returned. These functions generate a graphical outputs as
  described above.
}
\description{
These functions generate graphical outputs aiming to provide an 
  overview of the extents of different regulatory modes of gene expression. 
  The graphical outputs consist of: \enumerate{ \item For all selected 
  contrasts, P-value and FDR density plots indicating the distribution of 
  P-values and FDRs for all analyzes found in the Anota2seqDataSet object. 
  \item Fold change plots where total mRNA log2FC is compared to translated 
  mRNA (e.g. polysome-associated mRNA or RPF) log2FC (for all selected 
  contrasts). \item Per identifier and per treatment fitted regression curves
  between total and translated mRNA for all samples which helps evaluating 
  regulation of single identifiers. }
}
\details{
\code{anota2seqPlotFC}: if visualizeRegModes is set to "none", the 
  Anota2seqDataSet object is only required to contains an output of 
  anota2seqAnalyze. If visualizeRegModes is set to "translation" or 
  "buffering", the object should contain the output of anota2seqAnalyze and 
  anota2seqSelSigGenes for the corresponding analysis. Finally, if 
  visualizeRegModes is set to "all" (default), anota2seqRegModes should have 
  been called on the Anota2seqDataSet object so that identifiers can be 
  colored according to the three regulatory modes (i.e. change in mRNA 
  abundance or change in translational efficiency leading to altered protein 
  levels or buffering).
  
  \code{anota2seqPlotGenes}: requires an Anota2seqDataSet object containing the output
  of \code{\link{anota2seqAnalyze}} and \code{\link{anota2seqSelSigGenes}} for
  translation and/or buffering. In the graphical output of \code{anota2seqPlotGenes}, 
  the results for each significant gene is displayed on a separate row. The 
  first graph shows all samples and per treatment regression lines using the 
  common slope with different colors for each treatment. The magnitude of the
  common slope is indicated. The second graph shows key statistics for the 
  identifier without the RVM model for all contrasts analyzed when running 
  anota2seqAnalyze. The third graph is similar to the second but with RVM 
  statistics instead.
}
\examples{
\dontrun{
data(anota2seq_data)
# Initialize the Anota2seqDataSet
Anota2seqDataSet <- anota2seqDataSetFromMatrix(dataP = anota2seq_data_P[1:1000,],
                                      dataT = anota2seq_data_T[1:1000,],
                                      phenoVec = anota2seq_pheno_vec,
                                      dataType = "RNAseq",
                                      normalize = TRUE)
# Perform anota2seqRun function (performQC is set to FALSE here to limit running 
# time of this example but the model assumptions should be assessed (see help of
# anota2seqPerformQC))                               
Anota2seqDataSet <- anota2seqRun(Anota2seqDataSet, 
                                 performQC = FALSE,
                                 performROT = FALSE)
# Visualize results
anota2seqPlotPvalues(Anota2seqDataSet, plotToFile = FALSE, selContrast = 1)
anota2seqPlotFC(Anota2seqDataSet, plotToFile = FALSE, selContrast = 1)
anota2seqPlotGenes(Anota2seqDataSet, plotToFile = TRUE, analysis = "translation",
           selContrast = 1)
           }

}
