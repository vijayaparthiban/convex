\name{PlotCNVStats}
\alias{PlotCNVStats}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Creates plots of  (1) number of deletions and duplications, and (2) ratio of deletions/duplications in each sample  }
\description{
Plots CNV calls stats grom a given data frame of CNV calls.
}
\usage{
PlotCNVStats(CNVcallsAll,file_prefix="CNVstats",outlier_cutoff=2)
}
\arguments{
  \item{CNVcallsAll}{ CNV calls in all samples - output of GetAllCalls }
  \item{file_prefix}{ file name prefix for plots; default: CNVstats }
  \item{outlier_cutoff}{ The y-axis ratio threshold marking outliers in 'deletions/duplications' plot; default=2}
}
\details{
PlotCNVStats
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  # Output contains 2 plots in the current directory. You may like to use setwd("your_folder") to set the current directory 
}
\seealso{ 
GetCNVCalls
}
\examples{
> setwd("~/")
> calls = GetCNVCalls(CNVfiles);
> PlotCNVStats(calls);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PlotCNVStats }
