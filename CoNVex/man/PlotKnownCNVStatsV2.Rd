\name{PlotKnownCNVStatsV2}
\alias{PlotKnownCNVStatsV2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Number of Calls vs. Known (percentage) CNV Calls  }
\description{
Plots the #Calls vs. Known (percentage) for all samples and saves a call stats file
}
\usage{
PlotKnownCNVStatsV2(CallsAll,known_dels=data.frame(),known_dups=data.frame(),mark_outliers=0,outlier_from_median=5,file_prefix="CNVstats_NumCalls_Known")
}
\arguments{
  \item{CallsAll}{ CNV calls in all samples (data frame) - output of GetAllCalls }
  \item{known_dels}{ Known deletions in a data frame (default: uses the list supplied by CoNVex, if this option is not used) }
  \item{known_dups}{ Known duplications in a data frame (default: uses the list supplied by CoNVex, if this option is not used) }
  \item{mark_outliers}{ default=0 (don't marks outliers). =1 marks outliers; check the next option }
  \item{outlier_from_median}{ Marks outliers that are <outlier_from_median>\% away from the median. Default: N=5; (ignored if mark_outliers = 0)  }
  \item{file_prefix}{ Prefix for the file name (default: CNVstats_NumCalls_Known}
}
\details{
PlotKnownCNVStatsV2 uses the GetOverlap() function internally. This function calls the java program that gets the overlap. By default, temporary files re created
in /tmp folder before estimating the overlap. You can change this option in GetOverlap() function.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  # Output contains 1 plot in the current directory. You may like to use setwd("your_folder") to set the current directory, and an explicit
file name (should you avoid the default!)
}
\seealso{ 
PlotKnownCNVStatsV2, GetOverlap, getClassPath
}
\examples{
> setwd("~/")
> calls = GetCNVCalls(CNVfiles);
> PlotKnownCNVStatsV2(CallsAll=calls);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PlotKnownCNVStatsV2 }
