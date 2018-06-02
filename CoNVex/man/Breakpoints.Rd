\name{Breakpoints}
\alias{Breakpoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate Breakpoints using Number of Reads }
\description{
MAD of the log2 ratio varies within a set of regions based on the number of reads that hit the probe region. Based on this, MADs are calculated. 
For that purpose, probe regions are binned according to the number of reads making sure 
that each bin has <max_bin_size> number of probe regions. Default (1000) works fine in most cases. 
This function returns the max of the total reads for each bin (breakpoints using cut function).
}
\usage{
Breakpoints(RDfile_RepSample,max_bin_size=1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RDfile_RepSample}{ Read depth file of a representative sample in the analysis batch (any sample should work)  }
  \item{max_bin_size}{ Bin size for number of probe regions from which MAD is estimated (default=1000 should be enough for most cases) }
}
\details{
Breakpoints.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  Differences between exome samples generated with same experimental protocols/parameters should not make 
a huge difference between samples. Each sample can use the same breakpoints generated
from a representative sample. However, you may like to select a representative sample using the Number of total mapped reads
per sample. i.e. Consider a sample that has the median number of total reads across samples!
}
\seealso{ 
BreakpointsCall, BreakpointsCallCommands, cut
}
\examples{
# BP = Breakpoints(RDfile_RepSample=RDfiles[1],max_bin_size=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ breakpoints, median absolute deviation }
