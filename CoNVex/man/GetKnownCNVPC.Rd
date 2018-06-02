\name{GetKnownCNVPC}
\alias{GetKnownCNVPC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns known CNVs (percentage) in a data frame  }
\description{
For a given data frame of CNV calls, this function checks whether each CNV is known in the consensus CNV list bundled with
the CoNVex package. If so, it calculates the overlap (forward overlap - percentage of the length of the call present in the Consensus CNV list).
It returns a data frame with all columns in the input CallsAll data frame together with a percentage column at the end.
}
\usage{
GetKnownCNVPC(CallsAll,known_dels=data.frame(),known_dups=data.frame())
}
\arguments{
  \item{CallsAll}{ CNV calls in all samples - output of GetAllCalls }
  \item{known_dels}{ Known deletions, if you like to use a different list of known CNVs }
  \item{known_dups}{ Known duplications, if you like to use a different list of known CNVs }
}
\details{
GetKnownCNVPC uses GetOverlap() internally. The java binaries (jar files bundled with this R package) must be in classpath to execute this function.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  Please include CoNVex.jar and other jars in classpath if it is not configured previously (you might have already done this for executing 'java ReadDepth' program).
}
\seealso{ 
GetCNVCalls, GetOverlap, getClassPath
}
\examples{
> setwd("~/")
> calls = GetCNVCalls(CNVfiles);
> GetKnownCNVPC(CallsAll=calls)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GetKnownCNVPC, Known CNVs }
