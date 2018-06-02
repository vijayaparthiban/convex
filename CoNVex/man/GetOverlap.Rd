\name{GetOverlap}
\alias{GetOverlap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get overlap between two lists of chromosome regions }
\description{
This functions accepts two lists of chromosome regions. Each list should have the following columns (1-3): chr, start, end (chr must be 1,2,3,...,X,Y). It returns max percentage overlap
of each region in chr_list1 with any of the regions in chr_list2. i.e. chr_list2 is returned with an additional column at the end having the percentage.
}
\usage{
GetOverlap(chr_list1,chr_list2, tmp_folder="/tmp/")
}
\arguments{
  \item{chr_list1}{ Chr list 1: chr, start, end columns }
  \item{chr_list2}{ Chr list 2: chr, start, end columns }
  \item{tmp_folder}{ (Default: '/tmp/') (optional) Temporary files folder }
}
\details{
GetOverlap is internally used by other functions. The ConVex java binary (CoNVex.jar bundled with CoNVex R package) must be in classpath to execute this function.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  Please include CoNVex.jar in classpath if it is not configured previously (you might have already done this for executing 'java ReadDepth' program).
}
\seealso{ 
GetCNVCalls
}
\examples{
GetOverlap(chr_list1,chr_list2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GetOverlap, Chromosome Overlap }
