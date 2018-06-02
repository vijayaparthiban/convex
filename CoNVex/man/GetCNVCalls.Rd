\name{GetCNVCalls}
\alias{GetCNVCalls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns CNV calls on all samples }
\description{
Returns CNV calls on all samples. Takes CNVfiles (file names -- one for each sample) as input and returns a data frame.
}
\usage{
GetCNVCalls(CNVfiles, names=1)
}
\arguments{
  \item{CNVfiles}{ CNV file names contains deletions and duplications -- one for each sample }
  \item{names}{ (default=1) Assigns column names; no names otherwise; }
}
\details{
GetCNVCalls
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  # Data frame Columns: "chr","start","end","num_probes","convex_score","cnv_type","sample_id"
}
\seealso{ 
SWCNV, SWCNVCall, SWCNVCallCommands
}
\examples{
> 	CallsAll = GetCNVCalls(CNVfiles);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GetCNVCalls }
