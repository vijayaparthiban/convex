\name{GetGAMScores}
\alias{GetGAMScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Returns a list of all samples' log2 ratio and ADM scores }
\description{
Returns a list of all samples' log2 ratio and ADM scores. Returned list contains 2 matrices - one with log2 ratio, other with ADM scores
}
\usage{
GetGAMScores(GAMfiles)
}
\arguments{
  \item{GAMfiles}{ GAM file names - output files of GAMCorrection }
}
\details{
GetGAMScores
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  # Returned list contains 2 matrices - one with log2 ratio, other with ADM scores
}
\seealso{ 
GAMCorrection, GAMCorrectionCommands
}
\examples{
> 	L2R_ADM = GetGAMScores(GAMfiles);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GetGAMScores }
