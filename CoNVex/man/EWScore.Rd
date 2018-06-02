\name{EWScore}
\alias{EWScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate MAD-weighted Scores }
\description{
This function calculates error-weighted scores from corrected log2 ratio and regional MADs. 
The calculation is partly similar to Agilent's ADM2 algorithm, but the exome log2 ratio 
and MAD-based calculation of weights are unrelated to ADM2 which was developed for aCGH data.
Regional MADs are first calculated using number of reads and breakpoints. 
}
\usage{
EWScore(d,doNorm=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{ Data frame having the following columns: chr, start, end, log2 ratio (corrected), MAD (regional)  }
  \item{doNorm}{ doNorm=0 (default) does not do any normalization (recommended). Other values normalize (Z-score) the log2 ratio. 
  In future, this can be extended for better normalization if required. }
}
\details{
EWScore
}
\references{ Agilent G4175AA CGH Analytics 3.4 User Guide, ADM2 algorithm }
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  This function is internally called from GAMCorrection function. Changing parameters are recommended only 
for advanced users (mostly CoNVex authors)
}
\seealso{ 
GAMCorrection\cr
Agilent user guides: \url{http://www.chem.agilent.com/} \cr

}
\examples{
> ADMrep = EWScore(d, doNorm=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EWScore, median absolute deviation }
