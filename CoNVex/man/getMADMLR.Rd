\name{getMADMLR}
\alias{getMADMLR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ MAD of mean log2 ratio or error-weighted score of a CNV region  }
\description{
After CNV detection, MAD is calculated from mean log2 ratio (across probe regions) and mean error-weighted scores of a CNV region. This is used to distinguish rare and common CNVs and potential false positives.
}
\usage{
getMADMLR(Fd_rep,aX,aXscore)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Fd_rep}{ Data frame containing CNV regions with the following columns: chr, start, end, num_probes, convex_score, cnv_type, sample_id  }
  \item{aX}{ Data frame containing log2 ratio: rowcount, chr, start, end, GC, dG, Tm, [log2_ratio of all samples - one sample per column] }
  \item{aXscore}{ Data frame containing MAD-weighted scores: rowcount, chr, start, end, GC, dG, Tm, [MAD-weighted scores of all samples - one sample per column].  }
}

\details{
getMADMLR returns following columns: mean_ratio_sample, mad_mean_ratio, mad_mean_admscore for each CNV.
}
\references{ R Library path }
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{ MAD-weighted scores are also called ADM scores or error-weighted scores in some places
}
\seealso{ 
SWCNV
}
\examples{
> convex_path = getMADMLR()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ getMADMLR, sw_lin64, sw_osx }
