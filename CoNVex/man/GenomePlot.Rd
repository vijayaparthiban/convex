\name{GenomePlot}
\alias{GenomePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Creates a genome-wide plot of all exome bait regions in the GAMfile  }
\description{
Creates a genome-wide plot of all exome bait regions' ADM scores in the GAMfile. Input requires a single sample's GAMfile and its sample_id.
}
\usage{
GenomePlot(GAMfile,sample_id="",output_folder)
}
\arguments{
  \item{GAMfile}{ GAM file name of a sample - output from GAMCorrection }
  \item{sample_id}{ sample_id; default: '' - GAMfile used in place for plotting }
  \item{output_folder}{ Output folder for storing plots}
}
\details{
GenomePlot
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  # Output contains a genomeplot in the output folder. Use this function within a for() loop for multiple samples
}
\seealso{ 
GetCNVCalls
}
\examples{
>#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GenomePlot }
