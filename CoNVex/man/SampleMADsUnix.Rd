\name{SampleMADsUnix}
\alias{SampleMADsUnix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate median absolute deviation of log2 ratio  }
\description{
SampleMADsUnix calculates MAD (median absolute deviation) of the corrected log2 ratio for each sample. 
In addition, it also calculates SampleMedians for Autosomes and ChrX
These metrics are used sample/call QC as well as for plotting and visualisation.
}
\usage{
SampleMADsUnix(sample_ids, GAMfiles, regions_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_ids}{ Vector of Sample IDs }
  \item{GAMfiles}{ GAM file names - output files of GAMCorrection }
  \item{regions_file}{ Probe regions file}
}
\details{
SampleMADsUnix uses GAMfiles to get the corrected log2 ratio (4th column) and calculates MAD.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  This can be used together with SampleMeans() function that calculates mean depth to plot Sample Means vs. MADs plots. 
This would explain sample noise in the analysis batch. This function can be called using SampleMeansMads bundled in the Rbatch folder. 
(MAD = Median Absolute Deviation)
}
\seealso{ 
SampleMeans, SampleMADsUnixMads
}
\examples{
> SampleMADsUnix(sample_ids=sids, GAMfiles=GAMcorrectionFiles,regions_file="path/to/regions_file")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SampleMADsUnix, Mean Depth }
