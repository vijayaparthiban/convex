\name{SampleMeans}
\alias{SampleMeans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimate mean depth across probe regions  }
\description{
SampleMeans calculates mean depth across probe regions within each sample in the analysis batch. It returns a vector (equal to the length of 'sample_ids') of 
mean depths. This can be used for further plotting and visualisation. 
}
\usage{
SampleMeans(sample_ids, RDfiles, regions_file)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_ids}{ Vector of Sample IDs }
  \item{RDfiles}{ Vector of output read depth files }
  \item{regions_file}{ Probe regions file}
}
\details{
SampleMeans uses RDfiles (read depth files) generated using the java ReadDepth program for each sample. 
It simply calculates the mean of the read depths of probe regions within the sample.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  This can be used together with SampleMADsUnix() function that calculates MAD(log2 ratio) to plot Sample Means vs. MADs plots. 
This would explain sample noise in the analysis batch. This function can be called using SampleMeansMads bundled in the Rbatch folder. 
(MAD = Median Absolute Deviation)
}
\seealso{ 
SampleMADsUnix, SampleMeansMads
}
\examples{
> SampleMeans(sample_ids=sids, RDfiles=ReadDepthFiles,regions_file="path/to/regions_file")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SampleMeans, Mean Depth }
