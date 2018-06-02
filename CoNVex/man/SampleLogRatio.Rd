\name{SampleLogRatio}
\alias{SampleLogRatio}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ SampleLogRatio - log2 ratio for each sample  }
\description{
This calculates log2 ratio for each sample using a median reference estimated from all samples
}
\usage{
SampleLogRatio(RDfiles, sample_ids, sample_gender, regions_file, chrX=0, sample_means_file="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RDfiles}{ Vector of output read depth files }
  \item{sample_ids}{ Vector of Sample IDs }
  \item{sample_gender}{ Gender of the samples - 'M', 'F' and 'U' are allowed }
  \item{regions_file}{ Probe regions file}
  \item{chrX}{ [optional] cX=0 (default) excludes Chr X, cX=1 includes it }
  \item{sample_means_file}{ Default="" (don't save sample means) File to save sample means - this is useful for sample QC later. }
}
\details{
SampleLogRatio function calculates log2 ratio (of depth/sample-median) for each probe region. chrX=0 (default) excludes Chr X. chrX=1 includes it. 
If Chr X=1, all samples must have 'M' or 'F' (male or female) so that Chr X is compared within males and females separately.
While using chrX=0, you can use 'U' (unknown) for the sample_gender. ChrY is fully excluded as it has very few probe regions, 
and the number of regions are significantly decreasing in the newer versions of off-the-shelf exome libraries.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  SampleLogRatioCall batch script calls it for all samples. 
SampleLogRatioCallCommands generates Unix commands for batch execution.
}
\seealso{ 
SampleLogRatioCallCommands, SampleLogRatioCall
}
\examples{
# Check examples of SampleLogRatioCallCommands
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ log2 ratio, sample median reference }
