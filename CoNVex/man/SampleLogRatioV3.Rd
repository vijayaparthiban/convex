\name{SampleLogRatioV3}
\alias{SampleLogRatioV3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ SampleLogRatioV3 - Optimised version of SampleLogRatioV2 using Unix's cut command (depends *nix style environment)  }
\description{
SampleLogRatioV3 (version=3) is similar to SampleLogRatioV2. i.e. It calculates log2 ratio using a median reference estimated from within a highly correlated subset of samples 
(as opposed to all samples in the analysis batch - like SampleLogRatio). This is beneficial if you use a heterogenous mixture of samples 
in the analysis batch. e.g., using same exome library, but the upstream data analysis (assembly, filtering, etc.) differs. 
}
\usage{
SampleLogRatioV3(RDfiles, sample_ids, sample_gender, regions_file, chrX=0, min_samples=25, RPKM=0, cor_matrix_file="", sample_means_file="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RDfiles}{ Vector of output read depth files }
  \item{sample_ids}{ Vector of Sample IDs }
  \item{sample_gender}{ Gender of the samples - 'M', 'F' and 'U' are allowed }
  \item{regions_file}{ Probe regions file}
  \item{chrX}{ [optional] cX=0 (default) excludes Chr X, cX=1 includes it }
  \item{min_samples=25}{ Minimum #samples expected in the subset. If chrX=1, both #males and #females should independently have this many number of samples }
  \item{RPKM=0}{ [optional] RPKM=1 Use RPKM (aka FPKM) to correlate the samples. RPKM=0 (default) Use read depth instead. } 
  \item{cor_matrix_file}{ Default="" [optional] Correlation matrix file (correlation between samples) }
  \item{sample_means_file}{ Default="" (don't save sample means) File to save sample means - this is useful for sample QC later. }
}
\details{
SampleLogRatioV3 function calculates log2 ratio (of depth/sample-median) for each probe region. chrX=0 (default) excludes Chr X. chrX=1 includes it. 
If ChrX=1, all samples must have 'M' or 'F' (male or female) so that Chr X is compared within males and females separately.
While using chrX=0, you can use 'U' (unknown) for the sample_gender. ChrY is fully excluded as it has very few probe regions, 
and the number of regions are significantly decreasing in the newer versions of off-the-shelf exome libraries. SampleLogRatioV3 differs from SampleLogRatio
by estimating median reference only from a correlated subset of samples. This function depends on Unix-like environment and cut command that comes 
pre-installed in those systems (including OSX). 
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  SampleLogRatioCall batch script calls it for all samples (use 'version=3' to use SampleLogRatioV3). 
SampleLogRatioCallCommands generates Unix commands for batch execution. Don't use this function if you sample size is too small.
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
