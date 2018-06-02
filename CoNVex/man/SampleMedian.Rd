\name{SampleMedian}
\alias{SampleMedian}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ SampleMedian - Median Depth across Samples  }
\description{
This function calculates median depth of each probe region across samples and returns them in two columns for males and females
separately. Depth of autosomal probe regions are calculated across all samples (so, these are exactly same in two columns for ease of use), 
while those of chromosome X is calculated within males and females separately. If ChrX=0, only autosomal probe regions' medians are returned.
Returned data frame includes the following columns: rowcount, chr, start, end, GC, dG, Tm, Median_Male, Median_Female
}
\usage{
SampleMedian(RDfiles, sample_ids, sample_gender, regions_file, chrX=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RDfiles}{ Vector of output read depth files }
  \item{sample_ids}{ Vector of Sample IDs }
  \item{sample_gender}{ Gender of the samples - 'M', 'F' and 'U' are allowed }
  \item{regions_file}{ Probe regions file}
  \item{chrX}{ [optional] cX=0 (default) excludes Chr X, cX=1 includes it }
}
\details{
SampleMedian function calculates median depth of each probe region across samples. chrX=0 (default) excludes Chr X. chrX=1 includes it. 
If Chr X=1, all samples must have 'M' or 'F' (male or female) so that Chr X is compared within males and females separately.
While using chrX=0, you can use 'U' (unknown) for the sample_gender. ChrY is fully excluded as it has very few probe regions, 
and the number of regions are significantly decreasing in the newer versions of off-the-shelf exome libraries.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  SampleMedian is not required for CNV detection, as it's also a part of SampleLogRatio internally.
Batch scripts are not available for this function.
}
\seealso{ 
SampleMedianCallCommands
SampleMedianCall
}
\examples{
# sample_median = SampleMedian(RDfiles, sample_ids, sample_gender, regions_file, chrX=0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ log2 ratio, sample median reference }
