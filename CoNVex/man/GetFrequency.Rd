\name{GetFrequency}
\alias{GetFrequency}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get frequency within or between Chr regions }
\description{
This function accepts upto two lists (as data frames) of chromosome regions as CNV calls. Each list should have chr, start, end (chr must be 1,2,3,...,X,Y) as the first 3 columns. 
Additionally, sample_id must be provided in the 7th column as in the CoNVex output of CNV calls. It returns the number of times each region (including itself) 
is encountered in the Chr (CNV) regions list. If a region in one sample overlaps multiple (split) regions in another sample, it's counted only once. So, the 
count is the number of samples (including the sample to which the CNV belongs) that have at least one overlapping CNV region. If you divide the
count with the number of total samples in the list, you will get a proportion. chr_list1 is returned with two additional columns at the end having the frequency as count and percentage. 
Refer to the details section below.
}
\usage{
GetFrequency(chr_list1=data.frame(), chr_list2=data.frame(), sid_col1=7, sid_col2=7, ro_threshold=0, tmp_folder="/tmp/")
}
\arguments{
  \item{chr_list1}{ (mandatory) List of Chr regions (e.g., CNV calls) }
  \item{chr_list2}{ (optional) List of Chr regions to compare against chr_list1 }
  \item{sid_col1}{ (optional) Sample ID column (sorted) in chr_list1; Default: column 7 }
  \item{sid_col2}{ (optional) Sample ID column (sorted) in chr_list2; Default: column 7 }
  \item{ro_threshold}{ (optional) Reciprocal overlap threshold in percentage; Default: > 0 (any overlap) }
  \item{tmp_folder}{ (Default: '/tmp/') (optional) Temporary files folder }
}
\details{
GetFrequency calculates the number of times each region (including itself) is encountered in the Chr regions list.  If chr_list1 is the only given list,  
GetFrequency() calculates the internal frequency within the list through ANY overlap between regions.  If chr_list2 is ALSO given, this is the number of
samples in chr_list2 in which each CNV/Chr region from chr_list1 is seen.

Possible use cases:
(1) CNV calls' internal frequency (within chr_list1 only with sample_id)
(2) CNV calls (chr_list1 with sample_id) and the number of times each call is seen in AN EXTERNAL BATCH of CNV calls (chr_list2 with sample_id)
This function simplifies and conveniently encapsulates the complexity of CNV frequency estimation. 
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  Please include CoNVex.jar in classpath if it is not configured previously (you might have already done this for executing 'java ReadDepth' program).
}
\seealso{ 
GetCNVCalls
}
\examples{
GetFrequency(chr_list1) # or
GetFrequency(chr_list1,chr_list2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GetFrequency, CNV frequency }
