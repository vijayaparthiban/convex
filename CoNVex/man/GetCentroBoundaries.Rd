\name{GetCentroBoundaries}
\alias{GetCentroBoundaries}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get Centromere Boundaries in Probe Regions }
\description{
This function determines centromere boundaries within probe regions by binning them separately for each chromosome arm (before and after centromere). 
While using the Smith-Waterman algorithm, CNVs (change points) are not merged across the centromere. This function is used internally by SWCNV.
}
\usage{
GetCentroBoundaries(d,centromeres)
}
\arguments{
  \item{d}{ Data frame having the following columns: chr, start, end, log2 ratio (corrected), MAD (regional)  }
  \item{centromeres}{ Centromere co-ordinates for each chromosome (as in UCSC) are supplied in 3 columns: chr, start, end }
}
\details{
GetCentroBoundaries
}
\references{ UCSC genome browser }
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  This function is internally called from SWCNV function. For different builds (hg19, hg18), use different centromere regions file to supply the co-ordinates for the 'centromere'. hg19 gaps table has been included with CoNVex package. 
For different organisms and chromosome builds, you can download the gaps table from the UCSC genome browser (link below).
How to download the gaps table? 
[1] Open the URL below ('see also' section) in a web browser
[2] In the 'group:' option, select 'All tables' 
[3] In the 'table:' option, select 'gap'.
[4] In the 'output format:' option, select 'All fields from select table'
[5] In the 'output file:' option, enter a file name
[6] Click 'get output' button to download the file
[7] Use the file for the centromere_regions_file option in SWCNV (or in SWCNVCall, SWCNVCallCommands appropriately)
[8] This function is normally called by SWCNV
}
\seealso{ 
SWCNV, SWCNVCall, SWCNVCallCommands\cr
Download UCSC gaps table from here:\cr
\url{http://genome.ucsc.edu/cgi-bin/hgTables} \cr

}
\examples{
> centro_boundaries = GetCentroBoundaries(d,hg_cent);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ GetCentroBoundaries, centromere regions, UCSC gaps table }
