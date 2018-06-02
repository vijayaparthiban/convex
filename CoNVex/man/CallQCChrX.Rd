\name{CallQCChrX}
\alias{CallQCChrX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Apply ChrX Call QC and return filtered CNV calls   }
\description{
Apply ChrX Call QC and return filtered CNV calls - ChrX only.
}
\usage{
CallQCChrX(cnv_calls = data.frame(), sample_meansmads = data.frame(), sd_medians = 2, abs_medianX_cutoff = 0.5, outlier_plot_file="")
}
\arguments{
  \item{cnv_calls}{ CNV calls (data frame) with mean log2 ratio, MAD of log2 ratio, common forward/backward and internal frequencies at any and 50 percent Chr overlap }
  \item{sample_meansmads}{ SampleMeans, SampleMADs, SampleMediansX and SampleMediansAuto AFTER sample QC (data frame) - output of SMMCallCommand() function }
  \item{sd_medians}{ default=2, For detecting outlier samples for ChrX. 'N' standard deviations from the lowess of SampleMediansX and SampleMediansAuto}
  \item{abs_medianX_cutoff}{ Extreme outlier examples: to explicitly remove gender 'U' (unknown) samples }
  \item{outlier_plot_file}{ Default="" (this function always plots samples' medians (and marks outliers if any). This option is for specifying file name explicitly (or ignore it to use a nice default name) }
}
\details{
CallQCChrX() applies call QC and returns filtered CNV calls. 
INPUT:
 - Input CNV calls data frame with the following columns (with 'column names' in brackets):
 - SampleMediansX vs. SampleMediansAuto
 - CoNVex score (convex_score)
 - Mean log2 ratio (mean_l2r)
 - MAD of log2 ratio
 - Internal frequency proportion at 50% reciprocal overlap (rc50_internal_freq)
 - common_forward and common_backward overlap proportion (output of recipeB package OR CoNVex's ChrOverlap() function)
 - Required Column names: "convex_score","mean_l2r","mad_l2r","internal_freq","rc50_internal_freq","common_forward","common_backward"
 - RETURNS CHRX CALLS AFTER QC
QC CHECK FOR:
  CASE 1 - >10Mb ChrX calls in any sample (remove all ChrX calls in those samples)
  CASE 2 - >0.5Mb calls in outlier samples (2xSD from the Lowess of SampleMediansX vs. SampleMediansAuto)  
  CASE 3 - All calls in gender 'U' (unknown) samples - remove all ChrX calls 
  - Apply all other Autosomal filtering criteria in ADDITION - call CallQC() function from here
REMOVAL CRITERIA:
  - Remove CNVs of samples matching Case 1 & Case 3
  - Remove CNVs of samples matching Case 2 (both 2a and 2b)
[ Possible extension: Mark Case 2a samples' ChrX as only match as 'FAIL' ]
}
\references{ CallQC ChrX }
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  
For ChrX, CallQCAll() and CallQCChrX() are used together. CallQCAll() must be used for all chromosomes.
}
\seealso{ 
CallQCAll(), SampleFilters()
}
\examples{
# > CallQCChrX()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ CallQC, ChrX }
