\name{CallQCAll}
\alias{CallQCAll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Apply Call QC and return filtered CNV calls  }
\description{
Apply Call QC and return filtered CNV calls - - must be used on autosomes and ChrX.
}
\usage{
CallQCAll(cnv_calls)
}
\details{
CallQCAll() applies call QC and returns filtered CNV calls. # Input CNV calls data frame with the following columns (with 'column names' in brackets):
 - CoNVex score (convex_score)
 - Mean log2 ratio (mean_l2r)
 - MAD of log2 ratio
 - Internal frequency proportion at 50 percent reciprocal overlap and any overlap (rc50_internal_freq, internal_freq)
 - common_forward and common_backward overlap proportion (output of recipeB package OR CoNVex's ChrOverlap() function)
 - Required Column names: "convex_score","mean_l2r","mad_l2r","internal_freq","rc50_internal_freq","common_forward","common_backward"
QC RULES:
REMOVE CNVs matching this criteria (apply on all CNVs):
- cnv_type=='DUP'  & convex_score<7 # [you can set this in Step 6: swt_dup=7 in the last step]
- cnv_size>500000 & convex_score<10
- cnv_size>200000 & convex_score<10 & num_probes >= 10
- mean_l2r >= 1.5
- common_forward == 0 & rc50_internal_freq2 > 0.05

REMOVE CNVs matching this criteria (apply on RARE CNVs. RARE: common_forward<0.5):
- cnv_type=='DUP' & num_probes=1 & convex_score<10 - 
- convex_score < 10 & num_probes == 1 & cnv_size < 500
- mean_l2r > -0.5 & cnv_type == 'DEL' & convex_score < 10
- mean_l2r < 0.25 & cnv_type == 'DUP' & convex_score < 10

Marked as 'FAIL' criteria (these CNVs are NOT removed):
- cnv_size<15000 & convex_score<10 & internal_freq > 0.01
- cnv_type=='DUP' & num_probes > 5 & convex_score < 8  
}
\references{ R Library path }
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  
For ChrX, CallQCAll() and CallQCChrX() are used together. CallQCAll() must be used for all chromosomes
}
\seealso{ 
CallQCChrX()
}
\examples{
# > CallQCAll()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ CallQC }
