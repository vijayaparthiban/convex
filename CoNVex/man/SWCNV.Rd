\name{SWCNV}
\alias{SWCNV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Change point detection using Smith-Waterman algorithm }
\description{
SWCNV uses pre-compiled swa_* binaries to call CNVs using the Smith-Waterman algorithm. This function acts as a wrapper for binary and extracts the calls through pipe() function.
Calls are filtered as per the tdel and tdup values (t value thresholds for CoNVex score). Returns a data frame containing deletions and duplications. This function is run for each sample
separately (one at a time). SWCNVCall, SWCNVCallCommands act as wrappers for batch execution.
}
\usage{
SWCNV(p,tdel,tdup,dv,GAMfile,sample_id,\ 
centromere_regions_file,output_folder,sw_exec)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ p value for Smith-Waterman algorithm }
  \item{tdel}{ t value threshold of CoNVex score (s) for the selection of deletion calls (default=5)  }
  \item{tdup}{ t value threshold of CoNVex score (s) for the selection of duplication calls (default=5) }
  \item{dv}{ convex_score/(num_probes^dv) >= t is used for selecting CNV calls - num_probes to the power of dv is used (default=0.5)  }
  \item{GAMfile}{ File containing corrected log2 ratio and error-weighted scores of a specific sample }
  \item{sample_id}{ Sample ID }
  \item{centromere_regions_file}{ Centromere regions file downloaded from the UCSC browser. Default: hg19 build's file bundled with CoNVex }
  \item{output_folder}{	Output folder for storing CNV calls - CNVcalls subfolder is ideal }
  \item{sw_exec}{ Smith-Waterman execution binary for change point detection - works similar to SW-array for array data }
}
\details{
SWCNV is called by SWCNVCall R program available from Rbatch folder. SWCNVCallCommands creates one Unix command per sample to call CNVs.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  This can be easily used together with batch queuing (e.g., bsub) programs. This function usually one Unix command per sample in the analysis.
}
\references{ 
# Smith, Temple F.; and Waterman, Michael S. (1981).
"Identification of Common Molecular Subsequences".
Journal of Molecular Biology 147: 195Ð197.

# Price, T.S., Regan, R., Mott, R., Hedman, A., Honey,
B., Daniels, R.J., Smith, L., Greenfield, A., Tiganescu,
A., Buckle, V., et al.(2005) SW-ARRAY: A dynamic programming
solution for the identification of copy-number changes in 
genomic DNA using array comparative genome hybridization data. Nucleic Acids Res. 33:3455Ð3464.
}
\seealso{ 
SWCNVCall
SWCNVCallCommands
}
\examples{
# Check the example of SWCNVCallCommands and the source code of SWCNVCall
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ CNV calls, Smith-Waterman algorithm }
