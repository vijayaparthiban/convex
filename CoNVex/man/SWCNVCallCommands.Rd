\name{SWCNVCallCommands}
\alias{SWCNVCallCommands}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Unix Command Generator for using SWCNVCall R program }
\description{
Generates Unix command(s) (e.g., one per batch) to call CNVs using the Smith-Waterman algorithm
}
\usage{
SWCNVCallCommands(p,swt_del,swt_dup,dv,GAMfiles,CNVfiles,sample_id,\ 
centromere_regions_file,output_folder,sw_exec,Rbatch_folder="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ p value for Smith-Waterman algorithm }
  \item{swt_del}{ t value threshold CoNVex score (s) for the selection of deletion calls (default=5)  }
  \item{swt_dup}{ t value threshold CoNVex score (s) for the selection of duplication calls (default=5) }
  \item{dv}{ convex_score/(num_probes^dv) >= t is used for selecting CNV calls - num_probes to the power of dv is used (default=0.5)  }
  \item{GAMfiles}{	Vector of files containing corrected log2 ratio and error-weighted scores }
  \item{CNVfiles}{	Vector of output files to store CNV calls }
  \item{sample_ids}{ Vector of Sample IDs }
  \item{centromere_regions_file}{ Centromere regions file downloaded from the UCSC browser. Default: hg19 build's file bundled with CoNVex }
  \item{output_folder}{	Output folder for storing CNV calls - CNVcalls subfolder is ideal }
  \item{sw_exec}{ Smith-Waterman execution binary for change point detection - works similar to SW-array for array data }
  \item{Rbatch_folder}{	[optional] Rbatch folder bundled with CoNVex for batch execution. Rbatch_folder="" (default, works mostly) assumes that it's in R's library folder. You may copy it to somewhere else and modify this explicitly }
}
\details{
SWCNVCallCommands calls SWCNVCall R program available from Rbatch folder. It creates one Unix command per sample to call CNVs.
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
Journal of Molecular Biology 147: 195-197.

# Price, T.S., Regan, R., Mott, R., Hedman, A., Honey,
B., Daniels, R.J., Smith, L., Greenfield, A., Tiganescu,
A., Buckle, V., et al.(2005) SW-ARRAY: A dynamic programming
solution for the identification of copy-number changes in 
genomic DNA using array comparative genome hybridization data. Nucleic Acids Res. 33:3455-3464.
}
\seealso{ 
SWCNVCall, SWCNV\cr
SW-Array: \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1151590/} \cr
}
\examples{
> p=2; swt_del=5; swt_dup=5; dv=0.5; # Smith-Waterman algorithm / CoNVex parameters
> sw_exec = paste(getLibPath(),"/CoNVex/exec/swa_lin64",sep=""); # X86_64 linux distribution
> output_folder = "/home/pv1/UK10K/CNVcalls/"; # CNVcalls sub-folder is ideal
> sw_commands = SWCNVCallCommands(p,swt_del,swt_dup,dv,GAMfiles,sample_id,\ 
centromere_regions_file,output_folder,sw_exec,Rbatch_folder="");
> 
> sw_commands[1]
[1] "R --vanilla --slave --args '2,5,5,0.5,/home/pv1/UK10K/L2R/GAM_MOPD5095427.dat,\ 
MOPD5095427,/home/pv1/R/x86_64-linux/2.11/CoNVex/extdata/gaps_table_hg19.txt,\ 
/home/pv1/UK10K/CNVcalls/,\ 
/home/pv1/R/x86_64-linux/2.11/CoNVex/exec/swa_lin64' < \ 
/home/pv1/R/x86_64-linux/2.11/CoNVex/Rbatch/SWCNVCall.R"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ CNV calls, Smith-Waterman algorithm }
