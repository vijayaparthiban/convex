\name{CNV2BED}
\alias{CNV2BED}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Create UCSC uploadable bed file from CNV calls }
\description{
This functions accepts the data frame created by GetCNVCalls() function and creates a bed file in a given name. You can either upload this to UCSC genome 
browser for a custom track or use it for other purposes.
}
\usage{
CNV2BED(CNVcalls=data.frame(),Bedfile="CoNVexCNVs_UCSC.bed", db="hg19")
}
\arguments{
  \item{CNVcalls}{ CNV calls from all samples - output from GetAllCalls() }
  \item{Bedfile}{ Name of the bed file to be created. Default: CoNVexCNVs_UCSC.bed }
  \item{db}{ Genome build; Default: Assumes that the coordinates are in 'hg19' format}
}
\details{
GetCNVCalls() function returns columns in this format: "chr", "start", "end", "num_probes", "convex_score", "cnv_type", "sample_id". This 
format is expected. DELS and DUPS are listed separately and coloured in red and blue respectively.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  
After creating the basic bed file, you can open it and modify few things to suit your purpose. 
}
\seealso{ 
GetCNVCalls
}
\examples{
UK10KCNVCalls = GetCNVCalls(CNVfiles=cnv_call_files)
CNV2BED(CNVcalls=UK10KCNVCalls)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ CNV2BED, UCSC bed file, genome browser}
