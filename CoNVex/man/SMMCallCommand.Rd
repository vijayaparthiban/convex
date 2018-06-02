\name{SMMCallCommand}
\alias{SMMCallCommand}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Call SampleMeans and SampleMADsUnix through SampleMeansMads   }
\description{
SMMCallCommand creates R command with arguments to call SampleMeansMads. SampleMeansMads calculated sample mean depth and MAD(log ratio) using SampleMeans
and SampleMADs.  
}
\usage{
SMMCallCommand(sample_info_file,regions_file,output_file,Rbatch_folder="")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sample_info_file}{ Sample information file created by SampleInfoPrep (as part of you analysis) }
  \item{regions_file}{ Probe regions file }
  \item{output_file}{ Sample means and MADs will be stored in this file }
  \item{Rbatch_folder}{ [optional] Rbatch folder bundled with CoNVex for batch execution. Default='' }
}
\details{
SMMCallCommand uses sample information file (especially read depth files and GAM files) to retrieve read depth and log2 ratio 
and stores sample means and MADs in the output_file.
}
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  This is useful for SampleMeans vs. MADs plots that act as an visualisation of noise in the analysis batch.
}
\seealso{ 
SampleMADsUnix, SampleMeans, SampleMeansMads (in the Rbatch folder)
}
\examples{
> SMMCallCommand(sample_info_file="path/to/sample_info_file",regions_file="path/to/regions_file",output_file,Rbatch_folder="",)
> # execute the resulting command in the Unix prompt
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ SMMCallCommand, Mean Depth }
