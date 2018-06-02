\name{ReadDepthCommands}
\alias{ReadDepthCommands}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Unix Command Generator for ReadDepth java program }
\description{
Generates Unix commands (e.g., one per sample) to calculate read depth. 
}
\usage{
ReadDepthCommands(regions_file, sample_ids, bamfiles, \ 
bamindex_files=NULL, RDfiles, output_folder="", chr_prefix="", max_memory=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regions_file}{ Probe regions file}
  \item{sample_ids}{ Vector of Sample IDs }
  \item{bamfiles}{ Vector of BAM files - rows should match sample_ids }
  \item{bamindex_files}{ [optional] Vector of BAM index files - not required if they are present in the BAM file folder and use SAM/BAM format's accepted naming convention }
  \item{RDfiles}{ Vector of output read depth files }
  \item{output_folder}{	[optional] Output folder for read depth files }
  \item{max_memory}{ [optional] default=2 (for 2GB), Use according to the size of BAM files  }
}
\details{
Calculates (1) mean read depth of a probe region (any chr region should work) and 
		   (2) number of sequence reads falling in that region
}
\references{ Picard, args4j (links below) }
\author{ 
Parthiban Vijayarangakannan\cr
Wellcome Trust Sanger Institute\cr
Cambridge, UK 
}
\note{  CoNVex is a read depth based CNV detection algorithm. This function is useful to generate batch execution commands for a large scale project. The commands execute the java program ReadDepth which is bundled with CoNVex package.
CoNVex.jar, args4j-<version>.jar, and sam-<version>.jar should be in your CLASSPATH environment variable. 
}
\seealso{ 
\url{http://picard.sourceforge.net/} \cr
\url{http://java.net/projects/args4j/}\cr
\url{http://javarevisited.blogspot.co.uk/2011/01/how-classpath-work-in-java.html} \cr
\url{http://samtools.sourceforge.net/SAM1.pdf} \cr
}
\examples{
> rdc = ReadDepthCommands(regions_file=regions_file,sample_ids=sample_ids,\ 
bamfiles=bamfiles,RDfiles=RDfiles,chr_prefix="Chr",output_folder=output_folder,\ 
max_memory=2);
> rdc[1]
[1] "java -Xmx2g ReadDepth -bam_file /home/pv1/reads.474946.recal.bam \ 
 -chr_prefix Chr \ 
 -regions_file /home/pv1/R/x86_64-linux/2.11/CoNVex/extdata/SureSelect_50Mb.txt \ 
 -rd_file /home/pv1/UK10K/ProbeRD_UK10K12345.dat"
 >
 > system("java ReadDepth -help")
 -bam_file VAL      : bam or cram file location of the sample
 -bamindex_file VAL : (optional) Index file of the bam file, if exists in
                      different folder or in a different file name
 -chr_prefix VAL    :  (optional) Chr prefix in bam files [if you observe
                      'Chr1', use 'Chr' (without quotes). Default: no prefix
                      (expects 1,2,..X,Y)]
 -rd_file VAL       : Output file: No. of reads and mean depths of the given
                      probe regions
 -regions_file VAL  : File with Chr regions: chr (1,2,..X,Y), start, end
Required options are missing!
Usage: java ReadDepth -bam_file /path/to/file.bam -bamindex_file /path/to/bamindexfile.bai\ 
 -regions_file /path/to/file.txt \ 
 -rd_file /path/to/output/file/region_reads_depth.dat
Type -help or -<anyjunk> to display the options :)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ read depth, reads, read count, depth, picard }
