# TODO: Determines the breakpoints for #reads - to calculate region-wise MAD
# 
# Author: pv1
###############################################################################

`Breakpoints` <- function(RDfile_RepSample,max_bin_size=1000) {
	
	x <- read.table(file=RDfile_RepSample, header=FALSE, check.names=FALSE);
	xxc = data.frame(x[,1:3],x[,4]+0.01,x[,5]+0.01)
	
# Break number of reads into equal chunks of 1000
	max <- max_bin_size
	d1 <- split(sort(xxc[,4]), ceiling(seq_along(xxc[,4])/max))
	nmax <- c(1:length(d1))
	for(j in 1:length(d1)) { nmax[j] <- max(d1[[j]])  }
	nnmax <- c(0,nmax[1:length(nmax)-1],max(xxc[,4])+1)
	nnmax <- unique(nnmax)
	
# data frame of all required cols (l2r, reads, cut-info)
	xc <- cut(xxc[,4],breaks=nnmax)
	xxc_cut <- data.frame(xxc,xc)
	xcu <- sort(unique(xc))	
	Breakpoints <- data.frame(xcu,nnmax[1:length(nnmax)-1]); names(Breakpoints) = c("#ReadsInRanges","#Breakpoints");
	#write.table(Breakpoints, file=Breakpoints_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
	return(Breakpoints)
}