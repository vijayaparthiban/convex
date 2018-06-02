# TODO: SampleMeans
# 
# Author: pv1
###############################################################################

`SampleMeans` <- function(sample_ids, RDfiles, regions_file) {
	
	sortx = read.table(file=regions_file, header=TRUE); 
	allX = data.frame(seq(1:length(sortx[,1])), sortx); names(allX) = c("rowcount", "chr","start","end","GC","dG","Tm");
	
	# Separate Autosomes and X,Y
	allXS = allX[allX[,2] != 'Y',]; allXS = allXS[allXS[,2] != 'X',]; # Autosomes only - remove X and Y ######
	nrow_auto = nrow(allXS); rm(allX); gc();
	
	# Read read depth from each sample's files
	allXD = lapply(RDfiles, function(x) read.table(x, header=FALSE, check.names=FALSE, colClasses=c("character",rep("integer",3),"numeric"), nrows=nrow(allXS), comment.char="", sep="\t")[,5])
	allXS = data.frame(allXS,do.call('cbind',allXD)); names(allXS)[8:ncol(allXS)] = sample_ids; 
	rm(allXD); gc();
	
	SampleMeans = colMeans(allXS[,8:ncol(allXS)]);
	rm(sortx);  rm(allXS); # Cleanup
	gc();
	return(SampleMeans);
}

