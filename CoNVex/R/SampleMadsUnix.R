# TODO: Sample MADs, medians
# 
# Author: pv1
###############################################################################

`SampleMADsUnix` <- function(sample_ids, GAMfiles, regions_file) {
	
	sortx = read.table(file=regions_file, header=TRUE); 
	allX = data.frame(seq(1:length(sortx[,1])), sortx); names(allX) = c("rowcount", "chr","start","end","GC","dG","Tm");
	
	# Separate Autosomes and X,Y
	allXS = allX[allX[,2] != 'Y',]; num_AX = nrow(allXS) # Number of probes (rows) in Autosomes and ChrX
	num_X = length(which(allXS[,2] == 'X')); # Number of probes (rows) in Chr X only
	num_A = length(which(allXS[,2] != 'X')) # Number of probes (rows) in Autosomes only 
	rm(allX); rm(sortx); gc();
	
	# Read read depth from each sample's files
	allXD = lapply(GAMfiles, function(x) read.table(pipe(paste("cut -f4 ",x,"| head -n",num_AX,sep="")), header=FALSE, check.names=FALSE, colClasses=c("numeric"), nrows=nrow(allXS), comment.char="", sep="\t")[,1])
	allXS = data.frame(allXS,do.call('cbind',allXD)); names(allXS)[8:ncol(allXS)] = sample_ids; 
	rm(allXD); gc();
	
	SampleMads = apply(allXS[1:num_A,8:ncol(allXS)],2,mad,constant=1); # MADs in autosomes only
	SampleMediansAuto = apply(allXS[1:num_A,8:ncol(allXS)],2,median);
	SampleMediansX = apply(allXS[(num_A+1):num_AX,8:ncol(allXS)],2,median);
	SampleMadsMedians = data.frame(SampleMads,SampleMediansAuto,SampleMediansX)
    rm(allXS); gc(); # Cleanup
	return(SampleMadsMedians);
}
