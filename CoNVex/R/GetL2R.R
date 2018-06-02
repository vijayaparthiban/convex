# TODO: GetL2RADM
# 
# Author: pv1
###############################################################################

GetL2R <- function(features_file, GAMfiles) {
	allX = read.table(file=features_file, header=TRUE)[,1:4]; 
	GAMfiles_CutColumn = paste("cut -f4",GAMfiles); # column = "4,8"; 
	# Read L2R scores only from each sample's files
	allXD = lapply(GAMfiles_CutColumn, function(x) read.table(pipe(x), header=FALSE, check.names=FALSE, colClasses=c("numeric"), nrows=nrow(allX), comment.char="", sep="\t"));
	return(allXD);
}

