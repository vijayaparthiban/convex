# TODO: CNVfiles contain the list of file names (one per sample) in which each sample's CNVs are stored.
# This functions returns a data frame of all CNVs in all samples. 
# Data frame Columns: "chr","start","end","num_probes","convex_score","cnv_type","sample_id"
#
# Author: pv1
###############################################################################

`GetCNVCalls` <- function(CNVfiles, names=1) {
	CallsX = lapply(CNVfiles, function(x) read.table(x, header=FALSE, sep="\t")) #### This is for p3 calls only
	CallsAll = data.frame(do.call('rbind',CallsX)); rm(CallsX); gc();
	if(names==1) {
		names(CallsAll)[1:7] = c("chr","start","end","num_probes","convex_score","cnv_type","sample_id");	
	} else if(names==2) {
		names(CallsAll)[1:9] = c("chr","start","end","num_probes","convex_score","cnv_type","sample_id","mean_l2r","mad_l2r");	
	}
	return(CallsAll);
}