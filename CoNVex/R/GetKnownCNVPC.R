# TODO: Return Known (%) in the consensus CNV list
# 
# Author: pv1
###############################################################################

`GetKnownCNVPC` <- function(CallsAll,known_dels=data.frame(),known_dups=data.frame()) {
	
	calls_dups = CallsAll[CallsAll$cnv_type=='DUP',]; calls_dels = CallsAll[CallsAll$cnv_type=='DEL',]; 
	if(nrow(known_dels)==0) {
		known_dels = read.table(paste(getLibPath(),"/CoNVex/extdata/known_dels_ccl2_1.txt",sep=""), header=FALSE); 
	}
	if(nrow(known_dups)==0) {
		known_dups = read.table(paste(getLibPath(),"/CoNVex/extdata/known_dups_ccl2_1.txt",sep=""), header=FALSE); 
	}
	cdups_pc = GetOverlap(known_dups,calls_dups); cdels_pc = GetOverlap(known_dels,calls_dels);
	calls_all = rbind(cdups_pc,cdels_pc); 
	calls_all = calls_all[order(calls_all$sample_id),]
	return(calls_all);
}

