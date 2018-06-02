# TODO: Plot Known CNV stats
# 
# Author: pv1
###############################################################################

`PlotKnownCNVStats` <- function(CallsAll,known_dels=data.frame(),known_dups=data.frame(),mark_outliers=0,outlier_from_median=5,file_prefix="CNVstats_Known") {
	
	calls_dups = CallsAll[CallsAll$cnv_type=='DUP',]; calls_dels = CallsAll[CallsAll$cnv_type=='DEL',]; 
	if(nrow(known_dels)==0) {
		known_dels = read.table(paste(getLibPath(),"/CoNVex/extdata/known_dels_ccl2_1.txt",sep=""), header=FALSE); 
	}
	if(nrow(known_dups)==0) {
		known_dups = read.table(paste(getLibPath(),"/CoNVex/extdata/known_dups_ccl2_1.txt",sep=""), header=FALSE); 
	}
	cdups_pc = GetOverlap(known_dups,calls_dups); cdels_pc = GetOverlap(known_dels,calls_dels);
	calls_all = rbind(cdups_pc,cdels_pc); names(calls_all)[ncol(calls_all)] = c("known_forward");
	calls_all = calls_all[order(calls_all$sample_id),];
	
	sample_ids = unique(calls_all$sample_id);
	known_pc_all = NULL;
	for(s in 1:length(sample_ids)) {
		d = calls_all[calls_all$sample_id==sample_ids[s],];
		known_pc = (length(which(d$known_forward>0))*100)/ nrow(d);
		known_pc = round(known_pc,digits=2);
		known_pc_all = c(known_pc_all,known_pc);
	}
	png_file = paste(file_prefix,".png",sep="")	
	png(file=png_file,width=960,height=480)
	plot(known_pc_all,pch=19,col=8,main="Known CNVs (%) per Sample ",xlab="Samples",ylab="Known (%)", ylim=c(0,100))
	if(mark_outliers!=0) {
		ocutoff = median(known_pc_all)-outlier_from_median;
		outliers = which(known_pc_all<=ocutoff); colours = rainbow(length(outliers));
		points(outliers,known_pc_all[outliers],col=colours,pch=19); outlier_names = as.character(sample_ids[outliers]);
		legend("bottomright",legend=outlier_names,col=colours,pch=19,cex=0.7);
	}
	dev.off();
}