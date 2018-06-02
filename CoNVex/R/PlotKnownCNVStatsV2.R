# TODO: Plot Known CNV stats - This is to plot #calls vs known (%) in Consensus CNV list (bundled with CoNVex)
# 
# Author: pv1
###############################################################################

`PlotKnownCNVStatsV2` <- function(CallsAll,known_dels=data.frame(),known_dups=data.frame(),file_prefix="CNVstats_NumCalls_Known") {
	
	cdel = CallsAll[grep('DEL',CallsAll$cnv_type,ignore.case=TRUE),]; cdup = CallsAll[grep('DUP',CallsAll$cnv_type,ignore.case=TRUE),];	
	cdeltab = table(cdel$sample_id); cdeldf = data.frame(names(cdeltab),as.numeric(cdeltab)); names(cdeldf) = c("sample_id","num_dels");
	cduptab = table(cdup$sample_id); cdupdf = data.frame(names(cduptab),as.numeric(cduptab)); names(cdupdf) = c("sample_id","num_dups");	
	CallCounts = merge(x=cdeldf,y=cdupdf,by="sample_id");
	
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
	known_pc_df = data.frame(sample_ids,known_pc_all); names(known_pc_df) = c("sample_id","known_pc");

	CallCountsKnown = merge(x=CallCounts,y=known_pc_df,by="sample_id");
	allcalls = CallCountsKnown$num_dels+CallCountsKnown$num_dups;
	CCKA = data.frame(CallCountsKnown,allcalls); names(CCKA)[ncol(CCKA)] = "num_calls";
	
	callstats_file = paste(file_prefix,"_CallStats.txt",sep="");
	write.table(CCKA,file=callstats_file,row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
	
	png_file = paste(file_prefix,".png",sep="")	
	png(file=png_file,width=960,height=480)
	plot(allcalls,CallCountsKnown$known_pc,col=4,cex=0.7,main="#Calls vs Known CNVs (%) per Sample ",xlab="#Calls per sample",ylab="Known (%)", ylim=c(0,100))
	dev.off()
}