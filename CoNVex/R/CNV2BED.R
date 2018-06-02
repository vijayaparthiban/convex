# TODO: Create a bed file from CNVcalls -- to be uploaded as a custom track in the UCSC genome browser
# 
# Author: pv1
###############################################################################

`CNV2BED` <- function(CNVcalls=data.frame(),Bedfile="CoNVexCNVs_UCSC.bed", db="hg19") {
	
	if(nrow(CNVcalls) == 0) {
		print("Please provide CNV calls in a data frame");
	} 
	cnvcalls_names = c("chr","start","end","num_probes","convex_score","cnv_type","sample_id");
	nn = unique(names(CNVcalls)==cnvcalls_names);
	
	if(nn && length(nn)==1) {
		# Appropriate bed file columns
		excol = paste("CoNVex score = ",CNVcalls$convex_score,"; ","NumProbes = ",CNVcalls$num_probes,"; CNV: ",CNVcalls$cnv_type,sep="");
		bedgraph = data.frame(paste("chr",CNVcalls$chr,sep=""),CNVcalls$start,CNVcalls$end,CNVcalls$sample_id,CNVcalls$convex_score,excol,CNVcalls$cnv_type);
		names(bedgraph) = c("chr","start","end","sample_id","convex_score","ucsc_desc","cnv_type");
		
		bg_dels = bedgraph[bedgraph$cnv_type=='DEL',1:6]; bg_dups = bedgraph[bedgraph$cnv_type=='DUP',1:6];
		
		track_defaults = c("browser position chr1:14000-9650000","browser hide all","browser pack refGene encodeRegions")
		track_dels = paste("track type=bedDetail name='DELS' color=255,0,0 description='CoNVex Deletions' visibility=full db=",db,sep="")
		track_dups = paste("track type=bedDetail name='DUPS' color=0,0,255 description='CoNVex Duplications' visibility=full db=",db,sep="")
		
		# Write top lines
		lines_at_top = c(track_defaults,track_dels);
		write.table(lines_at_top, file = Bedfile, quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE);		
		# Write deletions only - append=TRUE
		write.table(bg_dels, file = Bedfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append=TRUE);		
		# Write duplication track name - append=TRUE
		write.table(track_dups, file = Bedfile, quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE, append=TRUE);		
		# Write duplication only - append=TRUE	
		write.table(bg_dups, file = Bedfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append=TRUE);		
	} else {
		print("names of the columns indicate that you don't have columns in the required format");
		print("Please use GetCNVCalls() to read the calls into a data frame!");
	}
}
