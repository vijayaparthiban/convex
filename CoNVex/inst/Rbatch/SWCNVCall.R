# TODO: Call CNVs from SW Array 
# 
# Author: pv1
###############################################################################

# Input values
# p, t, dv - SW Array parameters
require(CoNVex);

ca = commandArgs(trailingOnly=TRUE); fn <- strsplit(ca,"\\,");

p = as.numeric(fn[[1]][1]); tdel = as.numeric(fn[[1]][2]); tdup = as.numeric(fn[[1]][3]); dv = as.numeric(fn[[1]][4]);  # SW Array parameters
GAMfile = fn[[1]][5]; # GAM file
cnv_calls_file = fn[[1]][6]
sample_id = fn[[1]][7];
centromere_regions_file = fn[[1]][8]; 
output_folder = fn[[1]][9]; # Output folder
sw_exec = fn[[1]][10];
version = as.numeric(fn[[1]][11]);

# Get the required functions
setwd(output_folder);

# Test print: print(p); print(tdel); print(tdup); print(dv); print(ADM3_reports_file); print(sample_id); print(output_folder);

Fd_report = SWCNV(p,tdel,tdup,dv,GAMfile,sample_id,centromere_regions_file,output_folder,sw_exec);

#cnv_calls_file = paste(getwd(),"/CoNVex_",sample_id,"_p",p,"_tdel",tdel,"_tdup",tdup,"_dv",dv,"_.txt",sep=""); # CNV calls / output file names

if(nrow(Fd_report)==0) {
	if(version==1) {
		write(paste("Failed to detect deletions/duplications (either or both) in this sample:",sample_id), stderr());
		write("CNV detection failed! Check sample data (noise, quality, mean depth, etc.)", stderr());
		q("no", status=5, runLast=FALSE);
	} else {
		
	}
}

if(!exists("Fd_report")) {
	write("CoNVex CNV detection failed. Check input/output options!", stderr());
	q("no", status=5, runLast=FALSE);
}	
write.table(Fd_report, file=cnv_calls_file, row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t")
print(paste("Found",nrow(Fd_report), "CNVs; Written file:")) # Print message
print(cnv_calls_file); 



