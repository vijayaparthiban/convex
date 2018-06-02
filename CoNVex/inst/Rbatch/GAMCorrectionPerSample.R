# TODO: This is to do GAM correction per sample -- called by command(s)
# 
# Author: pv1
###############################################################################

require(CoNVex);
require(mgcv);

ca = commandArgs(trailingOnly=TRUE);

fn = strsplit(ca,"\\,")
L2Rfile = fn[[1]][1]
features_file = fn[[1]][2]
GAMfile = fn[[1]][3] # bed file
RDfile = fn[[1]][4]
BPfile = fn[[1]][5]
gc_only = fn[[1]][6]

print(paste("Input log2 ratio:",L2Rfile));
print(paste("Output file corrected log2 ratio:",GAMfile));

ADMRep = GAMCorrection(L2Rfile,features_file,RDfile,BPfile,doNorm=0,gc_only)

if(!exists("ADMRep")) {
	write("GAM correction failed. Check the input/output options!", stderr());
	q("no", status=4, runLast=FALSE);
}
L2RfileRows = nrow(read.table(L2Rfile));

if(nrow(ADMRep) != L2RfileRows) {
	write("GAM correction failed. Check the input/output options!", stderr());
	q("no", status=4, runLast=FALSE);
}

# Write GAM corrected log2 ratio, Error-weighted scores
write.table(ADMRep, file=GAMfile, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
