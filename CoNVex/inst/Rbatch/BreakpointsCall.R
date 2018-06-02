# TODO: Breakspoints using a representative sample in the batch -- any sample may do
# 
# Author: pv1
###############################################################################

require(CoNVex);
ca = commandArgs(trailingOnly=TRUE);

fn = strsplit(ca,"\\,")
RDfile_RepSample = fn[[1]][1]
max_bin_size = as.integer(fn[[1]][2]) # 7 columns - chr, start, end, GC, dG, Tm
BP_file = fn[[1]][3]

if(!file.exists(RDfile_RepSample)) {
	write("Representative sample's read depth file is not valid!", stderr());
	q("no", status=2, runLast=FALSE);
}
if(is.na(max_bin_size)) {
	write("Max bin size is not valid!", stderr());
	write(paste("Max bin size 1",max_bin_size), stderr());
	q("no", status=2, runLast=FALSE);
}
if(!is.integer(max_bin_size)) {
	write("Max bin size is not valid!", stderr());
	write(paste("Max bin size 2",max_bin_size), stderr());
	q("no", status=2, runLast=FALSE);
}

BP = Breakpoints(RDfile_RepSample, max_bin_size)

if(!exists("BP")) {
	write("Breakspoints're not generated. Check input files and options!", stderr());
	q("no", status=2, runLast=FALSE);
}
write.table(BP, file=BP_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
