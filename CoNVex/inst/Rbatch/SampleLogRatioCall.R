# TODO: This is to do GAM correction per sample for -- called by command(s)
# 
# Author: pv1
###############################################################################

# Input command line arguments:
# sampl2r - l2r file name (all file names w/ full path)
# allXSfeat - all features dG, Tm, GC, etc.
# gamfile - gam file name
require(CoNVex);
ca = commandArgs(trailingOnly=TRUE);

fn = strsplit(ca,"\\,")
sample_info_file = fn[[1]][1]
regions_file = fn[[1]][2] # 7 columns - chr, start, end, GC, dG, Tm
features_file = fn[[1]][3] # 7 columns - chr, start, end, GC, dG, Tm
cX = fn[[1]][4]
cor_matrix_file = fn[[1]][5]
version = as.numeric(fn[[1]][6])
rpkm = as.numeric(fn[[1]][7])
min_samp = as.numeric(fn[[1]][8])
sample_means_file = fn[[1]][9]

if(!file.exists(regions_file)) {
	write("Regions file is not valid!", stderr());
	q("no", status=2, runLast=FALSE);
}
if(!file.exists(sample_info_file)) {
	write("Sample info file is not valid!", stderr());
	q("no", status=2, runLast=FALSE);
}

baminfoRD = read.table(file=sample_info_file,header=FALSE)
sample_ids = as.character(baminfoRD[,1]);
gender = as.character(baminfoRD[,2]);
RDfiles = as.character(baminfoRD[,4]);
L2Rfiles = as.character(baminfoRD[,5]);

if(version==1) {
	SLR = SampleLogRatio(RDfiles, sample_ids, gender, regions_file, chrX=as.numeric(cX),sample_means_file)
} else if(version==2) {
	SLR = SampleLogRatioV2(RDfiles, sample_ids, gender, regions_file, chrX=as.numeric(cX), min_samples=min_samp, RPKM=rpkm, sample_means_file=sample_means_file)
} else if(version==3) {
	SLR = SampleLogRatioV3(RDfiles, sample_ids, gender, regions_file, chrX=as.numeric(cX), cor_matrix_file, min_samples=min_samp, RPKM=rpkm, sample_means_file=sample_means_file)
} 
if(!exists("SLR")) {
	write("Sample log ratio not generated. Check the input files and options!", stderr());
	q("no", status=2, runLast=FALSE);
}
if(length(sample_ids) != ncol(SLR[,8:ncol(SLR)])) {
	write("Data does not match with number of sample IDs!", stderr());
	q("no", status=2, runLast=FALSE);
}
write.table(SLR[,1:7], file = features_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Output log2 ratio of each sample are stored in a separate file
for(i in 1:length(sample_ids)) {
	write.table(SLR[,i+7], file = L2Rfiles[i], quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}
