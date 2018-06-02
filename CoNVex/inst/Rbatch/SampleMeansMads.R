# TODO: Call SampleMeans, SampleMADs
# 
# Author: pv1
###############################################################################

require(CoNVex);
ca = commandArgs(trailingOnly=TRUE);

fn = strsplit(ca,"\\,")
sample_info_file = fn[[1]][1]
regions_file = fn[[1]][2] # 7 columns - chr, start, end, GC, dG, Tm
output_file = fn[[1]][3]

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
GAMfiles = as.character(baminfoRD[,6]);

SampleMeans = SampleMeans(sample_ids, RDfiles, regions_file); SampleMeans = round(SampleMeans,digits=2);
SampleMadsMedians = SampleMADsUnix(sample_ids, GAMfiles, regions_file); SampleMadsMedians = round(SampleMadsMedians,digits=3);

SMM = data.frame(sample_ids,gender,SampleMeans,SampleMadsMedians);
names(SMM) = c("sample_id","Gender","SampleMeans","SampleMads","SampleMediansAuto","SampleMediansX")
		
write.table(SMM,file=output_file,row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE);
