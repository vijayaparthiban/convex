# TODO: Multi Panel Plot - Batch caller
# 
# Author: pv1
###############################################################################

require(CoNVex);
ca = commandArgs(trailingOnly=TRUE);

fn = strsplit(ca,"\\,")
CNVfile = fn[[1]][1]
sample_id = fn[[1]][2] 
sample_info_file = fn[[1]][3] 
features_file = fn[[1]][4] 
gene_file = fn[[1]][5] 
output_folder = fn[[1]][6] 
mother_id = ""; father_id="";

CNVcalls = read.table(CNVfile, header=TRUE, sep="\t"); ## All CNV calls
calls2plot = CNVcalls[CNVcalls$sample_id==sample_id,]; column_names = names(calls2plot); rm(CNVcalls); ## CNV calls of only one sample - sample_id

sample_info = read.table(sample_info_file, header=FALSE, sep="\t"); 
chosen_samples = c(sample_id, mother_id, father_id);
si = sample_info[sample_info[,1] %in% chosen_samples,];

if(nrow(sample_info) > 25) {
	GAMfiles = as.character(sample_info[1:25,6]); 
	sample_ids = as.character(sample_info[1:25,1]); 
} else {
	GAMfiles = as.character(sample_info[,6]); 
	sample_ids = as.character(sample_info[,1]); 
}
	GAMfiles = c(GAMfiles, as.character(si[,6]));
	sample_ids = c(sample_ids, as.character(si[,1]));

allXD = GetL2RADM(features_file, GAMfiles);
allX = read.table(features_file, header=TRUE, sep="\t")[,1:4]; allXscore = allX;
L2R <- lapply(allXD, function(z) z[c("V1")]); ADM <- lapply(allXD, function(z) z[c("V2")]); 
rm(allXD); gc(); # Clean up
allX = data.frame(allX,do.call('cbind',L2R)); names(allX)[5:length(allX)] = sample_ids;
allXscore = data.frame(allXscore,do.call('cbind',ADM)); names(allXscore)[5:length(allXscore)] = sample_ids;
rm(L2R); rm(ADM); gc(); # Clean up

ginfo = read.table(gene_file, header=TRUE);

setwd(output_folder);
for(i in 1:nrow(calls2plot)) {
	MultiPanelPlot(allX=allX,allXscore=allXscore,Fd_rep=calls2plot[i,],sample_id=sample_id,known_rstr="", gene_info=ginfo,misc_str="",mother_id=mother_id,father_id=father_id)
}

