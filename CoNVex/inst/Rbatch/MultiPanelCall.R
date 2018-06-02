# TODO: Multi Panel Plot - Batch caller
# 
# Author: pv1
###############################################################################

require(CoNVex);
ca = commandArgs(trailingOnly=TRUE);

fn = strsplit(ca,"\\,")
CNVfile = fn[[1]][1]
sample_id = fn[[1]][2] 
sample_column_id = fn[[1]][3] 
cor_matrix_file = fn[[1]][4] 
sample_info_file = fn[[1]][5] 
features_file = fn[[1]][6] 
gene_file = fn[[1]][7] 
output_folder = fn[[1]][8] 
mstr = fn[[1]][9] 
		
CNVcalls = read.table(CNVfile, header=TRUE, sep="\t"); ## All CNV calls
calls2plot = CNVcalls[CNVcalls$sample_id==sample_id,]; column_names = names(calls2plot); rm(CNVcalls); ## CNV calls of only one sample - sample_id

cor_matrix = read.table(cor_matrix_file, header=TRUE, sep="\t", check.names=FALSE, nrows=1); 
cor_matrix_names = names(cor_matrix); # rownames(cor_matrix) = cor_matrix_names;

unix_read_command = paste("cut -f",sample_column_id,cor_matrix_file)
cor_sample_us = read.table(pipe(unix_read_command), header=TRUE, sep="\t", check.names=FALSE)[,1]; 
rm(cor_matrix); gc(); names(cor_sample_us) = cor_matrix_names;
cor_sample = sort(cor_sample_us, decreasing=TRUE); 
cor_sample_sids = names(cor_sample);

mid = ""; fid = ""; mother_id = ""; father_id="";

if(any(column_names=="mother_id")) {
	mid = as.character(calls2plot$mother_id[1]);
}
if(any(column_names=="father_id")) {
	fid = as.character(calls2plot$father_id[1]);
}
highly_correlated = cor_sample_sids[1:25]; # top N excluding the first (same) sample 

sample_info = read.table(sample_info_file, header=FALSE, sep="\t"); GAMfiles = as.character(sample_info[,6]); 
sample_ids = as.character(sample_info[,1]); 

if(any(cor_sample_sids==mid)) {
	mother_id=mid; 	highly_correlated = c(highly_correlated,mid);
}
if(any(cor_sample_sids==fid)) {
	father_id=fid; 	highly_correlated = c(highly_correlated,fid);
}
sample_info_hc = sample_info[(sample_ids %in% highly_correlated),];

GAM2Read = as.character(sample_info_hc[,6]); sample_ids_read = as.character(sample_info_hc[,1]);

allXD = GetL2RADM(features_file, GAM2Read);
allX = read.table(features_file, header=TRUE, sep="\t")[,1:4]; allXscore = allX;
L2R <- lapply(allXD, function(z) z[c("V1")]); ADM <- lapply(allXD, function(z) z[c("V2")]); 
rm(allXD); gc(); # Clean up
allX = data.frame(allX,do.call('cbind',L2R)); names(allX)[5:length(allX)] = sample_ids_read;
allXscore = data.frame(allXscore,do.call('cbind',ADM)); names(allXscore)[5:length(allXscore)] = sample_ids_read;
rm(L2R); rm(ADM); gc(); # Clean up

ginfo = read.table(gene_file, header=TRUE);

setwd(output_folder);
for(i in 1:nrow(calls2plot)) {
	MultiPanelPlot(allX=allX,allXscore=allXscore,Fd_rep=calls2plot[i,],sample_id=sample_id,known_rstr="", gene_info=ginfo,misc_str=mstr,mother_id=mother_id,father_id=father_id)
}