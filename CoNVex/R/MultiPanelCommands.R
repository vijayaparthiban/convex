# TODO: Multipanel commands
# 
# Author: pv1
###############################################################################

`MultiPanelCommands` <- function(CNVcallsfile, cor_matrix, sample_info_file, features_file, gene_file, all_samples=0, Rbatch_folder="",output_folder=".", misc_str="") {
	
	if(CNVcallsfile == "") {
		print("Missing: CNV calls file!");
	}
	if(cor_matrix == "" & all_samples == 0) {
		print("Missing: correlation matrix file!");
	}
	
	d = read.table(CNVcallsfile, header=TRUE, sep="\t"); 
	sample_ids = unique(d$sample_id);
	
	dc_sids = data.frame();
	if(cor_matrix != "") {
		dc = read.table(cor_matrix, header=TRUE, sep="\t", check.names=FALSE, nrows=1);
		dc_sids = sample_ids[!(sample_ids %in% names(dc))];
	}
	if(output_folder==".") {
		output_folder = getwd();
	}
	
	if(length(dc_sids)>0) {
				write("The following sample IDs don't exist in the correlation matrix:", stderr());
				write(paste(dc_sids,collapse=", "));
	} else {
		if(Rbatch_folder=="") {Rbatch_folder = paste(getLibPath(),"/CoNVex/Rbatch/",sep="");}
		if(all_samples==0) {
			sample_col_ids = sapply(sample_ids,function(x) { return(which(names(dc)==x)) }) # in which column of cor_matrix does each sample_id exist?
			fn1 = paste(sep=",",CNVcallsfile,sample_ids,sample_col_ids,cor_matrix,sample_info_file,features_file,gene_file,output_folder,misc_str);
			MPPCommands = paste("R --no-save --no-restore --slave --args '",fn1,"' < ",Rbatch_folder,"MultiPanelCall.R", sep="")
		} else {
			fn2 = paste(sep=",",CNVcallsfile,sample_ids,sample_info_file,features_file,gene_file,output_folder,misc_str);
			MPPCommands = paste("R --no-save --no-restore --slave --args '",fn2,"' < ",Rbatch_folder,"MultiPanelCallAllSamp.R", sep="")
		}
	}
	return(MPPCommands);	
}
