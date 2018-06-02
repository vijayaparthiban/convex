# TODO: MeansMads commands
# 
# Author: pv1
###############################################################################

`MeansMadsCommands` <- function(CNVfiles, sample_ids, cor_matrix, sample_info_file, regions_file, features_file, all_samples=0, Rbatch_folder="",output_folder=".", cnv_files_header=0) {

	if(cor_matrix == "" & all_samples == 0) {
		print("Missing: correlation matrix file!");
	}
	if(Rbatch_folder=="") {Rbatch_folder = paste(getLibPath(),"/CoNVex/Rbatch/",sep="");}
	output_files = paste(substr(CNVfiles,1,nchar(CNVfiles)-4),"MeansMads.txt",sep="");
	
	if(all_samples==0) {
		fn1 = paste(sep=",",CNVfiles,sample_ids,cor_matrix,sample_info_file,regions_file,features_file,output_folder,output_files, cnv_files_header);
		MPPCommands = paste("R --no-save --no-restore --slave --args '",fn1,"' < ",Rbatch_folder,"MeansMadsCall.R", sep="")
	} else {
		fn2 = paste(sep=",",CNVfiles,sample_ids,sample_info_file,regions_file,features_file,gene_file,output_folder,output_files, cnv_files_header);
		MPPCommands = paste("R --no-save --no-restore --slave --args '",fn2,"' < ",Rbatch_folder,"MeansMadsCallAllSamp.R", sep="")
	}
	return(MPPCommands);	
}
