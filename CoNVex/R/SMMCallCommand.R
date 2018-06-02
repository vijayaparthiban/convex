# TODO: SMM Call Command - calls SampleMeansMads.R in the Rbatch_folder to calculate sample means and MADs
# 
# Author: pv1
###############################################################################

`SMMCallCommand` <- function(sample_info_file,regions_file,output_file,Rbatch_folder="") {
	fn = paste(sep=",",sample_info_file,regions_file,output_file);
	if(Rbatch_folder=="") {Rbatch_folder = paste(getLibPath(),"/CoNVex/Rbatch/",sep="");}
	SMMcommand = paste("R --no-save --no-restore --slave --args '",fn,"' < ",Rbatch_folder,"SampleMeansMads.R", sep="")
	return(SMMcommand);	
}
