# TODO: SampleLogRatio -- command generator
# 
# Author: pv1
###############################################################################

`SampleLogRatioCallCommands` <- function(sample_info_file,regions_file,features_file,cX=0,Rbatch_folder="",cor_matrix,version=1, RPKM=0, min_samples=25, sample_means_file="") {
	fn = paste(sep=",",sample_info_file,regions_file,features_file,cX,cor_matrix,version,RPKM,min_samples,sample_means_file);
	if(Rbatch_folder=="") {Rbatch_folder = paste(getLibPath(),"/CoNVex/Rbatch/",sep="");}
	SWACommands = paste("R  --no-save --no-restore --slave --args '",fn,"' < ",Rbatch_folder,"SampleLogRatioCall.R", sep="")
	return(SWACommands);	
}