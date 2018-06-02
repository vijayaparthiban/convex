# TODO: Smith-Waterman based CNV calling commands
# 
# Author: pv1
###############################################################################

`SWCNVCallCommands` <- function(p,swt_del,swt_dup,dv,GAMfiles,CNVfiles,sample_id,centromere_regions_file,output_folder,sw_exec,Rbatch_folder="",version=1) {
	fn = paste(sep=",",p,swt_del,swt_dup,dv,GAMfiles,CNVfiles,sample_ids,centromere_regions_file,output_folder,sw_exec,version);
	if(Rbatch_folder=="") {Rbatch_folder = paste(getLibPath(),"/CoNVex/Rbatch/",sep="");}
	SWACommands = paste("R --no-save --no-restore --slave --args '",fn,"' < ",Rbatch_folder,"SWCNVCall.R", sep="")
	return(SWACommands);	
}