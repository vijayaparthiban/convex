# TODO: Correct systematic errors caused by GC, dG and Tm
#
# Author: pv1
###############################################################################

`GAMCorrectionCommands` <- function(features_file,L2Rfiles,GAMfiles,RDfiles,BPfile,output_folder,sample_ids,gc_only=0,Rbatch_folder="") {	
	fn = paste(sep=",",L2Rfiles,features_file,GAMfiles,RDfiles,BPfile,gc_only);
	if(Rbatch_folder=="") {Rbatch_folder = paste(getLibPath(),"/CoNVex/Rbatch/",sep="");}
	SECcommands = paste("R --no-save --no-restore --slave --args '",fn,"' < ",Rbatch_folder,"GAMCorrectionPerSample.R", sep="")
	return(SECcommands);	
}