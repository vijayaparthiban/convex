# TODO: Correct systematic errors caused by GC, dG and Tm
#
# Author: pv1
###############################################################################

`BreakpointsCallCommands` <- function(RDfile_RepSample,max_bin_size=1000,BPfile,Rbatch_folder="") {	
	fn = paste(sep=",",RDfile_RepSample,max_bin_size,BPfile);
	if(Rbatch_folder=="") {Rbatch_folder = paste(getLibPath(),"/CoNVex/Rbatch/",sep="");}
	BPcommands = paste("R --no-save --no-restore --slave --args '",fn,"' < ",Rbatch_folder,"BreakpointsCall.R", sep="")
	return(BPcommands);
}

