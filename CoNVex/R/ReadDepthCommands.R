
`ReadDepthCommands` <- function(regions_file, sample_ids, bamfiles, bamindex_files=NULL, RDfiles, output_folder="", chr_prefix="", max_memory=2) {
	
	max_memory = paste("-Xmx",max_memory,"g",sep=""); current_folder = getwd();
	
	if(output_folder=="") { setwd(output_folder); output_folder = getwd(); } else { setwd(output_folder); output_folder = getwd(); setwd(current_folder); }
	
	if(is.null(bamindex_files)) {
		if(chr_prefix=="") {
			java_commands = paste("java ",max_memory," ReadDepth -bam_file ",bamfiles," -regions_file ",regions_file," -rd_file ",RDfiles,sep="");
		} else {
			java_commands = paste("java ",max_memory," ReadDepth -bam_file ",bamfiles," -chr_prefix ",chr_prefix," -regions_file ",regions_file," -rd_file ",RDfiles,sep="");
		}
	} else {
		if(chr_prefix=="") {
			java_commands = paste("java ",max_memory," ReadDepth -bam_file ",bamfiles," -bamindex_files ",bamindex_files," -regions_file ",regions_file," -rd_file ",RDfiles,sep="");
		} else {	
			java_commands = paste("java ",max_memory," ReadDepth -bam_file ",bamfiles," -chr_prefix ",chr_prefix," -bamindex_files ",bamindex_files," -regions_file ",regions_file," -rd_file ",RDfiles,sep="");
		}
	}
	return(java_commands);	
}
