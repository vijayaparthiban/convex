# TODO: SampleInfoPrep prepares the information for analysis:
# * Requires input info as shown in function options
# * Creates output_folder if it doesn't exist
# * Check whether the BAM files exist
# * Checks gender - only 'F', 'M' and 'U' are allowed
# * Checks whether sample_info_file exists already - overwrite=TRUE (default: FALSE) overwrites it without warning
# * Saves the sample_info_file in the specified output_folder
# Author: pv1
###############################################################################

`SampleInfoPrepInteractive` <- function(sample_ids,gender,bamfiles,RDfiles,L2Rfiles,GAMfiles,CNVfiles,sample_info_file,output_folder,overwrite=FALSE) {

	sample_ids = as.character(sample_ids); gender = as.character(gender); bamfiles =  as.character(bamfiles); 
	RDfiles = as.character(RDfiles); L2Rfiles = as.character(L2Rfiles); GAMfiles = as.character(GAMfiles); CNVfiles = as.character(CNVfiles);
	
	if(!file.exists(output_folder)) {
		dir.create(output_folder); 
	}
	# Create subfolders
	if(!file.exists(output_folder)) {
		write("Specified output folder was neither there no could be created", stderr());
		write("Check your disk access/availability!", stderr());
	} else {
		setwd(output_folder); 
	} 
	if(!file.exists("L2R")) {
		dir.create("L2R"); 
	}
	if(!file.exists("CNVcalls")) {
		dir.create("CNVcalls");
	}
	checkbams = file.exists(bamfiles); 
	if(length(which(checkbams==FALSE)) > 0) {
		write("BAM files of these samples are missing:", stderr());
		write(paste("Sample ID:",sample_ids[!checkbams],"-->",paste(bamfiles[!checkbams])), stderr());
		#q("no", status=1, runLast=FALSE);
	}
	allowed_gender = c('F','M','U');
	gender_check = is.element(gender,allowed_gender);
	
	if(any(gender_check==FALSE)) {
		write("Only 'F','M' and 'U' are allowed in gender -- for Female, Male and Unknown respectively!", stderr());
		#q("no", status=1, runLast=FALSE);
	}
	# Proper characters
	baminfoRD = data.frame(sample_ids,gender,bamfiles,RDfiles,L2Rfiles,GAMfiles,CNVfiles);
	if(!exists("baminfoRD")) {
		write("Sample info not generated. Check the input files and options!", stderr());
		write("Check: No. of sample_ids, bamfiles, etc. match each other!", stderr());
		#q("no", status=2, runLast=FALSE);
	}
	
	if(overwrite) {
		write.table(baminfoRD,file=sample_info_file,row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE);
	} else {
		if(file.exists(sample_info_file)) {
			write("Sample info file already exists! Make sure that you're not accidentally overwriting this", stderr());
			write("Use 'overwrite=TRUE' to overwrite! [OR] Move/rename/delete the existing file", stderr());
			#q("no", status=1, runLast=FALSE);			
		} else {
			write.table(baminfoRD,file=sample_info_file,row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE);
		}
	}
	return(baminfoRD)
}