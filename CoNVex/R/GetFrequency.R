# TODO: Get CNV frequency
# RETURNED: % overlap of chr_list2 
# Author: pv1
###############################################################################

`GetFrequency` <- function(chr_list1=data.frame(), chr_list2=data.frame(), sid_col1=7, sid_col2=7, ro_threshold=0, tmp_folder="/tmp/") {
	random_str = paste(format(Sys.time(), "%a%b%d_%H%M"),"_",sample(1:100000,1),"_",sample(1:100000,1),sep="");
	tmp_file1 = paste(tmp_folder,"chr_list1_",random_str,".txt",sep="");
	tmp_file2 = paste(tmp_folder,"chr_list2_",random_str,".txt",sep="");
	write.table(chr_list1,file=tmp_file1, row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t");
	
	if(nrow(chr_list2)==0) {
		OVgen = paste("java CNVfrequency -chr_list1",tmp_file1,"-sampleID_col1",sid_col1,"-ro_threshold",ro_threshold);
		OVgenPC = as.numeric(system(OVgen,intern=TRUE));
		OVgenDF = data.frame(chr_list1,OVgenPC,(OVgenPC/length(unique(chr_list1[,sid_col1])))); names(OVgenDF)[(ncol(OVgenDF)-1):ncol(OVgenDF)] = c("frequency1","frequency2");
	} else {
		write.table(chr_list2,file=tmp_file2, row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t");
		OVgen = paste("java CNVfrequency -chr_list1",tmp_file1,"-chr_list2",tmp_file2,"-sampleID_col1",sid_col1,"-sampleID_col2",sid_col2,"-ro_threshold",ro_threshold);
		OVgenPC = as.numeric(system(OVgen,intern=TRUE));
		OVgenDF = data.frame(chr_list1,OVgenPC,(OVgenPC/length(unique(chr_list2[,sid_col2])))); names(OVgenDF)[(ncol(OVgenDF)-1):ncol(OVgenDF)] = c("frequency1","frequency2");
	}
	return(OVgenDF);
}
