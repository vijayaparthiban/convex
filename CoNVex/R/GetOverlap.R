# TODO: Chr Overlap between two data frames
# RETURNED: % overlap of chr_list2 
# Author: pv1
###############################################################################

`GetOverlap` <- function(chr_list1, chr_list2, tmp_folder="/tmp/") {
	random_str = paste(format(Sys.time(), "%a%b%d_%H%M"),"_",sample(1:100000,1),"_",sample(1:100000,1),sep="");
	tmp_file1 = paste(tmp_folder,"chr_list1_",random_str,".txt",sep="");
	tmp_file2 = paste(tmp_folder,"chr_list2_",random_str,".txt",sep="");
	
	write.table(chr_list1[,1:3],file=tmp_file1, row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t");
	write.table(chr_list2[,1:3],file=tmp_file2, row.names=FALSE, col.names=FALSE, quote=FALSE,sep="\t");
	OVgen = paste("java ChrOverlap -chr_list1",tmp_file1,"-chr_list2",tmp_file2);
	OVgenPC = as.numeric(system(OVgen,intern=TRUE));
	OVgenDF = data.frame(chr_list2,OVgenPC); names(OVgenDF)[ncol(OVgenDF)] = "forward_overlap";
	return(OVgenDF);
}
