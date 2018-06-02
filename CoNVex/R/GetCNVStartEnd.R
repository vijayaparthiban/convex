# TODO: Add comment
# 
# Author: pv1
###############################################################################


`GetCNVStartEnd` <- function(d,Fd_report) {
	
	dd_swa = NULL; 
	for(i in 1:length(Fd_report[,1])) {
		tmp_swa = d[as.character(d[,2]) == as.character(Fd_report[i,1]) & d[,3] >= Fd_report[i,2] & d[,3] <= Fd_report[i,3],]
		tmp_swa = data.frame(tmp_swa[1,2],min(tmp_swa[,3]),max(tmp_swa[,4]))
		dd_swa = rbind(dd_swa,tmp_swa)
	}
	Fd_reportX = data.frame(dd_swa,Fd_report[,4:5]); names(Fd_reportX) = c("chr","Start","End");
	return(Fd_reportX);
}
