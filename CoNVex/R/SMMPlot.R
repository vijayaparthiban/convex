# TODO: Add comment
# 
# Author: pv1
###############################################################################


SMMPlot <- function(smm_file="",output_plot_file="SampleMeansMadsPlot.png",MAD_Cutoff=0.3) {
	
	if(smm_file=="") {
		print("Please specify the SampleMeans / MADs file - the output from SMMCallCommand()")
	} else {
		d = read.table(smm_file,header=TRUE);
		png(file=output_plot_file,width=720,height=720);
		plot(d[,2],d[,3],pch=19,col=4,cex=0.7,xlab="Sample Means", ylab="MADs", main="Sample Means vs. MADs");
		dd = d[d[,3]>MAD_Cutoff,];
		points(dd[,2],dd[,3],col=2,cex=0.8)
		text(x=dd[,2],y=dd[,3],label=as.character(dd[,1]),pos=4,cex=0.8,col=2)
		dev.off();
	}	
}
