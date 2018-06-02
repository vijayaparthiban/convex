# TODO: Genome plots
# 
# Author: pv1
###############################################################################

`GenomePlot` <- function(GAMfile,sample_id="",output_folder) {
	d = read.table(as.character(GAMfile)); sample_id = as.character(sample_id); if(sample_id=="") {sample_id=GAMfile;}
	plotfile = paste(output_folder,"/GenomePlot_",sample_id,".png",sep="");
	chrs = unique(d[,1]);
	ptitl = paste("Genome ADM scores: ",sample_id,sep="");
	lgnd = c(sample_id,paste("MAD =",round(mad(d[,4],constant=1),2)));
	
	png(file=plotfile, width=960,height=480);
	plot(d[,8], pch=19, cex=0.2, xlab="Probe regions in genome order", ylab="ADM scores", main=ptitl, xaxt="n")
	for (j in chrs) abline(v=max(which(d[,1]==j)), lty=2, col="red");
	for (j in chrs) text(labels=paste(j), x=median(which(d[,1]==j)), y=2, col="red",cex=0.8);
	legend("topright",legend=lgnd,cex=0.8);
	dev.off()
}


