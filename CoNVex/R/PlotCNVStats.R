# TODO: Add comment
# 
# Author: pv1
###############################################################################

`PlotCNVStats` <- function(CNVcallsAll,file_prefix="CNVstats",outlier_cutoff=2) {
	
	cps = paste(file_prefix,"_","CallsperSample.png",sep="");
	dvd = paste(file_prefix,"_","DelDupRatio.png",sep="");
	
	cdel = CNVcallsAll[grep('DEL',CNVcallsAll$cnv_type,ignore.case=TRUE),]; cdup = CNVcallsAll[grep('DUP',CNVcallsAll$cnv_type,ignore.case=TRUE),];	
	cdeltab = table(cdel$sample_id); cdeldf = data.frame(names(cdeltab),as.numeric(cdeltab)); names(cdeldf) = c("sample_id","num_dels");
	cduptab = table(cdup$sample_id); cdupdf = data.frame(names(cduptab),as.numeric(cduptab)); names(cdupdf) = c("sample_id","num_dups");	
	CallCounts = merge(x=cdeldf,y=cdupdf,by="sample_id");
	
	png(file=cps,width=960,height=480)
	matplot(CallCounts[,2:3],ylim=c(0,max(CallCounts[,2:3])),type="l",col=c(2,4),lty=6,lwd=2,main="Number of CNV calls",xlab="Samples",ylab="Number of CNV calls");
	legend("topright",legend=c("Deletions","Duplications"),col=c(2,4),lty=6,lwd=2);
	dev.off()
	
	ddr = CallCounts[,2]/CallCounts[,3]; ddr[is.nan(ddr)] = 0; ddr[ddr == Inf] = 0; 
	outliers = which(ddr >= outlier_cutoff);
	outlier_names = as.character(CallCounts[outliers,1]); colours = rainbow(length(outliers));
	yl = c(0,3); if(max(ddr) > 3) { yl = c(0,ceiling(max(ddr))); }
	png(file=dvd,width=960,height=480)
	plot(ddr,pch=19,col=8,main="Ratio: Deletions/Duplications",xlab="Samples",ylab="Ratio", ylim=yl)
	if(length(outliers)!=0) {
		points(outliers,ddr[outliers],col=colours,pch=19)
		legend("topright",legend=outlier_names,col=colours,pch=19,cex=0.7);
	} 
	dev.off()	
}