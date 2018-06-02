# TODO: Multi panel plot
# 
# Author: pv1
###############################################################################

`MultiPanelPlot` <- function(allX,allXscore,Fd_rep,sample_id,known_rstr="", gene_info,misc_str="",mother_id="",father_id="") {
	
	getGeneSymbols <- function(x) {	
		x1 = x[!is.na(x)];	x2 = unique(x1); g1 = paste(x2,collapse="",sep="");	return(g1);
	}
	
	getKnownCNVString <- function(CallsAll,known_dels=data.frame(),known_dups=data.frame()) {
		
		#calls_dups = CallsAll[CallsAll$cnv_type=='DUP',]; calls_dels = CallsAll[CallsAll$cnv_type=='DEL',]; 
		if(nrow(known_dels)==0) {
			known_dels = read.table(paste(getLibPath(),"/CoNVex/extdata/known_dels_ccl2_1.txt",sep=""), header=FALSE); 
		}
		if(nrow(known_dups)==0) {
			known_dups = read.table(paste(getLibPath(),"/CoNVex/extdata/known_dups_ccl2_1.txt",sep=""), header=FALSE); 
		}
		cdels_pc = GetOverlap(CallsAll,known_dels); known1 = cdels_pc[cdels_pc[,ncol(cdels_pc)]>0,];
		cdups_pc = GetOverlap(CallsAll,known_dups); known2 = cdups_pc[cdups_pc[,ncol(cdups_pc)]>0,];
		
		ccc = 1; 	
		if(nrow(known1) > 0) {
			known_str <- paste(known1[,2],"-",known1[,3],sep="",collapse=", "); known_str <- paste("Events: ",known_str,sep="");
			known1_len = known1[,3]-known1[,2]; 
			kstr1 = paste("Known Dels: ",length(known1[,1])," events, ranging from ",min(known1_len),"-",max(known1_len),"bp in size",sep="")
			krmat1 = data.frame(known1[,2],known1[,3],ccc)
		} else {
			kstr1 = "Known Dels: NA"
			krmat1 = data.frame(0,0,ccc)
		}

		ccc = 2; 
		if(nrow(known2) > 0) {
			knownd_str = paste(known2[,2],"-",known2[,3],sep="",collapse="; "); knownd_str <- paste("Events: ",knownd_str,sep="");
			known2_len = known2[,3]-known2[,2]; 
			kstr2 = paste("Known Dups: ",length(known2[,1])," events, ranging from ",min(known2_len),"-",max(known2_len),"bp in size",sep="")
			krmat2 = data.frame(known2[,2],known2[,3],ccc)
		} else {
			kstr1 = "Known Dups: NA"
			krmat2 = data.frame(0,0,ccc)
		}
		names(krmat1) = c("start","end","cnv")
		names(krmat2) = c("start","end","cnv")	
		return(rbind(krmat1,krmat2));
	}
	
	chr = Fd_rep$chr; gene_start = Fd_rep$start; gene_end = Fd_rep$end; no_of_exons <- Fd_rep$num_probes; score = round(Fd_rep$convex_score, digits=2)
	diff = gene_end-gene_start; cnv = as.character(Fd_rep$cnv_type)
	lgnd1 <- paste(gene_start,"-",gene_end," (",(diff+1),"bp) / ",no_of_exons," probe regions / CoNVex score = ",score,sep="") # Legend
	plotfile = paste("CoNVexMultiPanel_",cnv,"_",sample_id,"_Chr",chr,"_",gene_start,"_",gene_end,".png",sep="")
	
	if(cnv=="DUP") {
		ptitl <- paste("CoNVex duplication call: Chr",chr,"in Sample",sample_id)
		colour = "blue";
	} else {
		ptitl <- paste("CoNVex deletion call: Chr",chr,"in Sample",sample_id)
		colour = "red";
	}
	
	all100 <- allX[as.character(allX[,2])==as.character(chr) & (allX[,3]>(gene_start-diff)) & (allX[,4]<(gene_end+diff)),];
	all100call <- allX[as.character(allX[,2])==as.character(chr) & (allX[,3]>=gene_start) & (allX[,4]<=gene_end),]; 
	
	if(nrow(all100)==1) { # only one row
		rowtmp <- which(allX[,1]==all100[,1])
		if(rowtmp <= 5) {
			rowtmps <- 1
		} else {
			rowtmps <- rowtmp-5
		}
		all100 <- allX[rowtmps:(rowtmp+5),]
		all100 <- all100[as.character(all100[,2])==as.character(chr),]					
	}
	if(nrow(all100)<=5) { # Centre point is set at the 2nd row
		rowtmp <- which(allX[,1]==all100[2,1])
		if(rowtmp <= 5) {
			rowtmps <- 1
		} else {
			rowtmps <- rowtmp-5
		}
		all100 <- allX[rowtmps:(rowtmp+5),]
		all100 <- all100[as.character(all100[,2])==as.character(chr),]			
	}
	all100score <- allXscore[allXscore[,1] %in% all100[,1],]
	
	png(file=plotfile, width=960, height=960)
	colmat = rbind(c(1:4),c(552,26,525,403),c(9,8,7,6)) # colour matrix: Rows 1,2,3: cnvtype, colour, y-location in plot
	mat = matrix(c(1,2,3), nrow=3,ncol=1) # Multipanel stuff
	par(mar = c(0,4,2,1))
	layout(mat, heights=c(1/3,1/3,1/3))
	
	nm1 = all100[,c(3,5:length(all100))]; nm2 = all100[,4:length(all100)]; names(nm1)[1] = "chrpos"; names(nm2)[1] = "chrpos"; nm = rbind(nm1,nm2); nm = nm[with(nm,order(nm[,1])),];
	cnvrange = c(min(all100call[,3]), max(all100call[,4])); 
	if(length(which(nm[,1] < cnvrange[1]))==0) {cnvmaxrange1 = cnvrange[1]; } else {cnvmaxrange1 = max(nm[nm[,1] < cnvrange[1],1])}
	if(length(which(nm[,1] > cnvrange[2]))==0) {cnvmaxrange2 = cnvrange[2]; } else {cnvmaxrange2 = min(nm[nm[,1] > cnvrange[2],1])}
	cnvmaxrange =  c(cnvmaxrange1,cnvmaxrange2);
	
	matplot(nm[,1], nm[,2:length(nm)], type="l", lty=1, xlab="", ylab="log2 ratio", main=ptitl, col=8, xaxt="n",cex.axis=1.5,cex.lab=1.5,cex.main=1.8,pch=0)
	sid = which(names(nm) == sample_id);
	
	if(mother_id != "") { sid = c(sid,which(names(nm) == mother_id)); colour = c(colour,colors()[547]); }
	if(father_id != "") { sid = c(sid,which(names(nm) == father_id)); colour = c(colour,"green"); }
	for(si in 1:length(sid)) { lines(nm[,1], nm[,sid[si]], lty=1,col=colour[si]); }
	
	if(cnv=="DUP") {
		legend("bottomleft",legend=lgnd1, text.col="darkviolet",cex=1.7)
	} else {
		legend("topleft",legend=lgnd1, text.col="darkviolet",cex=1.7)
	}
	abline(v=cnvrange, col=9, lty=3)
	abline(v=cnvmaxrange, col=9, lty=2)
	
	nm1 = all100score[,c(3,5:length(all100score))]; nm2 = all100score[,4:length(all100score)]; names(nm1)[1] = "chrpos"; names(nm2)[1] = "chrpos"; nm = rbind(nm1,nm2); nm = nm[with(nm,order(nm[,1])),];
	
	par(mar = c(0,4,1,1))
	matplot(nm[,1], nm[,2:length(nm)], type="l", lty=1, xlab="", ylab="ADM3 scores", col=8, xaxt="n",cex.axis=1.5,cex.lab=1.5,cex.main=1.8,pch=0)
	
	for(si in 1:length(sid)) { lines(nm[,1], nm[,sid[si]], lty=1,col=colour[si]); }
	abline(v=cnvrange, col=9, lty=3)
	abline(v=cnvmaxrange, col=9, lty=2)
	
	par(mar = c(4,4,1,1))
	lv = sample(0:10,length(nm[,1]), replace=TRUE); lv[1] = 0; lv[2] = 10;
	plot(nm[,1],lv, xlab="Chr position", ylab="", yaxt="n",col="#FFFFFF", cex.axis=1.5,cex.lab=1.5)
	abline(v=cnvrange, col=9, lty=3)
	abline(v=cnvmaxrange, col=9, lty=2)
	if(nrow(gene_info)>0) {
		# Genes + Transcripts
		GT = gene_info[as.numeric(rownames(gene_info)) %in% all100call[,1],c("gene_symbol")]
		genes = getGeneSymbols(GT);
		gx = cnvrange[1]; gy = min(colmat[3,])-1;
		text(x=gx, y=gy, labels=genes, col=colors()[556], pos=4, cex=1.5)
	}
	known_rstr_df = data.frame();
	if(known_rstr == "") {
		known_rstr_df = getKnownCNVString(Fd_rep);
	}
	
	if(nrow(known_rstr_df) > 0) {
		for(cm in colmat[1,]) {
			cnve = known_rstr_df[known_rstr_df[,3]==cm,]; ccol = colors()[colmat[2,cm]]; yloc = as.numeric(colmat[3,cm]); 
			if(nrow(cnve)>0) { 
				for(jj in 1:nrow(cnve)) {
					ks = c(cnve[jj,1]:cnve[jj,2])
					lines(x=ks, y=rep(colmat[3,cm],length(ks)), col=ccol,lwd=2)		
				}
			}
		}
	}
	legend("bottomright",legend=c(misc_str,paste("Mother",mother_id),paste("Father",father_id),"Known_Dels","Known_Dups"), box.lty=0, lty = 1, col = colors()[c(3,547,254,colmat[2,1:2])], lwd=2, cex=1.3)
	dev.off()
}

