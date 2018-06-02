# TODO: SW CNV call function
# 
# Author: pv1
###############################################################################

`SWCNV` <- function(p,tdel,tdup,dv,GAMfile,sample_id,centromere_regions_file,output_folder,sw_exec) {
	
	d = read.table(GAMfile); d = data.frame(seq(1,nrow(d)),d); # Read the sample EW scores

	# get Centromere Boundaries - based on the centromere regions file
	hg_chr_gaps = read.table(file=centromere_regions_file, header=TRUE,comment.char = "");
	hg_cent = hg_chr_gaps[hg_chr_gaps$type=='centromere',];
	hg_cent = data.frame(substr(hg_cent[,2],4,5),hg_cent[,3],hg_cent[,4]); names(hg_cent) = c("chr","start","end");
	centro_boundaries = GetCentroBoundaries(d[,2:4],hg_cent);
	
	# Set the output folder
	setwd(output_folder); # GAMfiles folder
	
	#print(paste("Applying Smith-Waterman algorithm on sample ",sample_id, "with p =",p, " ....")) # Print message
	u = unique(d[,2]); del_report = NULL; dup_report = NULL;
	
	for (c in 1:length(u)) {
		
		SE = centro_boundaries[centro_boundaries[,1]==u[c],]; 
		for(cs in 1:nrow(SE)) { # left and right sides of centromeres separately
			defaults = paste("-s",format(SE[cs,2],scientific=FALSE),"-e",format(SE[cs,3],scientific=FALSE),"-r 8 -p",p) # Set SW-Array command line options
			sw_command_dup = paste(sw_exec,"-f",GAMfile,"-c",u[c],defaults,"-g 1"); sw_command_del = paste(sw_exec,"-f",GAMfile,"-c",u[c],defaults,"-g -1");
			tmp_dup = read.table(pipe(sw_command_dup), header=TRUE); tmp_del = read.table(pipe(sw_command_del), header=TRUE);
			
			if(nrow(tmp_dup) > 0) {
				tmp_dup = tmp_dup[order(tmp_dup[,2],tmp_dup[,3]),]; dup_report = data.frame(rbind(dup_report,tmp_dup));
			}
			if(nrow(tmp_del) > 0) {
				tmp_del = tmp_del[order(tmp_del[,2],tmp_del[,3]),]; del_report = data.frame(rbind(del_report,tmp_del));
			}
		}	
	} # For each chromosome
	
	Fdup_report = dup_report[((dup_report[,5]/(dup_report[,4])^dv)) >= tdup,]; Fdel_report = del_report[((del_report[,5]/(del_report[,4])^dv)) >= tdel,];
	Ndups = nrow(Fdup_report); Ndels = nrow(Fdel_report); Nall = Ndups+Ndels;
	
	if(Nall==0) {
		write(paste("Couldn't detect CNVs in this sample:",sample_id), stderr());
		write("Please take a look at the sample data if this is an error (noise, quality, etc.)", stderr());
		return(data.frame());
	} else if(Ndels==0) {
		write(paste("Couldn't detect deletions in this sample:",sample_id), stderr());
		write("Please take a look at the sample data if this is an error (noise, quality, etc.)", stderr());
		write("Returning only duplications...", stderr());
		Fdup_reportX = GetCNVStartEnd(d,Fdup_report);
		Fdup_CVscores = Fdup_reportX[,5]/(Fdup_reportX[,4]^dv); 
		Fdup_CVscores = round(Fdup_CVscores,digits=2);
		Fdup_reportX = data.frame(Fdup_reportX[,1:4],Fdup_CVscores,"DUP",sample_id); names(Fdup_reportX) = c("chr","start","end","num_probes","convex_score","cnv_type","sample_id");
		return(Fdup_reportX);
	} else if(Ndups==0) {
		write(paste("Couldn't detect duplications in this sample:",sample_id), stderr());
		write("Please take a look at the sample data if this is an error (noise, quality, etc.)", stderr());
		write("Returning only deletions...", stderr());
		Fdel_reportX = GetCNVStartEnd(d,Fdel_report);	
		Fdel_CVscores = Fdel_reportX[,5]/(Fdel_reportX[,4]^dv);
		Fdel_CVscores = round(Fdel_CVscores,digits=2);
		Fdel_reportX = data.frame(Fdel_reportX[,1:4],Fdel_CVscores,"DEL",sample_id); names(Fdel_reportX) = c("chr","start","end","num_probes","convex_score","cnv_type","sample_id");	
		return(Fdel_reportX);
	}
	else {
		Fdup_reportX = GetCNVStartEnd(d,Fdup_report);
		Fdel_reportX = GetCNVStartEnd(d,Fdel_report);	
		Fdup_CVscores = Fdup_reportX[,5]/(Fdup_reportX[,4]^dv); Fdel_CVscores = Fdel_reportX[,5]/(Fdel_reportX[,4]^dv);
		Fdup_CVscores = round(Fdup_CVscores,digits=2); Fdel_CVscores = round(Fdel_CVscores,digits=2);
		Fdup_reportX = data.frame(Fdup_reportX[,1:4],Fdup_CVscores,"DUP",sample_id); names(Fdup_reportX) = c("chr","start","end","num_probes","convex_score","cnv_type","sample_id");
		Fdel_reportX = data.frame(Fdel_reportX[,1:4],Fdel_CVscores,"DEL",sample_id); names(Fdel_reportX) = names(Fdup_reportX);	
		return(rbind(Fdup_reportX,Fdel_reportX));
	}
}

