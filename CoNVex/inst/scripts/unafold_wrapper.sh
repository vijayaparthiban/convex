#!/usr/local/bin/bash
# fires unafold - multiple processes 

# * $1 Prefix for the farm output file (dG, Tm values are stored here) – option 1
# * $2 Prefix for RNA fasta file – option 2
# * $3 Prefix for RevComp fasta file – option 3
# * $4 How many times (as per the number of input fasta files) should it run? – option 4

total=`expr $4 + 1`
count=1			

# Enable this to include UNAfold in your path. You may like to do it in a different way!
#export PATH=$PATH:/path/to/unafold/bin:

while [ $count -lt $total ]	

do
# Sanger additional catalog
#bsub -J UNAfoldA$count -o rnae_unafold_output_sanger$count.txt perl melt.pl --NA=RNA --temperature=65 --Ct=0.00001 exon_sanger_addl_rna$count.fasta exon_sanger_addl_revcomp$count.fasta

# Agilent Catalog
echo "perl melt.pl --NA=RNA --temperature=65 --Ct=0.00001 $2$count.fa $3$count.fa > $1$count.txt"
count=`expr $count + 1`
if [[ "$5" && $5 == 'wait' ]]; then
echo "wait" 
fi 
done	


