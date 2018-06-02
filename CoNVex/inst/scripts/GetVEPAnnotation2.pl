#!/usr/local/bin/perl
# Author: Parthiban Vijayarangakannan
# Date: 14.02.2013
# INPUT: MUST BE similar to CoNVex output: chr, start, end should in the first 3 columns of the input file; DEL/DUP in column 6 - should end with .txt extension / should be without HEADER
# As of now, this script uses the Ensembl build 70 at Sanger, but please edit the script to change it to [ YOUR LOCAL ENSEMBL API LOCATION/VERSION ] if needed
# This script outputs only the MOST SEVERE VEP consequence. The consequences are RANKED in the following order:
# 1.	transcript_ablation
# 2.	frameshift_variant
# 3.	coding_sequence_variant
# 4.	initiator_codon_variant
# 5.	inframe_deletion
# 6.	stop_lost
# 7.	transcript_amplification
# 8.	mature_miRNA_variant
# 9.	non_coding_exon_variant
# 10.	5_prime_UTR_variant
# 11.	3_prime_UTR_variant

use lib '/software/pubseq/PerlModules/Ensembl/www_70_1/ensembl/modules/'; # CHANGE THIS IF REQUIRED
use lib '/software/pubseq/PerlModules/Ensembl/www_70_1/ensembl-variation/modules/';

my $infile = $ARGV[0];
my $ext = ".";
my $outfile = substr $infile, 0, rindex($infile, $ext);
$outfile .= ".VEP.txt";

use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous'
);
my $slice_adaptor = $registry->get_adaptor( 'human', 'core', 'slice' );
my $svf_adaptor = Bio::EnsEMBL::Registry->get_adaptor("human","variation","structuralvariationfeature");

open(INFILE, "<$infile") or die("Couldn't read input file!\nUsage: perl GetVEPAnnotation2.pl /path/to/input_file.txt"); # INFILE as in CoNVex CNV calls output
open(OUTFILE, ">$outfile"); 

my ($chr, $start, $end, $ctype);
my $vname = "my_sv";
my %cnv_type = (DEL => "deletion", DUP => "duplication", TDUP => "tandem_duplication", INS => "insertion");

my @cons_rank = ("transcript_ablation","frameshift_variant","coding_sequence_variant","initiator_codon_variant","inframe_deletion","stop_lost","transcript_amplification","mature_miRNA_variant","non_coding_exon_variant","5_prime_UTR_variant","3_prime_UTR_variant","stop_gained","regulatory_region_variant","intron_variant","upstream_gene_variant","downstream_gene_variant","intergenic_variant","unknown_consequence");
my %cons_rank_hash;
@cons_rank_hash{@cons_rank} = (0..$#cons_rank);

my ($msocons, $sot, $msostr);
my $count = 0;

while(my $line = <INFILE>) {
	
	chomp($line);
	my @args = split('\t',$line);
	$chr = @args[0]; $start = @args[1]; $end = @args[2]; $ctype = @args[5];
	
	if ($chr==23)  {
		$chr="X";
	}
	elsif ($chr==24)  {
		$chr="Y";
	}
	my $slice = $slice_adaptor->fetch_by_region("chromosome", $chr);

	my $svf = Bio::EnsEMBL::Variation::StructuralVariationFeature->new_fast({
    	start          => $start,
    	end            => $end,
	    strand         => 1,
	    adaptor        => $svf_adaptor,
	    variation_name => $vname,
	    slice          => $slice,
	    class_SO_term  => $cnv_type{$ctype},    # possible values: 'insertion','deletion','tandem_duplication','duplication'
	});

	# get the list of the consequence types
	my $overlap_cons = $svf->get_all_StructuralVariationOverlaps;
	my %msocons_list; $msostr = "";
	my $mso_rank_tmp = 100, $mso_rank = 17;
	
	foreach my $con(@$overlap_cons) {
	    $msocons = $con->display_consequence(); 
	    # print "$chr:$start-$end >> $sot; Rank: $rank; Tier: $tier\n"; 
	    $msocons_list{$msocons}++;
	    $mso_rank_tmp = $cons_rank_hash{$msocons};
	    if($mso_rank_tmp < $mso_rank) {
	    	$mso_rank = $mso_rank_tmp
	    }
	}
	#for my $msk (keys %msocons_list) {
	#	$msostr .= "$msk ($msocons_list{$msk}); ";
	#}
	print OUTFILE "$line\t$cons_rank[$mso_rank]\n";
}
close(INFILE); close(OUTFILE);

print "Ensembl VEP consequences for the given CNV list:\n";
print "$outfile\n";

