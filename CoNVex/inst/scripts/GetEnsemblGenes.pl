
#!/usr/local/bin/perl
# Author: Parthiban Vijayarangakannan 
# Date: 14.11.2012 / Updated: 25.06.2013
# INPUT: chr, start, end should in the first 3 columns of the input file - should end with .txt extension / should be without HEADER
# OUTPUT: two columns (Genes, Transcripts) will be appended to ALL input columns - output file ends with .ENSEMBL.txt extension
# As of now, this script uses the Ensembl build 70, but please edit the script to change it to [ YOUR LOCAL ENSEMBL API VERSION ] if required

use lib '/software/pubseq/PerlModules/Ensembl/www_72_1/ensembl/modules/'; # CHANGE THIS IF REQUIRED
use lib '/software/pubseq/PerlModules/Ensembl/www_72_1/bioperl-live/';

my $infile = $ARGV[0];
my $ext = ".";
my $outfile = substr $infile, 0, rindex($infile, $ext);
$outfile .= "_GENES.txt";

use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
		-host => 'ensembldb.ensembl.org',
		-user => 'anonymous'
);

my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

# INFILE as in CoNVex CNV calls output
open(INFILE, "<$infile") or die("Couldn't read input file!\nUsage: perl GetEnsemblGenes.pl /path/to/input_file.txt"); # INFILE as in CoNVex CNV calls output
open(OUTFILE, ">$outfile");

my ($chr, $start, $end, $slice, $all_genes);
my $size = 0; my $nv = "";

while(my $line = <INFILE>) {

	chomp($line);
	my @args = split('\t',$line);
	$chr = @args[0]; $start = @args[1];	$end = @args[2];
	
	if ($chr==23)  {
		$chr="X";
	}
	elsif ($chr==24)  {
		$chr="Y";
	}	
	$slice = $slice_adaptor->fetch_by_region( 'chromosome', $chr, $start, $end );
	$all_genes = "";
	my %gene_hash  = ();
	my %trans_hash = ();

	foreach my $gene ( @{ $slice->get_all_Genes() } )
	{      # Get all gene symbols AND put them in a hash
		my $stable_id = $gene->external_name();
		if ( !$stable_id == '' ) {
			$gene_hash{$stable_id} = $nv;
		}

	}
	for my $gene_key ( sort keys %gene_hash )
	{      # Append all gene symbols to a single string, separated by ';'
		$all_genes .= $gene_key . ';';
	}

	$size = scalar( keys %gene_hash ); # If hash is empty, add 'NA' to the string -- this mean no known gene symbols associated with this region
	if ( $size == 0 ) {
		$all_genes = "NA";
	}

	#print OUTFILE "@args[0]\t$start\t$end\t$all_genes\n";    # Output file with gene symbols AND transcripts
	print OUTFILE "$line\t$all_genes\n";    # Output file with gene symbols AND transcripts
}
close(INFILE);
close(OUTFILE);

print "Ensembl genes/transcipts for the given CNV list:\n";
print "$outfile\n";
