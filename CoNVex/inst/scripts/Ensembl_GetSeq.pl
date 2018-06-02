#!/usr/local/bin/perl

### MODIFY THE LINES BELOW TO SPECIFY YOUR LOCAL ENSEMBL API LOCATION ###
use lib '/software/pubseq/PerlModules/Ensembl/www_70_1/ensembl/modules/';
use lib '/software/pubseq/PerlModules/Ensembl/www_70_1/bioperl-live/';

### USAGE ###
# perl Ensembl_GetSeq.pl <genomic_regions_input_file>

# GET THE DB ADAPTOR - THIS IS FOR HG19 BUILD
use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';
$registry->load_registry_from_db(
	-host => 'ensembldb.ensembl.org',
	-user => 'anonymous'
);
my @db_adaptors   = @{ $registry->get_all_DBAdaptors() };
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

# INPUT GENOMIC REGIONS FILE WITH 3 COLUMNS: CHR, START, END [CHR MUST BE IN THE FORMAT 1,2,...X,Y] - NO HEADER
my $regions_file = $ARGV[0];
my $fasta_file   = $regions_file . ".fa";

# OUTPUT FASTA FILE - ADD .fa EXTENSION TO THE INPUT FILE
open( INFILE, "<$regions_file" ) || die("Could not open genomic regions file!");
open( OUTFILE, ">$fasta_file" ) || die("Could not create output fasta file!");

my $line;
my $slice;
my $unmasked_DNAL;
my $args;
my $start;
my $end;

while (my $line = <INFILE>) {
	chomp($line);
	my @args = split( '\t', $line );

	if ( @args[0] == 23 ) {
		@args[0] = "X";
	}
	$start = @args[1];
	$end   = @args[2];

	$slice = $slice_adaptor->fetch_by_region( 'chromosome', @args[0], $start, $end );
	$unmasked_DNAL = $slice->seq();
	print OUTFILE ">", @args[0], ":", $start, ":", $end, "\n", $unmasked_DNAL,"\n";
}
close(INFILE)  || die("Could not close genomic regions file!");
close(OUTFILE) || die("Could not create/close fasta file!");

print "Fasta output file:\n";
print $fasta_file,"\n";
