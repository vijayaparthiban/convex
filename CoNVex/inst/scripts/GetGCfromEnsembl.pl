
#!/usr/local/bin/perl
# Author: Parthiban Vijayarangakannan
# Date: 14.02.2012 / Updated: 14.01.2013
# INPUT: chr, start, end should in the first 3 columns of the input file - should end with .txt extension / should be without HEADER
# OUTPUT: GC will be appended to the first 3 input columns  - output file ends with .GC.txt extension with 4 columns
# As of now, this script uses Emsembl build 69, but please edit the script to change it to [ YOUR LOCAL ENSEMBL API VERSION ] if required

use lib '/software/pubseq/PerlModules/Ensembl/www_69_1/ensembl/modules/'; # CHANGE THIS IF REQUIRED

use Bio::EnsEMBL::Registry;
my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
		-host => 'ensembldb.ensembl.org',
		-user => 'anonymous'
);

my @db_adaptors = @{ $registry->get_all_DBAdaptors() };
my $slice_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Slice' );

my $infile = $ARGV[0];
my $ext = ".";
my $outfile = substr $infile, 0, rindex($infile, $ext);
$outfile .= ".GC.txt";

open(INFILE, "<$infile");
open(OUTFILE, ">$outfile");

my ($slice, $chr, $start, $end, $regend, $args, $unmasked_DNAL, $GC_count, $GC_fraction);

while(my $line = <INFILE>) {

	chomp($line);	
	my @args = split('\t',$line);
	$chr = @args[0]; $start = @args[1]; $end = @args[2];
	
	if ($chr==23)  {
		$chr="X";
	}
	elsif ($chr==24)  {
		$chr="Y";
	}
	$slice = $slice_adaptor->fetch_by_region('chromosome', @args[0], $start, $end);
	$unmasked_DNAL = $slice -> seq();
	$GC_count = ($unmasked_DNAL =~ tr/gcGC//); # Count GC
	$GC_fraction = sprintf("%.2f",($GC_count / ($end - $start + 1)));
	print OUTFILE "@args[0]\t$start\t$end\t$GC_fraction\n"; 
}
close(INFILE);
close(OUTFILE);

print "Input regions (chr,start,end in tab-delim format):\n";
print "$infile\n";
print "Output file (chr,start,end,GC):\n";
print "$outfile\n";
#bsub -q normal -R'select[mem>500] rusage[mem=500]' -M500000 -o gcc.out 'perl /path/to/GetGCfromEnsembl.pl /path/to/input_chr_regions.txt'
