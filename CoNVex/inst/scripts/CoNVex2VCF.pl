#!/usr/bin/env perl
# CoNVex (tab-delimited) to VCF format converter EXAMPLE SCRIPT
# Please INLCUDE/EXCLUDE details in %INFO to match that of columns in tab-delimited txt file with header
# By default, ALL %INFO columns are expected to be present. If not, VCF format will be incorrect
# If a column is not present, VCF headers (including that of absent columns) will be created

use strict;
use warnings FATAL => 'all';

use Getopt::Long;
use Pod::Usage;
use Iterator::Simple::Util::CSV qw( icsv );

use DateTime;
use Const::Fast;
use Hash::MoreUtils qw( slice_def );
use Vcf;
use FaSlice;

const my $FILE_FORMAT => 'VCFv4.1';
const my @VCF_COLUMNS => qw( CHROM POS ID REF ALT QUAL FILTER INFO FORMAT );

my %INFO = (
	'MEANL2R' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Mean log2 ratio of cnv interval',
        accessor    =>  sub { $_[0]->{mean_l2r}; }
    },
    'MADL2R' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'MAD of log2 ratio across samples at cnv interval',
        accessor    =>  sub { $_[0]->{mad_l2r}; }
    },
    'SVTYPE' => {
        Type        => 'String',
        Number      => 1,
        Description => 'Type of structural variant',
        accessor    => sub { $_[0]->{cnv_type}; }
    },
    'SVLEN' => {
        Type        => 'Integer',
        Number      => 1,
        Description => 'Length of structural variant',
        accessor    => sub { $_[0]->{end}+0 - $_[0]->{start}+0;}
    },
    'END' => {
        Type        => 'Integer',
        Number      => 1,
        Description => 'End position of the variant described in this record',
        accessor    => sub { $_[0]->{end}; }
    },
    'NUMPROBES' => {
        Type        => 'Integer',
        Number      => 1,
        Description => 'Number of probes in CNV interval',
        accessor    => sub { $_[0]->{num_probes}; }
    },
    'CONVEXSCORE' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Call quality score',
        accessor    =>  sub { $_[0]->{convex_score}; }
    },
    'COMMONFORWARDS' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Forwards overlap with CNV consensus common events',
        accessor    =>  sub { $_[0]->{common_forward}; }
    },
     'COMMONBACKWARDS' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Backwards overlap with CNV consensus common events',
        accessor    =>  sub { $_[0]->{common_backward}; }
    },
     'RAREFORWARDS' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Forwards overlap with CNV consensus rare events',
        accessor    =>  sub { $_[0]->{rare_forward}; }
    },
     'RAREBACKWARDS' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Backwards overlap with CNV consensus rare events',
        accessor    =>  sub { $_[0]->{rare_backward}; }
    },
     'SEGDUPFORWARDS' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Forwards overlap with segmental duplications',
        accessor    =>  sub { $_[0]->{segdup_forward}; }
    },
     'SEGDUPBACKWARDS' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Backwards overlap with segmental duplications',
        accessor    =>  sub { $_[0]->{segdup_backward}; }
    },
     'INTERNALFREQ' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Internal frequency',
        accessor    =>  sub { $_[0]->{internal_freq}; }
    },
    'RC50INTERNALFREQ' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'RC50 Internal frequency',
        accessor    =>  sub { $_[0]->{rc50_internal_freq}; }
    },
    'CALLQC' => {
        Type        => 'String',
        Number      => 1,
        Description => 'Call QC',
        accessor    =>  sub { $_[0]->{call_qc}; }
    },
    'VEP' => {
        Type        => 'String',
        Number      => 1,
        Description => 'VEP consequenes, Ensembl 73',
        accessor    =>  sub { $_[0]->{vep_annotation}; }
    },
    'GENESYMBOL' => {
        Type        => 'String',
        Number      => '.',
        Description => 'Gene symbol (Ensembl 73)',
        accessor    =>  sub { my $gs = $_[0]->{gene_symbol}; $gs =~ s/;$//; $gs =~ s/;/,/g; return $gs; }
    },
    'CONVEX' => {
        Type        => 'Flag',
        Number      => 0,
        Description => 'CNV called by CoNVeX',
        accessor    => sub { 1 }
    }
);

my %FORMAT = (
    'INHERITANCE' => {
        Type        => 'String',
        Number      => 1,
        Description => 'Inheritance configuration',
        accessor    =>  sub { "." }
    },
    'INHERITANCEP' => {
        Type        => 'Float',
        Number      => 1,
        Description => 'Posterior probability inheritance configuration',
        accessor    =>  sub { "." }
    }
);


sub write_column_header {
    my $sample_id = shift;
    print "#" . join( "\t", @VCF_COLUMNS, $sample_id ) . "\n";

    return;
}

sub escape_meta {
    my $v = shift;

    if ( $v !~ m/[,="'\s]/ ) {
        return $v;
    }

    $v =~ s{\\}{\\\\}g;
    $v =~ s{"}{\\"}g;

    return sprintf( '"%s"', $v );
}

sub add_headers {
    my ( $vcf, $sample_id, $run_date, $pipeline_version ) = @_;

    $vcf->add_header_line( { key => 'fileformat', value => $FILE_FORMAT } );
    $vcf->add_header_line( { key => 'rundate',    value => $run_date->ymd } );

    if ( defined $pipeline_version ) {
        $vcf->add_header_line( { key => 'source', value =>  "CNV calls from CoNVex version $pipeline_version" } );
    }
    else {
        $vcf->add_header_line( { key => 'source', value => 'CNV calls from CoNVex' } );
    }

    for my $id ( sort keys %INFO ) {
        $vcf->add_header_line( { key => 'INFO', ID => $id, slice_def $INFO{$id}, qw( Number Type Description ) } );
    }

    for my $id ( sort keys %FORMAT ) {
        $vcf->add_header_line( { key =>  'FORMAT', ID => $id, slice_def $FORMAT{$id}, qw( Number Type Description ) } );
    }
    
    $vcf->add_header_line( { key =>  'FILTER', ID => 'FAIL', Description => 'Fail call QC' } );

    # Dirty hack, abusing Vcf class encapsulation
    $vcf->{columns} = [ @{ $vcf->{mandatory} }, 'FORMAT', $sample_id ];
}

sub format_info {
    my ( $key, $variant ) = @_;

    my $meta = $INFO{$key};

    return unless $meta and $meta->{accessor};

    my $v = $meta->{accessor}->( $variant );

    return unless defined $v;

    if ( $meta->{Type} eq 'Flag' ) {
        return $v ? $key : ();
    }

    return sprintf( '%s=%s', $key, $v );
}

sub variant_info {
    my $variant = shift;

    my @info;
    for my $key ( keys %INFO ) {
        my $f = $INFO{$key}{accessor}
            or next;
        if ( my $v = format_info( $key, $variant ) ) {
            push @info, $v;
        }
    }
    return join ';', @info;
}


sub variant_format {
    my $variant = shift;

    my (@keys, @vals);

    while ( my ( $key, $meta ) = each %FORMAT ) {
        my $f = $meta->{accessor}
            or next;
        defined( my $v = $f->($variant) )
            or next;
        push @keys, $key;
        push @vals, $v;
    }

    return ( join( ':', @keys ), join( ':', @vals ) );
}


sub variant_id {
    my $variant = shift;

    my $id = $variant->id || '.';

    $id =~ s/,/;/g;

    return $id;
}

sub write_variant {
    my ($variant, $ref) = @_;

    #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT

    print join( "\t",
                $variant->{'chr'},
                $variant->{'start'},
                '.',
                $ref->get_base($variant->{'chr'}, $variant->{'start'}),
                '<' . $variant->{'cnv_type'} . '>',
                $variant->{convex_score},
                $variant->{call_qc},
                variant_info( $variant ),
                variant_format( $variant )
            ) . "\n";
}

GetOptions(
    'help'        => sub { pod2usage( -verbose => 1 ) },
    'man'         => sub { pod2usage( -verbose => 2 ) },
    'sample_id=s' => \my $sample_id,
    'convex_calls_file=s' => \my $file_name,
    'ref=s' => \my $reference
) or pod2usage( 2 );

pod2usage( "--sample-id & --file-name & --ref are required" ) unless defined $sample_id && $file_name && $reference;

my $convex_version = "0.5+"; # CHANGE IF REQUIRED
my $run_date = DateTime->from_epoch( epoch => (stat $file_name)[9] );

my $vcf = Vcf->new();
add_headers( $vcf, $sample_id, $run_date, $convex_version);
print $vcf->format_header();

my $convex_results =icsv( $file_name, use_header=>1, sep_char => "\t" );
my $ref = FaSlice->new(file=>$reference);

# $vcf (Vcf.pm) is explicitly not used for this below - as the proper format seems to have issues with VEP (for instance)
while ( my $result = $convex_results->next ) { 
    write_variant( $result, $ref );
};

__END__

=pod

=head1 NAME

CoNVex2VCF

=head1 SYNOPSIS

  CoNVex2VCF --sample_id=MYSAMP_12345 --convex_calls_file=MYSAMP_12345_CoNVexCNVs.txt

=head1 DESCRIPTION

Generate a VCF containing the convex results for the
specified sample. --convex_calls_file requires calls with all QC metrics with the following column headers.

QC metrics' column headers:
"chr","start","end","num_probes","convex_score","cnv_type""convex_score",
"mean_l2r","mad_l2r","internal_freq","rc50_internal_freq","common_forward",
"common_backward","call_qc","vep_annotation","gene_symbol","segdup_forward",
segdup_backward"

=cut

