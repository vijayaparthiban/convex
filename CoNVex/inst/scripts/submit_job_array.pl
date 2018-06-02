#!/usr/bin/env perl
use strict;
use warnings;
#use version;

use Carp qw(croak);
use English qw(-no_match_vars);
use Fatal qw(:void open close);

# Documentation is at the __END__

# $Id: submission.pl,v 1.3 2008/03/25 14:21:24 io1 Exp $
# $Date: 2008/03/25 14:21:24 $
# $Author: io1 $
# $Revision: 1.3 $

our $VERSION = ( '$Revision:$' =~ m/ (\d+.\d+) /xms )[0];

# Applications implementation goes here ...

# Get the command_file and LSB_JOBINDEX
my $current_index = $ENV{'LSB_JOBINDEX'};

# Get the commands
# This modification allows us to pass the commands fromm a pipe
my @commands = <>;

# Run the commands
system $commands[ $current_index - 1 ];

exit 0;

__END__

=head1 NAME

submit_job_array - submit jobs to the farm


=head1 VERSION

This document describes submit_job_array version 0.0.1


=head1 USAGE

bsub -J"jobname[1-<job_count>]%<concurrent_jobs>" -q <queue> -o <out_file> submit_job_array.pl <command_file>

  
=head1 DESCRIPTION

A simple wrapper script for submitting jobs to the farm as job arrays. See the
LSF documentation at L<< http://www.wtgc.org/IT/ISG/courses/FarmCourse2007.html >>
and L<< http://www.wtgc.org/IT/ISG/lsf/job-arrays.shtml >> for more details.


=head1 REQUIRED ARGUMENTS 

=over 4

=item C<< command_file >>

A file containing all the commands in the job array, one command per line.

=back


=head1 OPTIONS

None.


=head1 DIAGNOSTICS

=over

=item C<< Please provide a command file at %s line 23 >>

You have forgotten to provide a command file. Supply one and try again.

=back


=head1 EXIT STATUS

Not applicable.


=head1 CONFIGURATION AND ENVIRONMENT

submit_job_array requires the LSB_JOBINDEX environment variable to be set. See
LSF documentation for full details.


=head1 DEPENDENCIES

This application requires the CPAN module version.pm


=head1 INCOMPATIBILITIES

None reported.


=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to the author at C<< <<io1@sanger.ac.uk>> >>.


=head1 AUTHOR

Nelo Onyiah  C<< <<io1@sanger.ac.uk>> >>


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2008, Nelo Onyiah C<< <<io1@sanger.ac.uk>> >>. All rights reserved.

This application is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.
