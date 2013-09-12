#!/usr/bin/env perl
# $Id: psiblast.pl 2447 2013-01-25 10:57:28Z hpm $
# ======================================================================
# 
# Copyright 2008-2013 EMBL - European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# ======================================================================
# PSI-BLAST Perl client.
#
# Requires SOAP::Lite. Tested with versions 0.60, 0.69 and 0.71.
#
# See:
# http://www.ebi.ac.uk/Tools/webservices/services/psiblast
# http://www.ebi.ac.uk/Tools/webservices/clients/psiblast
# http://www.ebi.ac.uk/Tools/webservices/tutorials/soaplite
# ======================================================================
print STDERR <<EOF
=============================================================================
NB: the service used by this client was decommisioned on Monday 28th March 
2011. See http://www.ebi.ac.uk/Tools/webservices/ for replacement services.
=============================================================================
EOF
;
# ======================================================================
# WSDL URL for service
my $WSDL = 'http://www.ebi.ac.uk/Tools/webservices/wsdl/WSPSIBlast.wsdl';

# Enable Perl warnings
use strict;
use warnings;

# Load libraries
use SOAP::Lite;
use LWP::Simple;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;

# Set interval for checking status
my $checkInterval = 15;

# Output level
my $outputLevel = 1;

# Process command-line options
my $numOpts = scalar(@ARGV);
my (
	$outfile, $outformat, $help,  $async,    $polljob,
	$status,  $jobid,     $trace, $sequence, $quiet,
	$verbose, $preselseq, $checkpoint, $resultTypes
);
my ( @content, @preselseqcontent, @checkpointcontent );
my %params = (
	'async' => '1',    # Use async mode and simulate sync mode in client
);
GetOptions(
	"database|D=s"  => \$params{'database'},     # Database to search
	"matrix|m=s"    => \$params{'matrix'},       # Scoring matrix to use
	"opengap|o=i"   => \$params{'opengap'},      # Open gap penalty
	"extendgap|e=i" => \$params{'extendgap'},    # Gap extension penalty
	"exp|E=f"       => \$params{'exp'},          # E-value threshold
	"filter|F=s"    => \$params{'filter'},       # Low complexity filter
	"expmulti|x=f"  =>
	  \$params{'expmulti'},    # E-value threshold for selecting sequences
	"dropoff|X=i"      => \$params{'dropoff'},         # Dropoff score
	"finaldropoff|Z=i" => \$params{'finaldropoff'},    # Final dropoff score
	"scores|s=i"       => \$params{'scores'},          # Number of scores
	"align|A=i"        => \$params{'numal'},           # Number of alignments
	"preselseq=s"  => \$preselseq,   # Selected sequence file for next iteration
	"checkpoint=s" => \$checkpoint,  # Checkpoint file for next iteration
	"sequence=s"   => \$sequence,    # Query sequence file or DB:ID
	"outfile|O=s"  => \$outfile,     # Output file name
	"outformat=s"  => \$outformat,   # Output file type
	"help|h"       => \$help,        # Usage help
	"async|a"      => \$async,       # Asynchronous submission
	"polljob"      => \$polljob,     # Get results
	"status"       => \$status,      # Get status
	"jobid|j=s"    => \$jobid,       # JobId
	'resultTypes'  => \$resultTypes, # List result formats
	"email|S=s" => \$params{email},  # E-mail address
	'quiet'     => \$quiet,          # Decrease output level
	'verbose'   => \$verbose,        # Increase output level
	"trace"     => \$trace,          # SOAP message debug
);
if ($verbose) { $outputLevel++ }
if ($quiet)   { $outputLevel-- }

# Get the script filename for use in usage messages
my $scriptName = basename( $0, () );

# Print usage and exit if requested
if ( $help || $numOpts == 0 ) {
	&usage();
	exit(0);
}

# If required enable SOAP message trace
if ($trace) {
	print STDERR "Tracing active\n";
	SOAP::Lite->import( +trace => 'debug' );
}

# Namespace and endpoint need to be used instead of the WSDL. By default
# these are extracted from the WSDL.
my ( $serviceEndpoint, $serviceNamespace ) = &from_wsdl($WSDL);

# Create the service interface, setting the fault handler to throw exceptions
my $soap = SOAP::Lite->proxy(
	$serviceEndpoint,
	timeout => 6000,    # HTTP connection timeout
	     #proxy => ['http' => 'http://your.proxy.server/'], # HTTP proxy
  )->uri($serviceNamespace)->on_fault(
	# Map SOAP faults to Perl exceptions (i.e. die).
	sub {
		my $soap = shift;
		my $res  = shift;
		if ( ref($res) eq '' ) {
			die($res);
		}
		else {
			die( $res->faultstring );
		}
		return new SOAP::SOM;
	}
);

# Extract the service namespace and endpoint from the service WSDL document 
# for use when creating the service interface.
sub from_wsdl {
	my (@retVal) = ();
	my $wsdlStr = get($WSDL); # Get WSDL using LWP.
	if(!defined($wsdlStr) || $wsdlStr eq '') {
		die "Error: unable to get WSDL document from $WSDL";
	}
	# Extract service endpoint.
	if ( $wsdlStr =~ m/<(\w+:)?address\s+location=["']([^'"]+)['"]/ ) {
		push( @retVal, $2 );
	}
	# Extract namespace.
	if ( $wsdlStr =~
		m/<(\w+:)?definitions\s*[^>]*\s+targetNamespace=['"]([^"']+)["']/ )
	{
		push( @retVal, $2 );
	}
	return @retVal;
}

# Client-side poll
sub clientPoll($) {
	my $jobid  = shift;
	my $result = 'PENDING';

	# Check status and wait if not finished
	#print STDERR "Checking status: $jobid\n";
	while ( $result eq 'RUNNING' || $result eq 'PENDING' ) {
		my $response = $soap->checkStatus($jobid);
		$result = $response->valueof('//status');
		if ( $outputLevel > 0 ) {
			print STDERR "$result\n";
		}
		if ( $result eq 'RUNNING' || $result eq 'PENDING' ) {

			# Wait before polling again.
			sleep $checkInterval;
		}
	}
}

# Get the results for a jobid
sub getResults($) {
	my $jobid = shift;

	# Check status, and wait if not finished
	clientPoll($jobid);

	# Use JobId if output file name is not defined
	unless ( defined($outfile) ) {
		$outfile = $jobid;
	}

	# Get list of data types
	my $response = $soap->getResults($jobid);
	my $resultTypes = $response->valueof('//result');

	# Get the data and write it to a file
	if ( defined($outformat) ) {    # Specified data type
		my $selResultType;
		foreach my $resultType (@$resultTypes) {
			if ( $resultType->{type} eq $outformat ) {
				$selResultType = $resultType;
			}
		}
		if ( defined($selResultType) ) {
			$response = $soap->poll( $jobid, $selResultType->{type} );
			my $res = $response->valueof('//result');
			if ( $outfile eq '-' ) {
				write_file( $outfile, $res );
			}
			else {
				write_file( $outfile . '.' . $selResultType->{ext}, $res );
			}
		}
		else {
			die "Error: unknown result format \"$outformat\"";
		}
	}
	else {    # Data types available
		      # Write a file for each output type
		for my $resultType (@$resultTypes) {
			if ( $outputLevel > 1 ) {
				print STDERR "Getting $resultType->{type}\n";
			}
			$response = $soap->poll( $jobid, $resultType->{type} );
			my $res = $response->valueof('//result');
			if ( $outfile eq '-' ) {
				write_file( $outfile, $res );
			}
			else {
				write_file( $outfile . '.' . $resultType->{ext}, $res );
			}
		}
	}
}

# Read a file
sub read_file($) {
	my $filename = shift;
	my ( $content, $buffer );
	if ( $filename eq '-' ) {
		while ( sysread( STDIN, $buffer, 1024 ) ) {
			$content .= $buffer;
		}
	}
	else {    # File
		open( FILE, $filename )
		  or die "Error: unable to open input file $filename ($!)";
		while ( sysread( FILE, $buffer, 1024 ) ) {
			$content .= $buffer;
		}
		close(FILE);
	}
	return $content;
}

# Write a result file
sub write_file($$) {
	my ( $filename, $data ) = @_;
	if ( $outputLevel > 0 ) {
		print STDERR 'Creating result file: ' . $filename . "\n";
	}
	if ( $filename eq '-' ) {
		print STDOUT $data;
	}
	else {
		open( FILE, ">$filename" )
		  or die "Error: unable to open output file $filename ($!)";
		syswrite( FILE, $data );
		close(FILE);
	}
}

# List result formats for job
if ( defined($resultTypes) && defined($jobid) ) {
	if ( $outputLevel > 1 ) {
		print STDERR "Getting output formats for job $jobid\n";
	}
	if($outputLevel > 0) {
		print "Type\tExtension\n==============================\n";
	}
	my $response = $soap->getResults($jobid);
	my $resultTypes = $response->valueof('//result');
	foreach my $resultType (@$resultTypes) {
		print $resultType->{'type'}, "\t", $resultType->{'ext'}, "\n";
	}
}

# Print usage if bad argument combination
elsif (   !( $polljob || $status )
	&& !( defined( $ARGV[0] ) || defined($sequence) ) )
{
	print STDERR 'Error: bad option combination', "\n";
	&usage();
	exit(1);
}

# Poll job and get results
elsif ( $polljob && defined($jobid) ) {
	if ( $outputLevel > 1 ) {
		print "Getting results for job $jobid\n";
	}
	&getResults($jobid);
}

# Job status
elsif ( $status && defined($jobid) ) {
	if ( $outputLevel > 0 ) {
		print STDERR "Getting status for job $jobid\n";
	}
	my $response = $soap->checkStatus($jobid);
	my $result = $response->valueof('//status');
	print "$result\n";
	if ( $result eq 'DONE' && $outputLevel > 0 ) {
		print STDERR "To get results: $scriptName --polljob --jobid $jobid\n";
	}
}

# Submit a job
else {
	my $content;
	if ( defined( $ARGV[0] ) ) {    # Bare option
		if ( -f $ARGV[0] || $ARGV[0] eq '-' ) {    # File
			$content =
			  { type => 'sequence', content => &read_file( $ARGV[0] ) };
		}
		else {                                     # DB:ID or sequence
			$content = { type => 'sequence', content => $ARGV[0] };
		}
	}
	if ($sequence) {                               # Via --sequence
		if ( -f $sequence || $sequence eq '-' ) {    # File
			$content = { type => 'sequence', content => &read_file($sequence) };
		}
		else {                                       # DB:ID or sequence
			$content = { type => 'sequence', content => $sequence };
		}
	}
	my (@content) = ();
	push @content, $content;

	if ($preselseq) {
		if ( -f $preselseq ) {
			$content =
			  { type => 'preselseq', content => read_file($preselseq) };
		}
	}
	push @preselseqcontent, $content;

	if ($checkpoint) {
		if ( -f $checkpoint ) {
			$content =
			  { type => 'checkpoint', content => read_file($checkpoint) };
		}
	}
	push @checkpointcontent, $content;

	my $paramsData = SOAP::Data->name('params')->type( map => \%params );
	my $contentData = SOAP::Data->name('content')->value( \@content );
	my $preselseqcontentData =
	  SOAP::Data->name('preselseqcontent')->value( \@preselseqcontent );
	my $checkpointcontentData =
	  SOAP::Data->name('checkpointcontent')->value( \@checkpointcontent );

	my $response =
		  $soap->runPSIBlast( $paramsData, $contentData, $preselseqcontentData,
			$checkpointcontentData );
	$jobid = $response->valueof('//jobid');

	if ( defined($async) ) {
		print STDOUT $jobid, "\n";
		if ( $outputLevel > 0 ) {
			print STDERR
			  "To check status: $scriptName --status --jobid $jobid\n";
		}
	}
	else {    # Synchronous mode
		if ( $outputLevel > 0 ) {
			print STDERR "JobId: $jobid\n";
		}
		sleep 1;
		&getResults($jobid);
	}
}

# Print program usage
sub usage {
	print STDERR <<EOF
PSI-BLAST
==========
   
Rapid sequence database search programs utilizing the PSI-BLAST algorithm
    
For more detailed help information refer to 
http://www.ebi.ac.uk/Tools/psiblast/psiblast_help_frame.html

[Required]

  -D, --database      : str  : database to search
  seqFile             : file : query sequence ("-" for STDIN)
 
[Optional]

  -m, --matrix	      : str  : scoring matrix
  -E, --exp           : real : 0<E<= 1000. Statistical significance threshold 
                               for reporting database sequence matches.
  -x, --expmulti      : real : 0<E<= 1000. Statistical significance threshold
                               for selecting hits for next iteration
  -A, --align	      : int  : number of alignments to be reported
  -s, --scores	      : int  : number of scores to be reported
  -n, --numal	      : int  : Number of alignments
  -o, --opengap	      : int  : Gap open penalty
  -e, --extendgap     : int  : Gap extension penalty
  -s, --fliter        : str  : Low complexity filte
  -X, --dropoff       : int  : Dropoff scor
  -Z, --finaldropoff  : int  : Final dropoff score
      --preselseq     : file : selected sequences for next iteration
      --checkpoint    : file : checkpoint file for further iteration

[General]	

  -h, --help          :      : prints this help text
  -a, --async         :      : forces to make an asynchronous query
  -S, --email	      : str  : e-mail address 
      --polljob       :      : Poll for the status of a job
  -j, --jobid         : str  : jobid that was returned when an asynchronous job 
                               was submitted.
      --resultTypes   :      : list result formats for a job
  -O, --outfile       : str  : file name for results (default is jobid;
                               "-" for STDOUT)
      --outformat     : str  : result format to retrieve
      --quiet         :      : decrease output
      --verbose       :      : increase output
      --trace	      :      : show SOAP messages being interchanged
   
Synchronous job:

  The results/errors are returned as soon as the job is finished.
  Usage: $scriptName --email <your\@email> [options...] seqFile
  Returns: results as an attachment

Asynchronous job:

  Use this if you want to retrieve the results at a later time. The results 
  are stored for up to 24 hours. 	
  Usage: $scriptName --async --email <your\@email> [options...] seqFile
  Returns: jobid

  Use the jobid to query for the status of the job. If the job is finished, 
  it also returns the results/errors.
  Usage: $scriptName --polljob --jobid <jobId> [--outfile string]
  Returns: string indicating the status of the job and if applicable, results 
  as an attachment.

EOF
}
