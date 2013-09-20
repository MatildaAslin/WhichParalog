#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Slurp;
use Bio::TreeIO;
use IO::String;
use Bio::AlignIO;
use Bio::Align::ProteinStatistics;
use Bio::Tree::DistanceFactory;
use Bio::Tree::TreeI;
use Bio::SimpleAlign;

#Declaring variables
my $input;
my $seq;
my @name; #get the name from the seq
my @all; #all names

#Saving input parameters
GetOptions ("i|input=s" => \$input);

#Load alignment file
my $str = Bio::AlignIO->new(-file =>$input, -format => 'fasta');
my $aln = $str->next_aln();

$aln->sort_alphabetically;

my $alnio = new Bio::AlignIO(-file => '>sortedAln.fasta', -format => 'fasta');
$alnio->write_aln($aln);
