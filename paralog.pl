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
my $leaf;

#Saving input parameters
GetOptions ("i|input=s" => \$input);

#Load tree file
my $str = Bio::TreeIO->new(-file =>$input, -format => 'nexus');
my $tree = $str->next_tree();

my @leaves = $tree->get_leaf_nodes();

#Loops over all leaves
foreach (@leaves){
	
	#Converts $seq to string and prints
	print $_->id . "\n"
}
