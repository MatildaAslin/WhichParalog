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
my $treefile;
my $groupfile;
my $leaf;
my %group;

#Saving input parameters
GetOptions ("t|treefile=s" => \$treefile, "g|groupfile=s" => \$groupfile);

open(INFILE, $groupfile) or die "Can't open file: $!\n";

foreach my $line (<INFILE>) {
	print $line;
}

#Make hash out of groupfile


#Load tree file
my $str = Bio::TreeIO->new(-file =>$treefile, -format => 'nexus');
my $tree = $str->next_tree();

my @leaves = $tree->get_leaf_nodes();

#Loops over all leaves
foreach (@leaves){
	
	#Converts $seq to string and prints
	#print $_->id . "\n"
}
