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
my %groups;
my $speciesName;
my %leavesName;
my @paralogSpecies;

#Saving input parameters
GetOptions ("t|treefile=s" => \$treefile, "g|groupfile=s" => \$groupfile);

open(INFILE, $groupfile) or die "Can't open file: $!\n";

#Make hash out of groupfile
foreach my $line (<INFILE>) { 

	chomp($line);
		
	my ($groupname, $protein) = split("\t", $line);
	
	my @species = split("_", $protein);

	#$groups{$species[0] . "_" . $species[1]} = $groupname;

	$groups{$species[0]} = $groupname;
}

#Load tree file
my $str = Bio::TreeIO->new(-file =>$treefile, -format => 'nexus');
my $tree = $str->next_tree();

#Puts all leaves in an array
my @leaves = $tree->get_leaf_nodes();

#Loops over all leaves and puts the speciesnames in hash
foreach (@leaves){
	my $proteinID = $_->id;
	my @id = split("_", $proteinID);
	$speciesName = $id[0] . "_" . $id[1];
	$leavesName{$speciesName}++;
}

#Check if there are >1 of one species in tree
for (keys %leavesName){
	if($leavesName{$_} > 1){
		push(@paralogSpecies,$_);
	}
}

for(@paralogSpecies){
	print $_ . "\n";
}
