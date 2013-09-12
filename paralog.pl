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
use Bio::Tree::TreeFunctionsI;
use Bio::SimpleAlign;

#Declaring variables
my $treefile;
my $groupfile;
my $leaf;
my %groups;
my $speciesName;
my %leavesName;
my @paralogSpecies;
my @bacteria;
my @reBacteria;

#Saving input parameters
GetOptions ("t|treefile=s" => \$treefile, "g|groupfile=s" => \$groupfile);

open(INFILE, $groupfile) or die "Can't open file: $!\n";

#Make hash out of groupfile
foreach my $line (<INFILE>) { 

	chomp($line);
		
	my ($groupname, $protein) = split("\t", $line);
		
	my @species = split("_", $protein);

	$groups{$species[0]} = $groupname;
}

#Load tree file
my $str = Bio::TreeIO->new(-file =>$treefile, -format => 'nexus');
my $tree = $str->next_tree();

#Puts all leaves in an array
my @leaves = $tree->get_leaf_nodes();

for (@leaves){
	if( checkIfBacteria($_) ){
		push(@bacteria, $_);
	}
}

#re-root tree with bacteria (the first one in the array)
$tree->reroot($bacteria[0]); 

#get leaves of re-rooted tree
my @reLeaves = $tree->get_leaf_nodes();

#gets bacterial nodes of re-rooted tree
for (@reLeaves){
	if( checkIfBacteria($_) ){
		push(@reBacteria, $_);
	}
}

my $lca = $tree->get_lca(-nodes => \@reBacteria);

print $lca . "\n";
	

#Loops over all leaves and puts the speciesnames in hash
foreach (@leaves){
	
	my $proteinID = $_->id;
	
	my @id = split("_", $proteinID) ;
	
	#Checks that all species names coorrespond to the ones in the groupfile
	typoLocator($id[0]) or die "ERROR: Typo (" . $id[0] . ") in treefile: " . $treefile . "\n Please correct this before running Tree Doctor\n" ;
	
	$speciesName = $id[0] . "_" . $id[1];
	
	$leavesName{$speciesName}++;
}

#Check if there are >1 of one species in tree
for (keys %leavesName){
	if($leavesName{$_} > 1){
		push(@paralogSpecies,$_);
	}
}

print "Paralog species:\n";
for(@paralogSpecies){
	print $_ . "\n";
}

#For debugging
#foreach my $key (keys %leavesName){
#	my @damp = split("_", $key); 
#	print $damp[0] . " => " . $groups{$damp[0]} . "\n";
#}

sub checkIfBacteria{
	my $nodeName = $_[0]->id;
	my @splitSpecies = split("_", $nodeName);
	
	if ($groups{$splitSpecies[0]} eq "bact"){
		return 1;
	}
}

sub typoLocator{
	
	if ($groups{$_[0]}){
		return 1;
	}
}
