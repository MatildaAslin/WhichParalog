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

#Make hash out of groupfile (Key: Species, Value: Group) 
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


#Put all bacterial nodes in an array
for (@leaves){
	if( checkIfBacteria($_) ){
		push(@bacteria, $_);
	}
}

#Number of bacteria in the tree
my $noBacteria = @bacteria;

#Checks if there are 3 or more bacteria in the tree
if ($noBacteria < 3){
	die "ERROR: Tree has less than 3 bacteria.\n";
}

#Checks if the bacterial species in the tree are monophyletic 
if (checkIfMonophyletic(\@bacteria)){
	#re-root tree with bacteria (the first one in the array)
	$tree->reroot($bacteria[0]); 
}
#Temporary else-statement
else {
	die "Bacteria in tree are not monophyletic!";
}

#get leaves of re-rooted tree
my @reLeaves = $tree->get_leaf_nodes();

print "Number of leaves after re-rooting: " . @reLeaves . "\n";

#gets bacterial nodes of re-rooted tree
for (@reLeaves){
	if( checkIfBacteria($_) ){
		push(@reBacteria, $_);
	}
}


#Loops over all leaves and puts the speciesnames in hash
foreach (@leaves){
	
	my $firstName = getName($_, 1);
	
	#Checks that all species names coorrespond to the ones in the groupfile
	typoLocator($firstName) or die "ERROR: Typo (" . $firstName . ") in treefile: " . $treefile . "\n Please correct this before running Tree Doctor\n" ;
	
	#first and second part of id
	$speciesName = getName($_, 2);
	
	$leavesName{$speciesName}++;
}

#Check if there are >1 of one species in tree
for (keys %leavesName){
	if($leavesName{$_} > 1){
		push(@paralogSpecies,$_);
	}
}

my @nodesToRemove;

foreach (@paralogSpecies){
	my @paralogNodes;
	my $paralog = $_;
	
	for(@leaves){
		
		if(getName($_, 2) eq $paralog){
		push(@paralogNodes, $_);
		}
	}
	
	my %branchLengths;
	for(@paralogNodes){
		$branchLengths{$_->id} = $_->branch_length();
	}

	
	#Node with shortest branch is detected
	my $min_key;
	my $min_value = 100;
	while ((my $key, my $value) = each %branchLengths) {
  		if ($value < $min_value) {
    		$min_value = $value;
    		$min_key = $key;
  		}
	}

	for (keys %branchLengths){
		if ($_ ne $min_key){
		
			$tree->remove_Node($tree->find_node(-id => $_))
			
		}
	}

}


my @paraFreeLeaves = $tree->get_leaf_nodes();

print "Number of leaves after paralogs are removed: " . @paraFreeLeaves . "\n";

#			 	
# Subroutines
#

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

sub checkIfMonophyletic{
	my @nodes = @{$_[0]};
	my $lca = $tree->get_lca(-nodes => \@nodes);
	my @children;
	
	foreach ($lca->get_all_Descendents()){
		if	($_->is_Leaf()){
			push(@children, $_);
		}
	}	
			
	if (checkIfWithinTheSameGroup(\@children)){
		return 1;
	}		
	else {
		return 0;
	}	

}	

sub checkIfWithinTheSameGroup{
	my @nodeArray = @{$_[0]};
	my @compareArray;	
		
	for(@nodeArray){
		
		push(@compareArray, getName($_, 1));
		
	}
	
	my $groupName = $groups{$compareArray[0]};
	
	foreach (@compareArray){
	
		if ($groups{$_} ne $groupName){
			return 0;
		}
	}
	
	return 1;	
		
}	 

# First parameter: Node
#Second parameter: 1 = return only first part of name. 2 = return first and second part
sub getName{
	my $node = $_[0];
	my $namePos = $_[1];
	my $outName;
	
	my @splitSpecies = split("_", $node->id); 
	
	if ($namePos == 1){
		return $splitSpecies[0];		
	}
	else {	
		return $splitSpecies[0] . "_" . $splitSpecies[1];
	}
	
}