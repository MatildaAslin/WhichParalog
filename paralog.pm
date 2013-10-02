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


sub paralogRemover{

my $treefile = $_[0];
my $groupfile = $_[1];

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

open (MYFILE, '>>Paralog-LOG.txt');

#Saving input parameters
#"GetOptions ("t|treefile=s" => \$treefile, "g|groupfile=s" => \$groupfile);

open(INFILE, $groupfile) or die "Can't open file: $!\n";

#Make hash out of groupfile (Key: Species, Value: Group) 
foreach my $line (<INFILE>) { 

	chomp($line);
		
	my ($groupname, $protein) = split("\t", $line);
	
	my @species = split("_", $protein);
	
	if ($species[0] eq "Ca"){
		
		$groups{$species[0] . "_" . $species[1]} = $groupname;
	
	} 
	
	else {
	
	$groups{$species[0]} = $groupname;
	
	}
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

#Checks if there are less than 3 bacteria in the tree
if ($noBacteria < 3){
	die "ERROR: Tree has less than 3 bacteria.\n";
}

#Checks if the bacterial species in the tree are monophyletic 
if (checkIfMonophyletic(\@bacteria)){
	#re-root tree with bacteria (the first one in the array)
	$tree->reroot($tree->get_lca(-nodes => \@bacteria)); 
}
#Temporary else-statement
else {
	die "Bacteria in tree are not monophyletic!";
}

#my $treeio1 = new Bio::TreeIO(-file => '>outputTree-ReRoot.newick', -format => 'newick');

#$treeio1->write_tree($tree);

#Get leaves of re-rooted tree
my @reLeaves = $tree->get_leaf_nodes();

print MYFILE "Number of leaves after re-rooting: " . @reLeaves . "\n";

#gets bacterial nodes of re-rooted tree
for (@reLeaves){
	if( checkIfBacteria($_) ){
		push(@reBacteria, $_);
	}
}


#Loops over all leaves and puts the speciesnames in hash
foreach (@reLeaves){
	
	my $keyName = getName($_, 1);
	
	if(getName($_, 1) eq "Ca"){
		$keyName = getName($_, 2)
	}
	
	#Checks that all species names coorrespond to the ones in the groupfile
	typoLocator($keyName) or die "ERROR: Typo (" . $keyName . ") in treefile: " . $treefile . "\n Please correct this before running Tree Doctor\n" ;
	
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
  
foreach (@paralogSpecies){
	my @paralogNodes;
	my $paralog = $_;
	
	for(@reLeaves){
		
		if(getName($_, 2) eq $paralog){
			push(@paralogNodes, $_)
		}
	}
		
	print MYFILE "\n" . "ParalogSpecies: " .  $paralog . "\n";
	
	my @paralogClade = getClade(\@paralogNodes);
	
#	print MYFILE "ParalogClade:\n"; 
#	for (@paralogClade){
#		print MYFILE $_->id . "\n"; 
#	}
	
	print MYFILE "Group: " . $groups{getName($paralogNodes[0], 1)} . "\n";
	my $startIndex = checkMonophyleticIndex(\@paralogClade, $groups{getName($paralogNodes[0], 1)});
	my @newParalogClade;
	
	print MYFILE "StartIndex: " . $startIndex . "\n"; 
	
	unless ($startIndex == 1){
		
		print MYFILE "Initializing monophyletic optimizing!\n";
		
		print MYFILE "ParalogNodes:\n"; 
		for (@paralogNodes){
			print MYFILE $_->id . "\n"; 
		}
		
		my $maxIndex = $startIndex;
		my @damp = @paralogNodes;
		my $candidate;
		
		for my $paralogNode (@paralogNodes){
		
			my @removedParalog;
			my $monoIndex;
					
			for my $node (@damp){
				unless ($paralogNode->internal_id eq $node->internal_id){
				
					my @kids;
					my $ancestor;
					my @decendents;
					
					if (@paralogNodes == 2){
						while (@kids < 2){
							$ancestor = $node->ancestor;
							@decendents = $ancestor->get_all_Descendents();
					
							for (@decendents){
								if ($_->is_Leaf && $_->id && $_->id ne $paralogNode->id){
									push(@kids, $_) 
								}
							}
						
						$node = $ancestor;
	
						}	
					}
					else{
					
						$ancestor = $node->ancestor;
						@decendents = $ancestor->get_all_Descendents();
					
						for (@decendents){
							if ($_->is_Leaf && $_->id && $_->id ne $paralogNode->id){
								push(@kids, $_) 
							}
						}
					
					}		
							
					push(@removedParalog, @kids);
					
				}
			
			
			}
		
#			print MYFILE "before print MYFILEing removedPAralog\n";
#			for (@removedParalog){
#				print MYFILE $_->id . "\n";
#			}
		
#			print MYFILE "after print MYFILEing removedParalog\n";
			
			print MYFILE "ParalogNode: " . $paralogNode->id . "\n";
				
			print MYFILE "parlogGroup: " . $groups{getName($paralogNode, 1)} . "\n";
			
			@newParalogClade = getClade(\@removedParalog);
			 
			$monoIndex = checkMonophyleticIndex(\@newParalogClade, $groups{getName($paralogNodes[0], 1)});
		
			print MYFILE "MonoIndex: " . $monoIndex . " MaxIndex: " . $maxIndex . "\n";
			if ($monoIndex > $maxIndex){
				$maxIndex = $monoIndex;
				$candidate = $paralogNode;
				
			}
			
			if ($candidate){ 
				print MYFILE "Candidate: " . $candidate->id . "\n";
			}
		}
		
		print MYFILE "maxIndex: " . $maxIndex . " > " . "StartIndex: " . $startIndex . "\n" ;
		if ($maxIndex > $startIndex){
			print MYFILE "The following OTU is removed to improve monophyly: " . $candidate->id . "\n";
			$tree->remove_Node($candidate);
			
			#Update ParalogNodes
			my @temp;
			for (@paralogNodes){
				unless ($_->id eq $candidate->id){
					push(@temp, $_);
				}
			}
			
			@paralogNodes = @temp;
		} 
	
	}
	
	
	
	print MYFILE "Initializing treeTrimming!\n";
	
	my %branchLengths;
	for(@paralogNodes){
		
		$branchLengths{$_} = $_->branch_length();
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



	for my $key (keys %branchLengths){
		    
		    for my $leaf (@reLeaves){
		    	if ($leaf eq $key && $leaf ne $min_key ){
					$tree->remove_Node($leaf);
				}
			}        
			
			my @paraFreeLeaves = $tree->get_leaf_nodes();
			
			for(@paraFreeLeaves){

				if(!$_->id){
				
					$tree->remove_Node($_);
				
				}

			}
			
			
		
	}

}


my @paraFreeLeaves = $tree->get_leaf_nodes();

print MYFILE "Number of leaves after paralogs are removed: " . @paraFreeLeaves . "\n";

my $treeio = new Bio::TreeIO(-file => '>outputTree.newick', -format => 'newick');

$treeio->write_tree($tree);

}

#			  # 	
# Subroutines #
#			  #

sub checkIfBacteria{

	my $nodeName = getName($_[0], 1);
	
	if ($nodeName eq "Ca"){
		$nodeName = getName($_[0], 2);
	}
	
	if ($groups{$nodeName} eq "bact"){
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
#	my $lca = $tree->get_lca(-nodes => \@nodes);
#	my @children;
	
#	foreach ($lca->get_all_Descendents()){
#		if	($_->is_Leaf()){
#			push(@children, $_);
#		}
#	}	
	my @children = getClade(\@nodes);
			
	if (checkIfWithinTheSameGroup(\@children)){
		return 1;
	}		
	else {
		return 0;
	}	

}	

#Input: Array with nodes. Returns true if all nodes belong to the same group, otherwise false.
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

#Input: Array with nodes. Returns all leaf nodes within clade, where clade is defined as all the decendents of the last common ancestor
#of all nodes in @nodes.
sub getClade{
	
	my @nodes = @{$_[0]};
	my $lca = $tree->get_lca(-nodes => \@nodes);
	my @children;
	
	foreach ($lca->get_all_Descendents()){
		if	($_->is_Leaf() && $_->id){
			push(@children, $_);
		}
	}
	
	return @children;
}

#Input: Array with nodes and groupname
sub checkMonophyleticIndex{
	
	my @children = @{$_[0]};
	my $groupName = $_[1];
		
	my $total = @children;
	my $counter = 0;
	
	
	
	foreach my $child (@children){
	
		if ($groups{getName($child, 1)} eq $groupName){
			$counter++;
		}
	}
	
	return $counter/$total;

}

# First parameter: Node
#Second parameter: 1 = return only first part of name. 2 = return first and second part
sub getName{
	my $node = $_[0];
	my $namePos = $_[1];
	my $outName;
	
	my @splitSpecies = split("_", $node->id); 
	
	if ((($namePos != 1 ) || ($splitSpecies[0] eq "Ca")) && (!$splitSpecies[1] eq "")){
		return $splitSpecies[0] . "_" . $splitSpecies[1];		
	}
	else {	
		return $splitSpecies[0];
	}
	
}