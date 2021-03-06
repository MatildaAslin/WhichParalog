#!/usr/bin/perl
#Main author: Matilda Aslin

package Paralog;

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

our $VERSION = '1.00';
use base 'Exporter';

our @EXPORT = qw(paralogRemover);

sub paralogRemover{


	#Declaring variables
	my $treefile = $_[0];
	my $groupfile = $_[1];
	my $leaf;
	my %groups;
	my $speciesName;
	my %leavesName;
	my @paralogSpecies;
	my @bacteria;
	my @reBacteria;

	#Load tree file
	my $str = Bio::TreeIO->new(-file =>$treefile, -format => 'nexus');
	my $tree = $str->next_tree();

	#Opens filehandle to LOG file
	open (MYFILE, '>>Paralog-LOG.txt');
	
	#Opens filehandle to groupfile 
	open(INFILE, $groupfile) or die "Can't open file: $!\n";

	#Make hash out of groupfile (Key: Species, Value: Group) 
	foreach my $line (<INFILE>) { 

		chomp($line);
		
		my ($groupname, $protein) = split("\t", $line);
	
		my @species = split("_", $protein);
	
		#Special case for Candidatus sepcies
		if ($species[0] eq "Ca"){
		
			$groups{$species[0] . "_" . $species[1]} = $groupname;
	
		} 
	
		else {
	
		$groups{$species[0]} = $groupname;
	
		}
	}



	#Puts all leaves in an array
	my @leaves = $tree->get_leaf_nodes();


	#Put all bacterial nodes in an array
	for (@leaves){
	
		
		if( checkIfBacteria($_, \%groups) ){
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
	if (checkIfMonophyletic(\@bacteria, $tree, \%groups)){
		
		#re-root tree with bacteria (the first one in the array)
		$tree->reroot($tree->get_lca(-nodes => \@bacteria)); 
	}
	else {
		die "Bacteria in tree are not monophyletic!";
	}

	#Get leaves of re-rooted tree
	my @reLeaves = $tree->get_leaf_nodes();

	print MYFILE "Number of leaves after re-rooting: " . @reLeaves . "\n";

	#gets bacterial nodes of re-rooted tree
	for (@reLeaves){
		if( checkIfBacteria($_, \%groups) ){
			push(@reBacteria, $_);
		}
	}


	#Loops over all leaves and puts the speciesnames in hash
	foreach (@reLeaves){
	
		my $keyName = getName($_, 1);
	
		if(getName($_, 1) eq "Ca"){
			$keyName = getName($_, 2)
		}
	
		#Checks that all species names correspond to the ones in the groupfile
		typoLocator($keyName, \%groups) or die "ERROR: Typo (" . $keyName . ") in treefile: " . $treefile . "\n Please correct this before running Tree Doctor\n" ;
	
		#first and second part of id
		$speciesName = getName($_, 2);
	
		#Species names are put in a hash
		$leavesName{$speciesName}++;
	}

	#Checks if there are >1 of one species in tree
	#They are identified as a species containing paralogs
	for (keys %leavesName){
		if($leavesName{$_} > 1){
			push(@paralogSpecies,$_);
		}
	}
  
	foreach (@paralogSpecies){
		my @paralogNodes;
		my $paralog = $_;
		
		#Goes through all leaves and get all OTUs for $paralog	
		for(@reLeaves){
			
			if(getName($_, 2) eq $paralog){
				push(@paralogNodes, $_)
			}
		}
		
		print MYFILE "\n" . "ParalogSpecies: " .  $paralog . "\n";
		
		my @paralogClade = getClade(\@paralogNodes, $tree);
		
		print MYFILE "Group: " . $groups{getName($paralogNodes[0], 1)} . "\n";
	
		#Monophyletic index for paralogclade	
		my $startIndex = checkMonophyleticIndex(\@paralogClade, $groups{getName($paralogNodes[0], 1)}, \%groups);
	
		my @newParalogClade;
	
		print MYFILE "StartIndex: " . $startIndex . "\n"; 
	
		#If monophyletic index is 1, it can not be improved, 
		#so monophyletic optimizing is only performed is if $startIndex < 1
		unless ($startIndex == 1){
		
			print MYFILE "Initializing monophyletic optimizing!\n";
		
			print MYFILE "ParalogNodes:\n"; 
			for (@paralogNodes){
				print MYFILE $_->id . "\n"; 
			}
		
			my $maxIndex = $startIndex;
			my @copyNodes = @paralogNodes;
			my $candidate;
		
			for my $paralogNode (@paralogNodes){
		
				my @removedParalog;
				my $monoIndex;
					
				#Nodes are gathered to be able to get new clade without $paralogNode
				#They are put into @removedParalog
				for my $node (@copyNodes){
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
			
				print MYFILE "ParalogNode: " . $paralogNode->id . "\n";
				
				print MYFILE "parlogGroup: " . $groups{getName($paralogNode, 1)} . "\n";
			
				@newParalogClade = getClade(\@removedParalog, $tree);
			 
			 	#Checks monophyletic index for new clade
				$monoIndex = checkMonophyleticIndex(\@newParalogClade, $groups{getName($paralogNodes[0], 1)}, \%groups);
		
				print MYFILE "MonoIndex: " . $monoIndex . " MaxIndex: " . $maxIndex . "\n";
				
				#If monoIndex is greater than maxIndex then
				#the current $paralogNode is set as candidate for removal
				if ($monoIndex > $maxIndex){
					$maxIndex = $monoIndex;
					$candidate = $paralogNode;
				
				}
			
				if ($candidate){ 
					print MYFILE "Candidate: " . $candidate->id . "\n";
				}
			}
		
			print MYFILE "maxIndex: " . $maxIndex . " > " . "StartIndex: " . $startIndex . "\n" ;
			
			#If the final $MaxIndex is greater than the $startIndex --> $candidate is removed
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
	
	
	
		print MYFILE "Initializing branchTrimming!\n";
	
		my %branchLengths;
		
		#Hash is made with node as key and branch length as value
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


		#All paralogs are removed except the one the the shortest branch
		for my $key (keys %branchLengths){
			
				for my $leaf (@reLeaves){
				
					if ($leaf eq $key && $leaf ne $min_key ){
					
 						print MYFILE "Removing the following OTU because of long branch length: " . $leaf->id . "\n"; 
					
						$tree->remove_Node($leaf);
					}
				}		 
			
				
				my @paraFreeLeaves = $tree->get_leaf_nodes();
			
				#Leaf nodes with no IDs are removed
				for(@paraFreeLeaves){

					if(!$_->id){
				
						$tree->remove_Node($_);
				
					}

				}
			
			
		
		}

	}


	my @paraFreeLeaves = $tree->get_leaf_nodes();

	print MYFILE "Number of leaves after paralogs are removed: " . @paraFreeLeaves . "\n";
	
	close (MYFILE);
	
	#Final tree is returned
	return $tree;
	
}

#			  # 	
# Subroutines #
#			  #

sub checkIfBacteria{

	my $nodeName = getName($_[0], 1);
	my %groups = %{$_[1]};
	
	if ($nodeName eq "Ca"){
		$nodeName = getName($_[0], 2);
	}
	
	if ($groups{$nodeName} eq "bact"){
		return 1;
	}
}

sub typoLocator{

	my %groups = %{$_[1]};
	
	if ($groups{$_[0]}){
		return 1;
	}
}


sub checkIfMonophyletic{
	my @nodes = @{$_[0]};
	my $tree = $_[1];
	my %groups = %{$_[2]};
	
	my @children = getClade(\@nodes, $tree);
			
	if (checkIfWithinTheSameGroup(\@children, \%groups)){
		return 1;
	}		
	else {
		return 0;
	}	

}	

#Input: Array with nodes and group hash. Returns true if all nodes belong to the same group, otherwise false.
sub checkIfWithinTheSameGroup{
	my @nodeArray = @{$_[0]};
	my %groups = %{$_[1]};
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

#Input: Array with nodes and tree. Returns all leaf nodes within clade, 
#where clade is defined as all the decendents of the last common ancestor
#of all nodes in @nodes.
sub getClade{
	
	my @nodes = @{$_[0]};
	my $tree = $_[1];
	my $lca = $tree->get_lca(-nodes => \@nodes);
	my @children;
	
	foreach ($lca->get_all_Descendents()){
		if	($_->is_Leaf() && $_->id){
			push(@children, $_);
		}
	}
	
	return @children;
}

#Input: Array with nodes, groupname and group hash.
#Counts the fraction of groupname among nodes
sub checkMonophyleticIndex{
	
	my @children = @{$_[0]};
	my $groupName = $_[1];
	my %groups = %{$_[2]};
		
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