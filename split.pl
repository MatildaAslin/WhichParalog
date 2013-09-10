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

foreach $seq($aln->each_seq){
	@name = split(/_/, $seq->id);
	push(@all, $name[0].$name[1]);
}

for (my $i=1; $i<scalar(@all); $i++){
	if ($all[$i] eq $all[$i-1]){
		#go through the genes and check if they can be split-genes, 
		#save candidat incase more than one copie
		#If it is a candidate, save the number i so it can be 
		#removed it it actually is a split gene
		print $all[$i]."\n";
	}
}

#Usefull code:
#  	add seq					$myalgin->add_seq($newseq);
# 	remove seq				$aln->remove_seq($seq);
#	get seq by pos			$seq=$aln->get_seq_by_pos($i):
#	number of genes			no_sequences()
#	new alignment of seq 1-5			$aln2=$aln->select(1,5); 
#	write new file			$aln->write_fasta(\*OUTPUT)
#   					print $aln->length, "\n";
#   					print $aln->no_residues, "\n";
#   					print $aln->is_flush, "\n";
#   					print $aln->no_sequences, "\n";
#	$factory = Bio::Tools::Run::Alignment::Muscle->new(-format =>  
#'fasta',  -verbose=>'', -quiet=>'', -log='inv.log');
#	aln=$factory->align($inputfilename);

#To-do list
# Check if there are genes from the same organism.
# Check these genes if they seem to be split. Check the length of gene. 
# See how many gaps there are (some kind of count for where the gap is).
# If gap on different places, put them together. Realigned.