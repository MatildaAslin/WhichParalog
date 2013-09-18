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
	if ($all[$i] eq $all[$i+1]) {
	    
		#go through the genes and check if they can be split-genes, 
		#save candidat incase more than one copie
		#If it is a candidate, save the number i so it can be 
		#removed it it actually is a split gene
		#print $all[$i]."\n";
	}
}

my @alnArray = $aln->each_seq;
for (my $i=1; $i<scalar(@alnArray) - 1; $i++) {

    #compare if species is the same, i.e. the names are equal
    my @nameCurrent = split(/_/, $alnArray[$i]->id);
    my @nameNext = split(/_/, $alnArray[$i+1]->id);
    #species are defined by the first two words of their id
	if ($nameCurrent[0].$nameCurrent[1] eq $nameNext[0].$nameNext[1]) {
	    (my $currentAln = $alnArray[$i]->seq) =~ s/-/ /g;
        (my $nextAln = $alnArray[$i+1]->seq) =~ s/-/ /g;
        
        #make bitwise or to find overlaps
        my $or = $currentAln | $nextAln;
        (my $overlap = $or) =~ s/[^A-Z]//g;
        
        if (length $overlap == 0) {
            print "$nameCurrent[0].$nameCurrent[1], $nameNext[0].$nameNext[1]:\n";
            print "OR: $or\n";
        }
	}
}

# make sure input is a sequence object. Count gaps in sequence
sub countGaps {
    my $sequence = $_[0];
    #split String to get each character individually
    my @splitSeq = split(//, $sequence);

    #counters for current number of gaps and maximum number of gaps
    my $gapCounter = 0;
    my $maxGaps = 0;

    foreach my $char (@splitSeq) {
        if ($char eq '-') { #count number of gaps in a row
            ++$gapCounter;
        } elsif ($gapCounter > $maxGaps) {
            $maxGaps = $gapCounter;
            $gapCounter = 0;
        } else {
            $gapCounter = 0;
        }
    }
    if ($gapCounter > $maxGaps) {
            $maxGaps = $gapCounter;
    }
    
    return $maxGaps;
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
