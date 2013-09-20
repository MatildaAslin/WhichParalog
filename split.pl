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

for (my $i=1; $i<scalar(@all) - 1; $i++){
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
	    my $currentAln = $alnArray[$i]->seq; # This removing the '-' charachter
        my $nextAln = $alnArray[$i+1]->seq;  # and then adding it again is not very pretty :(
        
        #make bitwise or to find overlaps
        my $or = $currentAln | $nextAln;
        $or =~ s/ /-/g;
        (my $overlap = $or) =~ s/[^A-Z]//g;
        
        if (length $overlap <= 10) { # What overlap threshold is good?
        # calculate mean of alignment
            if (calcMean($currentAln) < calcMean($nextAln)) {
                $alnArray[$i]->seq($currentAln . $nextAln);
                $aln->remove_seq($alnArray[$i+1]);
            } else {
                $alnArray[$i]->seq($nextAln . $currentAln);
                $aln->remove_seq($alnArray[$i+1]);
            }
        }
	}
}

#removing gaps manually, remove_gaps() gives weird results.
foreach $seq($aln->each_seq()){        
    (my $noGaps = $seq->seq) =~ s/-//g;
    $seq->seq($noGaps);
}

# At the moment the new aligment is written to an output file
# For future: realing the file and maybe give out the alignment object
my $alnio = new Bio::AlignIO(-file => '>outputAln.fasta', -format => 'fasta');
$alnio->write_aln($aln);

# calculate the mean of a sequence
# sum up the positions of the sequence that are not a gap
# then divied it by the number of positions
sub calcMean {
    my $sequence = $_[0];
    my @splitSeq = split(//, $sequence);
    
    my $mean = 0;
    my $counter = 0;
    
    for (my $i=1; $i<scalar(@splitSeq)-1; $i++) {
        if ($splitSeq[$i] ne '-') {
            $mean = $mean + $i;
            $counter++;
        }
    }
    return ($mean/$counter);
}
