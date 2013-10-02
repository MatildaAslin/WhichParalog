#!/usr/bin/perl

package Split;

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
use Bio::Tools::Run::Alignment::Muscle;

our $VERSION = '1.00';
use base 'Exporter';

our @EXPORT = qw(split_gene);

sub split_gene {
    #Declaring variables
    my $seq;
    my $factory = Bio::Tools::Run::Alignment::Muscle->new();

    my $aln = $_[0];

    $aln->sort_alphabetically;

    my $found_split = 0; # boolean flag for checking if file needs to be changed
    my @alnArray = $aln->each_seq;
    for (my $i=1; $i<scalar(@alnArray) - 1; $i++) {

        #compare if species is the same, i.e. the names are equal
        my @nameCurrent = split(/_/, $alnArray[$i]->id);
        my @nameNext = split(/_/, $alnArray[$i+1]->id);
        #species are defined by the first two words of their id
	    if ($nameCurrent[0].$nameCurrent[1] eq $nameNext[0].$nameNext[1]) {
	        my $currentAln = $alnArray[$i]->seq;
            my $nextAln = $alnArray[$i+1]->seq;
            
            #make bitwise or to find overlaps
            my $or = $currentAln | $nextAln;
            $or =~ s/ /-/g;
            (my $overlap = $or) =~ s/[^A-Z]//g;
            
            if (length $overlap <= 5) {
            $found_split = 1;
            # calculate mean of alignment
                if (calcMean($currentAln) < calcMean($nextAln)) {
                    $alnArray[$i]->seq($currentAln . $nextAln);
                    $aln->remove_seq($alnArray[$i+1]);
#                    print $nameCurrent[0].$nameCurrent[1] . ":\n" . $currentAln . "\n";
 #                   print $nameNext[0].$nameNext[1] . ":\n" . $nextAln . "\n";
                } else {
                    $alnArray[$i]->seq($nextAln . $currentAln);
                    $aln->remove_seq($alnArray[$i+1]);
  #                  print $nameCurrent[0].$nameCurrent[1] . ":\n" . $currentAln . "\n";
   #                 print $nameNext[0].$nameNext[1] . ":\n" . $nextAln . "\n";
                }
            }
	    }
    }
    
    my $align = $aln;
    
    if ($found_split) { # if there are no split genes in the file, nothing has to be done
        #removing gaps manually, remove_gaps() gives weird results.
        foreach $seq($aln->each_seq()){        
            (my $noGaps = $seq->seq) =~ s/-//g;
            $seq->seq($noGaps);
        }

        # At the moment the new aligment is written to an output file for realignment
#        my $alnio = new Bio::AlignIO(-file => '>outputAln.fasta', -format => 'fasta');
 #       $alnio->write_aln($aln);
		
		my @seq_array;
		
        for my $seq ($aln->each_seq){
        	push(@seq_array, $seq);
        }
        
        # Pass the factory a list of sequences to be aligned.
  #      my $inputfilename = 'outputAln.fasta';
        # $aln is a SimpleAlign object.
        $align = $factory->align(\@seq_array);
    }
    return $align;
}

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
