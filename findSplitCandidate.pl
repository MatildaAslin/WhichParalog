#!/usr/bin/perl

use strict;
use warnings;
use Bio::Seq;


my $seq_obj = Bio::Seq->new(-seq => "MGVDDRIEELEEKLKKTPVNKATEKERGRLKSQIAQLKEEKQKK--QKGTGETSGYAVEKTGDATVALVGFPSVGKSTLLNELTNTEVDTGAYEFTTLEVNPGVLKYKGANIQILDVPGLIGGAADGRGGGTQVLSVVRNADLILVVLDPEEL--RDDEIRDEVYKACLRTDSKPPNMKIEKKDRVG--------------------------------------------------------------------------------------------------YHLQLILILKRRQ-------------------------------------------------------------L",
                         -alphabet => 'protein' );
                         
#get sequence as String object
my $sequence = $seq_obj->seq;

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

print $maxGaps . "\n";
