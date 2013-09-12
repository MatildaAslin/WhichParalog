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
my @len;
my @name;
my @id;


#Saving input parameters
GetOptions ("i|input=s" => \$input);

#Load alignment file
my $str = Bio::AlignIO->new(-file =>$input, -format => 'fasta');
my $aln = $str->next_aln();
my $str2 = Bio::AlignIO->new(-file =>$input, -format => 'fasta');
my $aln2 = $str2->next_aln();

$aln->sort_alphabetically;
$aln2->sort_alphabetically;

foreach $seq($aln2->each_seq){
	my $tmp = $seq->seq();
	$tmp =~ s/-//g;
	$seq->seq($tmp);
	push(@len, length($seq->seq));
	@id = split(/_/, $seq->id);
	push(@name, $id[0]." ".$id[1]);
}

for(my $i=0; $i<$aln->no_sequences; $i++){
	if($len[$i]/$aln->length()<0.7){
		print "Need new blast for sequences number: ".($i+1)."\n";
		push(my @bad, $name[$i]);
		my @protein=split(//, $aln->get_seq_by_pos($i+1));
		my $pos=0;
		my $start=0;
		for (my $j=0;$j<scalar(@protein);$j++){
			if($protein[$j] ne "-"){
				$pos++;
			}
			if($protein[$j] eq "-"){
				$start=$pos;
			}
		}
	}
}


# 		my $gene=$aln->get_seq_by_pos($i+1);
# 		print $gene->seq;
# 	print $aln->length()."\n";
# 	print "@len"."\n";