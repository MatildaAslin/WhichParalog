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
my @blastinfo;
my @pos=(0,0);

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
		$seq = $aln->get_seq_by_pos($i+1);
		my @protein=split(//, $seq->seq);
		my $start=0;
		my $stop=0;
		my $longest=$stop-$start;
		my $seclong=0;
		my @secpos=(0,0);
		for (my $j=0;$j<scalar(@protein);$j++){
 			if($protein[$j] eq "-"){
	 			$start=$j;
	 			while(($j<scalar(@protein))&&($protein[$j] eq "-")) {
		 			$j++;
		 		}
		 		$stop=$j;
		 		if($longest<$stop-$start){
			 		$longest=$stop-$start;
			 		@pos=($start,$stop);
			 	}
			 	if(($seclong<$stop-$start)&&($longest>$stop-$start)){
				 	$seclong=$stop-$start;
			 		@secpos=($start,$stop);
			 	}
	 		}
		}
		if (($secpos[0]-$pos[1]<10) && ($secpos[0]!=0)){
			@pos=($pos[0],$secpos[1]);
		}
		push(@blastinfo, $name[$i], $pos[0], $pos[1]);
	}
}

my $blastaln=$aln->slice($blastinfo[$k*3+1],$blastinfo[$k*3+2]);

print $blastinfo[0]."\n";
print $blastinfo[1]."\n";
print $blastinfo[2]."\n";
print $blastinfo[3]."\n";
print $blastinfo[4]."\n";
print $blastinfo[5]."\n";
print $blastinfo[6]."\n";
print $blastinfo[7]."\n";
print $blastinfo[8]."\n";
print $blastinfo[9]."\n";
print $blastinfo[10]."\n";
print $blastinfo[11]."\n";

# 		my $factory = Bio::Tools::Run::AnalysisFactory::Pise->new();
#         my $psiblast = $factory->program('psiblast');
# 		my $gene=$aln->get_seq_by_pos($i+1);
# 		print $gene->seq;
# 	print $aln->length()."\n";
# 	print "@len"."\n";