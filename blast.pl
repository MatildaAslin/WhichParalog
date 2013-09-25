#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;
use File::Slurp;
use File::Temp;
use Bio::TreeIO;
use IO::String;
use Bio::AlignIO;
use Bio::Align::ProteinStatistics;
use Bio::Root::Root;
use Bio::Tools::GuessSeqFormat;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Tree::DistanceFactory;
use Bio::Tree::TreeI;
use Bio::SimpleAlign;
use Bio::SearchIO;
use Bio::SearchIO::Writer::HitTableWriter;

#Declaring variables
my $alignment;
my $db_dir;
my $seq;
my @len;
my @name;
my @id;
my @blastinfo;
my @pos=(0,0);

#my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(-program  => 'psiblast', -DB_NAME => 'Geoarchaeon_NAG1.faa', -DB_DIR => '/Users/Matilda/WhichParalog/testdb');

#Saving input parameters
GetOptions ("a|alignment=s" => \$alignment, "db|db_dir=s" => \$db_dir );

#Load alignment file
my $str = Bio::AlignIO->new(-file =>$alignment, -format => 'fasta');
my $aln = $str->next_aln();
my $str2 = Bio::AlignIO->new(-file =>$alignment, -format => 'fasta');
my $aln2 = $str2->next_aln();

$aln->sort_alphabetically;
$aln2->sort_alphabetically;

foreach $seq($aln2->each_seq){ #Take away the gaps to count the length
	my $tmp = $seq->seq();
	$tmp =~ s/-//g;
	$seq->seq($tmp);
	push(@len, length($seq->seq));
	@id = split(/_/, $seq->id);
	push(@name, $id[0]."_".$id[1]); #Get the name of the organism
}

for(my $i=0; $i<$aln->no_sequences; $i++){
	if($len[$i]/$aln->length()<0.7){ #See if the sequence cover less than 70% of alignment
		print "Need new blast for sequence: ".$name[$i]."\n"; #Just to check if it is working
		$seq = $aln->get_seq_by_pos($i+1); #get sequence with gap
		my @protein=split(//, $seq->seq); #put sequence in array
		my $start=0;
		my $stop=0;
		my $longest=$stop-$start;
		my $seclong=0;
		my @secpos=(0,0);
		for (my $j=0;$j<scalar(@protein);$j++){	#Find the longest gap (and the seconde longest)
 			if($protein[$j] eq "-"){
	 			$start=$j;
	 			while(($j<scalar(@protein))&&($protein[$j] eq "-")) {
		 			$j++;
		 		}
		 		$stop=$j;
		 		if($longest<$stop-$start){	#Gives the positions for the longest gap
			 		$longest=$stop-$start;
			 		@pos=($start,$stop);
			 	}
			 	if(($seclong<$stop-$start)&&($longest>$stop-$start)){ #Gives the positions for the secondlongest gap
				 	$seclong=$stop-$start;
			 		@secpos=($start,$stop);
			 	}
	 		}
		}
		if (($secpos[0]-$pos[1]<10) && ($secpos[0]!=0)){ #If the longest and secondlongest gap is close together we merge them
			@pos=($pos[0],$secpos[1]);
		}
		my $blastaln=$aln->slice($pos[0]+1,$pos[1]); #Take the alignment for the gap and use for blast
#		my $species = $name[$i]; # The first two things in the name
		my $species = "Geo";
		my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(-program  => 'psiblast', -DB_NAME => $db_dir . "/" . $species);
		my $psiblast = $factory->psiblast(-query => $blastaln);
		#print $psiblast;
		my $writer = Bio::SearchIO::Writer::HitTableWriter->new();
		my $blio = Bio::SearchIO->new( -file => ">damp", -format=>'blast', -writer => $writer );
		$blio->write_result($psiblast);
	}
}
		

#Things to think about
#Should we compare strings with the same name first and if one of them are longer than 70% we don't have to look on the other?




# 		my $gene=$aln->get_seq_by_pos($i+1);
# 		print $gene->seq;
# 	print $aln->length()."\n";
# 	print "@len"."\n";