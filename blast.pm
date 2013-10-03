#!/usr/bin/perl

package Blast;

use base 'Exporter';
use Bio::Align::ProteinStatistics;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::Root::Root;
use Bio::SearchIO;
use Bio::SearchIO::Writer::HitTableWriter;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::Tools::GuessSeqFormat;
use Bio::Tools::Run::Alignment::MAFFT;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Tree::DistanceFactory;
use Bio::TreeIO;
use Bio::Tree::TreeI;
use IO::String;
use File::Slurp;
use File::Temp;
use Getopt::Long;
use Storable 'dclone';
use strict;
use warnings;

our $VERSION = '1.00';
our @EXPORT = qw(splitBlast);

sub splitBlast{
	
	#Declaring variables
	my $aln = $_[0];
	my $aln2 = dclone($_[0]);
	my $db_dir = $_[1];
	my $seq;
	my @len;
	my @name;
	my @id;
	my @blastinfo;
	my @pos=(0,0);
	my $realn;

	#Sort alignment alphabetically
	$aln->sort_alphabetically;
	$aln2->sort_alphabetically;
	
	#Take away the gaps to count the length of each sequence
	foreach $seq($aln2->each_seq){ 
		my $tmp = $seq->seq();
		$tmp =~ s/-//g;
		$seq->seq($tmp);
		push(@len, length($seq->seq));
		@id = split(/_/, $seq->id);
		push(@name, $id[0]."_".$id[1]); #Get the name of the organism
	}

	#### Finding sequences with big gaps ####
	for(my $i=0; $i<@name; $i++){
		if($len[$i]/$aln->length()<0.7){ #See if the sequence coverage is less than 70%
			$seq = $aln->get_seq_by_pos($i+1);	#Get the aligned sequence
			my @protein=split(//, $seq->seq);	#Split sequence into array
			my $start=0;
			my $stop=0;
			my $longest=$stop-$start;
			my $seclong=0;
			my @secpos=(0,0);
			
			#Find the longest gap (and the second longest)
			for (my $j=0;$j<scalar(@protein);$j++){	
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
					if(($seclong<$stop-$start)&&($longest>$stop-$start)){ #Gives the positions for the second longest gap
						$seclong=$stop-$start;
						@secpos=($start,$stop);	
					}
				}
			}
			
			if (($secpos[0]-$pos[1]<10) && ($secpos[0]!=0)){ #If the longest and second longest gap is close together they are merged
				@pos=($pos[0],$secpos[1]);
			}
			
			#### BLAST ####
			my $blastaln=$aln->slice($pos[0]+1,$pos[1]); #Take the alignment covering the gap and use for blast
			my $species = $name[$i]; #Name of the database that will be used
			my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(-program  => 'psiblast', -DB_NAME => $db_dir . "/" . $species);
			my $psiblast = $factory->psiblast(-query => $blastaln); #Doing a psiblast on the sliced alignment
			
			## Uncomment if you wish to get a blastresult file ##
#			my $writer = Bio::SearchIO::Writer::HitTableWriter->new();
#			my $blio = Bio::SearchIO->new( -file => ">" . $name[$i] . "_BlastResult", -format=>'blast', -writer => $writer );
#			$blio->write_result($psiblast); #Writing a file with the blast result
		
			my @hits = $psiblast->hits;
			my @hitsplit = split('\|', $hits[0]->name);
			my $hitID = $hitsplit[1];	#Get a uniqe id to use for finding the gene
			my @gene;
			my $newgene;
			my @seq_array;
			
			my $file = "testdb/" . $name[$i]; # CHANGE THE DIRECTORY TO THE FILE (or always question about it)
			open FILEHANDLE, $file or die $!;
			my @organism = <FILEHANDLE>;
		
			for(my $q=0;$q<scalar(@organism);$q++) { #Goes trough the organism file
				@gene=split(/ /,$organism[$i]);	
				if($organism[$q] =~ m/$hitID/){	#Find the right gene in the organism file
					while($organism[$q+1] !~ m/>/){ #Get the sequence
						$newgene=$newgene.$organism[$q+1];
						$q++;
					}
					$newgene =~ s/\n//g; #removes enter in the sequence
				}
			}
			
			my $newseq=Bio::LocatableSeq->new(-id=>">".$name[$i], -seq=>$newgene); #Create sequence of the gene
			
			#If sequence can fit in gap, realign and get the aligned sequence
			if ($newseq->length()< $blastaln->length()){ 
				$aln->add_seq($newseq);
				for my $seq ($aln->each_seq){ #make a sequence array object that is used for aligning
					push(@seq_array, $seq);
				}
			
				#Download Mafft: http://mafft.cbrc.jp/alignment/software/
				my $alnfactory=Bio::Tools::Run::Alignment::MAFFT->new(); ## Can choose parameters ##
				$realn = $alnfactory->align(\@seq_array);	#Realign with new sequence
			}
		}
	}
	
	#If new sequence is added, return alignment, else return 0
	if (@name != $realn->no_sequences){
		return $realn;
	}
	else{
		return 0;
	}	
}