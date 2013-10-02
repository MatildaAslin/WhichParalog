#!/usr/bin/perl
use strict;
use warnings;
use Bio::Align::ProteinStatistics;
use Bio::AlignIO;
use Bio::LocatableSeq;
use Bio::Root::Root;
use Bio::SearchIO;
use Bio::SearchIO::Writer::HitTableWriter;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::Tools::GuessSeqFormat;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Tree::DistanceFactory;
use Bio::TreeIO;
use Bio::Tree::TreeI;
use IO::String;
use File::Slurp;
use File::Temp;
use Getopt::Long;
use Bio::Tools::Run::Alignment::MAFFT;

our $VERSION = '1.00';
use base 'Exporter';

our @EXPORT = qw(splitBlast);

sub splitBlast{
	#Declaring variables
	my $alignment;
	my $db_dir;
	my $seq;
	my @len;
	my @name;
	my @id;
	my @blastinfo;
	my @pos=(0,0);

	#Saving input parameters
	#GetOptions ("a|alignment=s" => \$alignment, "db|db_dir=s" => \$db_dir );

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

	#for(my $i=0; $i<$aln->no_sequences; $i++){
	for(my $i=0; $i<@name; $i++){
		if($len[$i]/$aln->length()<0.7){ #See if the sequence cover less than 70% of alignment
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
			my $species = $name[$i]; #namnet på databasen
			my $factory = Bio::Tools::Run::StandAloneBlastPlus->new(-program  => 'psiblast', -DB_NAME => $db_dir . "/" . $species);
			my $psiblast = $factory->psiblast(-query => $blastaln);
			my $writer = Bio::SearchIO::Writer::HitTableWriter->new();
			my $blio = Bio::SearchIO->new( -file => ">" . $name[$i] . "_BlastResult", -format=>'blast', -writer => $writer );
			$blio->write_result($psiblast);
		
			my @hits = $psiblast->hits;
	
			my @hitsplit = split('\|', $hits[0]->name);
		 
			my $hitID = $hitsplit[1];
		
			print $name[$i] . " => " . $hitID . "\n"; 
	
			my @gene;
			my $newgene;
	# 		
		
	# 		#Maybe you should ask them to choose the path were they have these files instead?
			my $file = "testdb/" . $name[$i];
			open FILEHANDLE, $file or die $!;
			my @organism = <FILEHANDLE>;
	# 		

		
			for(my $q=0;$q<scalar(@organism);$q++) { 
				@gene=split(/ /,$organism[$i]);	#Goes trough the organism file
				if($organism[$q] =~ m/$hitID/){	#Find the right gene in the organism file
					while($organism[$q+1] !~ m/>/){ #Get the sequence
						$newgene=$newgene.$organism[$q+1];
						$q++;
					}
					$newgene =~ s/\n//g; #removes enter in the sequence
				}
			}
			my $newseq=Bio::LocatableSeq->new(-id=>">".$name[$i], -seq=>$newgene); #Create sequence
			print $newseq->seq . "\n";

			my @seq_array;
			if ($newseq->length()< $blastaln->length()){ #If sequence can fit in gap, realign and get the aligned sequence
				$blastaln->add_seq($newseq);
				for my $seq ($blastaln->each_seq){
					push(@seq_array, $seq);
				}
			
				#Download Mafft: http://mafft.cbrc.jp/alignment/software/
				my $alnfactory=Bio::Tools::Run::Alignment::MAFFT->new(); #Don't know which parameters that should be here
				my $blastrealn = $alnfactory->align(\@seq_array);

				$newseq=$blastrealn->get_seq_by_pos($blastrealn->no_sequences);
				$aln->add_seq($newseq); #Put the aligened sequence in the alignment and will be treated as a split gene

			
			}
		
		}
	}

	return $aln;
		
}