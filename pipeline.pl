#!/usr/bin/perl
use strict;
use warnings;
use File::Slurp;
use File::Temp;
use Getopt::Long;
use Bio::TreeIO;
use IO::String;
use Bio::AlignIO;
use Bio::Align::ProteinStatistics;
use Bio::Tree::DistanceFactory;
use Bio::Tree::TreeI;
use Bio::SimpleAlign;
use Bio::Tools::Run::Alignment::MAFFT;
#use Paralog;
use Split;
#use Blast;

my $alignment;
my $db_dir;

GetOptions ("a|alignment=s" => \$alignment, "db|db_dir=s" => \$db_dir );

my $str = Bio::AlignIO->new(-file =>$alignment, -format => 'fasta');
my $aln = $str->next_aln();

#run splitgene

my $splitAlign = split_gene($aln);

my $out = new Bio::AlignIO(-file => '>outputAln', -format => 'fasta');

$out->write_aln($splitAlign);


for my $seq ($splitAlign->each_seq){
	print $seq->seq . "\n";
}

#run blast

#run splitgene

#make tree

#run paralog 