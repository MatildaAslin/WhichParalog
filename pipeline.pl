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
use Paralog;
use Split;
use Blast;

my $alignment;
my $db_dir;
my $treefile;
my $groupfile;

GetOptions ("a|alignment=s" => \$alignment, "db|db_dir=s" => \$db_dir, "t|treefile=s" => \$treefile, "g|groupfile=s" => \$groupfile);

my $str = Bio::AlignIO->new(-file =>$alignment, -format => 'fasta');
my $aln = $str->next_aln();

#run splitgene

my $splitAlign = split_gene($aln);

#run blast

if(my $blastAln = splitBlast($splitAlign, $db_dir)){
	my $splitAlign2 = split_gene($blastAln);
	my $out = new Bio::AlignIO(-file => '>splitBlastAln.phy', -format => 'fasta');
	$out->write_aln($splitAlign2);
}
else {
	my $out = new Bio::AlignIO(-file => '>splitAln.phy', -format => 'fasta');
	$out->write_aln($splitAlign);
}

#make a tree of the alignment file "splitAlign2

#run paralog 

my $tree = paralogRemover($treefile, $groupfile);

my $treeio = new Bio::TreeIO(-file => '>outputTree.newick', -format => 'newick');

$treeio->write_tree($tree);


