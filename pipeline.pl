#!/usr/bin/perl

=pod

=head1 SYNOPSIS

Tree Doctor pipeline. It uses 3 modules in the following order: Split.pm, Blast.pm and Paralog.pm.

B<Split.pm> - Identifies "split gene" candidates and merge them if the overlap is less than a certain threshold. 
If two sequences are merged, the sequences are re-aligned with MAFFT (with default parameters). 

B<Blast.pm> - There can be a split gene present in an alignment even though only one part is present in the alignment. 
This module aims to find the (potentially) missing part by using psiblast. 
If a sequence is added by this module, Split.pm will re-run and (potentially) merge the sequence with its "split gene partner".

B<Paralog.pm> - This modules removes paralogs in two steps.
First it checks if a paralog decreases the monophyly of the group the species belong to. It is then removed. 
If there are any paralogs left, they are removed based on branch length, only keeping the one with the shortest branch length.

  
=head1 USAGE

perl pipeline.pl -a <alignmentfile> -db <path to blastdatabases> -t <treefile> -g <groupfile>

=head1 INPUT

=head2 -a | --alignment

Alignment file in fasta format.

=head2 -db | --db_dir

Path to folder where blastdatabses for the organisms in the alignment/tree are present.

=head2 -t | --treefile

Treefile in nexus format.

=head2 -g | --groupfile

Tab-delimited file with group and OTU name. 

For example: E<10>
eury	Aciduliprofundum_boonei_T469 

=head1 OUTPUT

=head2 Alignment file
 New alignment in fasta. If it been modified 

=head1 AUTHORS

Jessika Nordin (Jessika.Nordin.0726@student.uu.se) E<10>
Veronika Scholz (Veronika.Scholz.0998@student.uu.se)E<10>
Matilda Aslin (Matilda.Aslin.3036@student.uu.se)E<10>

=head1 DATE

Thu  3 Oct 2013 11:19:44 CEST

=cut

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

if(my $blastAln = splitBlast($aln, $db_dir)){
	my $splitAlign2 = split_gene($blastAln);
	my $out = new Bio::AlignIO(-file => '>splitBlastAln.fasta', -format => 'fasta');
	$out->write_aln($blastAln);
}
else {
	my $out = new Bio::AlignIO(-file => '>splitAln.fasta', -format => 'fasta');
	$out->write_aln($splitAlign);
}

#make a tree of the alignment file "splitAlign2

#run paralog 

my $tree = paralogRemover($treefile, $groupfile);

my $treeio = new Bio::TreeIO(-file => '>outputTree.newick', -format => 'newick');

$treeio->write_tree($tree);


