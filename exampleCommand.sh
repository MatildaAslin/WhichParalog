#split.pl command
perl split.pl -i example/arCOG00358_pruned.fasta

#blast.pl command
perl blast.pl -a example/arCOG00358_pruned.fasta -db /Users/Matilda/WhichParalog/testdb

#paralog.pl command
perl paralog.pl -t example/arCOG00415_colored.nex -g example/organisms_groups.tab
