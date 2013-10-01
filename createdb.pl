#!/usr/bin/perl
use strict;
use warnings;
use Bio::Tools::Run::StandAloneBlastPlus;

$fac = Bio::Tools::Run::StandAloneBlastPlus->new(
	-db_name => 'Ca_Nanosalina',
	-db_data => 'Ca_Nanosalina_sp_J07AB43.faa',
	-create => 1);