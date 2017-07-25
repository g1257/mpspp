#!/usr/bin/perl
=pod
Copyright (c) 2009-2015, UT-Battelle, LLC
All rights reserved

[MPS++, Version 0.]

*********************************************************
THE SOFTWARE IS SUPPLIED BY THE COPYRIGHT HOLDERS AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED.

Please see full open source license included in file LICENSE.
*********************************************************

=cut
use warnings;
use strict;

use lib "../../PsimagLite/scripts";
use Make;

my ($arg) = @ARGV;

if (defined($arg) and -r "$arg" and $arg ne "Config.make") {
	my $cmd = "cp Config.make Config.make.bak";
	system($cmd);
	$cmd = "cp $arg Config.make";
	system($cmd);
}

my %provenanceDriver = (name => 'Provenance', aux => 1);
my %progGlobalsDriver = (name => 'ProgramGlobals', aux => 1);
#my %restartDriver = (name => 'RestartStruct', aux => 1);
#my %finiteLoopDriver = (name => 'FiniteLoop', aux => 1);
#my %utilsDriver = (name => 'Utils', aux => 1);
#my %su2RelatedDriver = (name => 'Su2Related', aux => 1);
#my %toolboxDriver = (name => 'toolboxdmrg',
#                     dotos => 'toolboxdmrg.o ProgramGlobals.o Provenance.o Utils.o');
my $dotos = "Provenance.o";
#my %observeDriver = (name => 'observe', dotos => $dotos);

my @drivers = (\%provenanceDriver,\%progGlobalsDriver);
#\%restartDriver,\%finiteLoopDriver,\%utilsDriver,
#\%observeDriver,\%toolboxDriver);

$dotos = "mpspp.o Provenance.o"; # RestartStruct.o FiniteLoop.o Utils.o ";
$dotos .= " ProgramGlobals.o";

#my $templates = DmrgDriver::createTemplates();

#for (my $i = 0; $i < $templates; ++$i) {
#	my $name = "DmrgDriver$i";
#	my %dmrgDriver = (name => $name, aux => 1);
#	push @drivers,\%dmrgDriver;
#	$dotos .= " $name.o ";
#}

my %dmrgMain = (name => 'mpspp', dotos => $dotos);

push @drivers,\%dmrgMain;

createMakefile();

sub createMakefile
{
#	unlink("Engine/Version.h");
	Make::backupMakefile();
	if (!(-r "Config.make")) {
		my $cmd = "cp Config.make.sample Config.make";
		system($cmd);
		print STDERR "$0: Executed $cmd\n";
	}

	my $fh;
	open($fh,">Makefile") or die "Cannot open Makefile for writing: $!\n";

	my %args;
	$args{"code"} = "MPS++";
	Make::newMake($fh,\@drivers,\%args);
	print STDERR "File Makefile has been written\n";
}

