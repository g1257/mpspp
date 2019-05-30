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
my %mpspp = (name => 'mpspp', dotos => "mpspp.o Provenance.o ProgramGlobals.o");
my %mpopp = (name => 'mpopp', dotos => "mpopp.o Provenance.o ProgramGlobals.o");
my @drivers = (\%provenanceDriver,\%progGlobalsDriver, \%mpspp, \%mpopp);

createMakefile();

sub createMakefile
{
	Make::backupMakefile();
	if (!(-r "Config.make")) {
		my $cmd = "cp Config.make.sample Config.make";
		system($cmd);
		print STDERR "$0: Executed $cmd\n";
	}

	my $fh;
	open($fh, ">", "Makefile") or die "Cannot open Makefile for writing: $!\n";

	my %args;
	$args{"code"} = "MPS++";
	Make::newMake($fh,\@drivers,\%args);
	print STDERR "File Makefile has been written\n";
}

