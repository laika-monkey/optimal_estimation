#!/usr/bin/perl

#_ or qdel -u wsessions

use strict;

chomp(my @jobs = `qstat`);

for my $job (@jobs) {
	my @cols = split(/\s+/, $job);

	my $pid = @cols[1];
	my $scr = @cols[3];

	system "qdel $pid " if ($scr eq 'run_lbldis'); 
}
