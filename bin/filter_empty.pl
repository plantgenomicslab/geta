#!/opt/miniconda3/bin/perl
use strict;
use warnings;
$/="\n>";
while (<>) {
	s/>//g;
	my ($id, $seq) = split (/\n/, $_);
	print ">$_" if (length $seq)
}
