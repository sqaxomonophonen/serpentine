#!/usr/bin/env perl
use strict;
use warnings;

my $usage = "usage: $0 <name> <outfile> <vert infile> <frag infile>\n";

$ARGV[2] or die($usage);

my $name = shift;
my $outfile = shift;

my $out = "";

for my $type ("vert", "frag") {
	my $infile = shift;
	-e $infile or die("no such file: $infile\n");

	$out .= "static const char* ${name}_${type}_src =\n";

	open IN, $infile;
	while (<IN>) {
		chop;
		$out .= "\t\"$_\\n\"\n";
	}
	$out .= ";\n";
	close IN;
}


open OUT, ">$outfile";
print OUT $out;
close OUT;
