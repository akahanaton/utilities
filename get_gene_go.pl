#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <gene.id> <go_file>
	Assume first colomne is the reads' name
};

getopts('', \%opts);
die($usage) if (@ARGV != 2);
my ($gene_id, $go_file) = @ARGV;

my %go;
open(GO, $go_file) || die("** Fail to open '$go_file'.\n");
while ( my $line = <GO>){
	my @arr = split(/\s+/, $line);
	$go{$arr[0]} = $line;
}
close GO;

open(Gene, $gene_id) || die("** Fail to open '$gene_id'.\n");
while( my $line = <Gene>)
{
	my @arr = split(/\s+/, $line);
	if(defined($go{$arr[0]})){
		print $go{$arr[0]};
	}else{
		printf "%s\n", $arr[0];
	}
}
close Gene;
