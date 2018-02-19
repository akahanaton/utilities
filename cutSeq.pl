#!/usr/bin/perl
#*Copyright 2005  BGI
#*All right reserved
#*Filename:  cutSeq.pl
#*Description:  cut the segment from the parent sequence
#*Version : 1.0
#*Progremmer: Li ruiqiang(lirq@genomics.org.cn)
#*Time:    2005.04.20
#

use Getopt::Long;

#**** USAGE *************************************
my $usage=<<"USAGE";
Usage : echo new_ID  Seq_ID  Start   Length | $0 Input.seq  STDOUT

notes:

USAGE
die $usage if (@ARGV<1);
#************************************************

open (Seq, $ARGV[0]) ||die;

my %S;
my $id;

while (<Seq>) {
	chomp;
	if (/^>(\S+)/) {
		$id = $1;
	}
	else {
		$S{$id} .= $_;
	}
}

#
my $out = 60;

while (<STDIN>) {
	chomp;
	my @a = split(/\s+/, $_);
	my $i = substr($S{$a[1]}, $a[2], $a[3]);
	print ">$a[0]\n";
	my $t=0;
	my $len = length $i;
	while ($len > $t) {
		my $e=substr($i, $t, $out);
		print "$e\n";
		$t += $out;
	}
}
