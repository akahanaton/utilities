#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use MyMod::Bio::Tools::NW;
use Bio::SeqIO;

my %opts;
my $usage = qq{
	Usage: $0 <fa.file> <length_cutoff>
};

getopts('', \%opts);
die($usage) if (@ARGV != 2);
my ($fas, $lenCut) = @ARGV;

my %seqInfo;
my $seqin=Bio::SeqIO->new(-file=>$fas,-format=>"fasta");
while(my $seq=$seqin->next_seq()){
	my @arr = split(/X+/, $seq->seq());
	foreach my $tmpSeq ( @arr){
		if (length $tmpSeq > $lenCut){
			$seqInfo{$tmpSeq} = $seq->id();
		}
	}
}

foreach my $seq(sort keys %seqInfo){
	print ">$seqInfo{$seq}\n";
	print "$seq\n";
}

