#!/usr/bin/perl -w

my $t_len = 0;
my @lens  = ();

my $file = shift or die("Usgae: $0 fasta_seq_file min_length percent\n");
my $min_size = shift || 100;
my $percent = shift || 50;

open(IN, $file) or die("failed to open file $file\n");
my $len = 0;
my $seq_num = 0;
my $mean = 0;
while(<IN>){
	if(/^>/){
        $seq_num += 1;
		if($len >= $min_size){
			push(@lens, $len);
			$t_len += $len;
		}
		$len = 0;
	} else {
		chomp;
		$len += length($_);
	}
}
if($len >= $min_size){
	push(@lens, $len);
	$t_len += $len;
}
close IN;

@lens = sort {$b <=> $a} @lens;

$mean = $t_len / $seq_num;

my $sum = $t_len/100*$percent;
my $nSome;
foreach my $l (@lens){
	$s += $l;
	if($s>=$sum){
		$nSome = $l;
		last;
	}
}

print "TotalBase= $t_len TotalSeq= $seq_num len_cutoff= $min_size Longest= $lens[0] Shortest= $lens[-1] N$percent= $nSome Mean= $mean\n";
