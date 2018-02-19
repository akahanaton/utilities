#!/usr/bin/perl
use strict;
use Carp;
use Getopt::Std;
use vars qw( $opt_r $opt_i $opt_b $opt_e $opt_t $opt_s $opt_p);

# Usage
my $usage = "
RGA - Reference Guided Assember. It does:
  1. calcuate the real coverage of Solexa microreads to the reference sequence(s)
  2. report the coordinates of all the maximum contigs and length
  3. output density graphes of chromosome size

Usage: perl $0 options
 required:
  -r	reference sequence(s) in a single FASTA file
  -i	input file (sequences) in a single FASTA file
 optional:
  -e	error_rate [optional, by default: 0.00]; calculate: (# of mismatches) / (read length)
  -b	BLAT results file [optional]
  -t	BLAT tileSize [optional, range in [6, 15], by default: 10; if -b is provided, this will be ignored]
  -s	BLAT stepSize [optional, range in [1, tileSize], by default: 1; if -b is provided, this will be ignored]
  -p	switch of plotting or not [by default: not plot] (need gnuplot v4.2.2 or later in search path)
  
";

# command line processing 
getopts('r:i:b:e:t:s:p');
die $usage unless ($opt_r and $opt_i);

my ($ref, $src, $blat, $error_rate, $step_size, $tile_size);

$ref	= $opt_r if $opt_r;
$src	= $opt_i if $opt_i;
$blat   = $opt_b if $opt_b;
$tile_size	= $opt_t ? $opt_t : 10; 
$step_size	= $opt_s ? $opt_s : 1; 
$error_rate	= $opt_e ? $opt_e : 0.00;

my $pid = getpgrp(0);
my $blat_cwd;

unless ($blat) {
	my @src_name_array;
	if ($src =~ /\//) {
		@src_name_array = split(/\//, $src);
		$blat = $src_name_array[$#src_name_array] . ".blatI0G0T10S$step_size.pslx";
	} else {
		$blat = $src. ".blatI0G0T10S$step_size.pslx";
	}

	print "\nBLAT...... \n!!! Caution: It may take long time (upto weeks!!!) if both reference and source files are large !!! \n";
	# print "blat -maxGap=0 -maxIntron=0 -stepSize=$step_size -minScore=18 -noHead -out=pslx -tileSize=$tile_size -ooc=/home/bpp/shenr/bin/10.ooc $ref $src $blat";
	`blat -maxGap=0 -maxIntron=0 -stepSize=$step_size -minScore=18 -noHead -fine -out=pslx -tileSize=$tile_size $ref $src $blat`;
	print "\tDone\n";
	$blat_cwd = $blat;
} else {
	my @blat_name_array;
	if ($blat =~ /\//) {
		@blat_name_array = split(/\//, $blat);
		$blat_cwd = $blat_name_array[$#blat_name_array];
	} else {
		$blat_cwd = $blat;
	}
}

# print "blat $blat_cwd\n"; # for testing

my %blat_line_hash = ();
my %coord_hash1 = ();
my $chr_size = 0;
open (BLAT, "$blat") or die "$!: $blat\n";

print "\nLoad BLAT results file...\n";
while (<BLAT>) {
	chomp;
	if (/^\d+?\t\d+?\t.*/) { # exclude first several lines of column definitions
		my @blat_field_array = split(/\t/, $_);
		# 10-31-2007 modified to reduce memory usage by > 60%
		
		# $blat_field_array[13] is T name (chr)
		# $blat_field_array[0]  is match
		# $blat_field_array[10] is Q size
		# $blat_field_array[11] is Q start
		# $blat_field_array[12] is Q end
		# $blat_field_array[14] is T size
		# $blat_field_array[15] is T start
		# $blat_field_array[16] is T end
		
		if ( ($blat_field_array[10]-$blat_field_array[0])/$blat_field_array[10] <= $error_rate ) { # for control of mismatches
			push ( @{ $blat_line_hash{$blat_field_array[13]} }, "$blat_field_array[0]\t$blat_field_array[10]\t$blat_field_array[11]\t$blat_field_array[12]\t$blat_field_array[14]\t$blat_field_array[15]\t$blat_field_array[16]");
			
			# in blat	0	10	11	12	14	15	16
			# in hash	0	1	2	3	4	5	6
		}
	}
}
close BLAT;
print "\tDone\n\n";

my $coverage_file = $blat_cwd.'.Err_'.$error_rate.'.RGA.'.$pid.'.txt';
# print "cov $coverage_file\n"; # for testing
unlink $coverage_file if (-e $coverage_file);
open (COVR, ">$coverage_file") or die "$!: $coverage_file\n";
	
foreach my $chr ( sort {$a cmp $b} keys %blat_line_hash ) {
	print ">$chr\n";
	print COVR ">$chr\n";
	my $density = $blat_cwd.'.Err_'.$error_rate.'.'.$chr.'.density.'.$pid.'.txt';
	my @line_array = @{ $blat_line_hash{$chr} }; 
	my %coord_array_hash = ();
	foreach my $line (@line_array) {
		my @hash_field_array = split(/\t/, $line);
		$chr_size = $hash_field_array[4];
		my ($start, $end);
		
		# BLAT T_start + 1 = BLAST T_start; but T_end is correct
		$start = ($hash_field_array[5] + 1) - $hash_field_array[2];
		$end = $start + $hash_field_array[1] - 1;
		# print "density $start\t$end\n";
		# 10-31-2007 fixed 
		if ($start <= 0) {
			# print "$_\n";
			# print "Wrong\t$start\t$end\t";
			$start = 1;
			
			# to solve those cases
			# 31      0       0       0       0       0       0       0       +       4529813|13220972-13220939|length=32     32      1       32      super_69        8507    0       31      1      31,      1,      0,      cgggacccgaaagatggtgaactatgcctga,        cgggacccgaaagatggtgaactatgcctga,
			# Wrong   0       31      Wrong2  31      0       1       32
			# 30      0       0       0       0       0       0       0       +       4538399|4421042-4421009|length=32       32      2       32      super_337       3590    0       30      1      30,      2,      0,      cgagaggaaccgttgattcacacaattggt, cgagaggaaccgttgattcacacaattggt,
			# Wrong   -1      30      Wrong2  30      0       2       32
		}
		if ($end > $hash_field_array[4]) {
			# print "$_\n";
			# print "Wrong\t$start\t$end\t";
			$end = $hash_field_array[4];
			
			# to solve those cases
			# 31      0       0       0       1       1       0       0       -       4511393|2510083-2510050|length=32       32      0       32      super_209       4301    4270    4301    2      24,7,    0,25,   4270,4294,      ggatagtagacagggacagtggga,tctcgtt,       ggatagtagacagggacagtggga,tctcgtt,
			# Wrong   4271    4302    Wrong2  31      0       0       32
			# 30      0       0       0       0       0       0       0       +       4523014|15757086-15757053|length=32     32      0       30      super_283       3559    3529    3559    1      30,      0,      3529,   caagtcagacgaacgatttgcacgtcagta, caagtcagacgaacgatttgcacgtcagta,
			# Wrong   3530    3561    Wrong2  30      0       0       30
			
		}
		# print "$start\t$end\n";
		if ( $start < $end ) {
			push ( @{$coord_array_hash{$start}}, $end ); # for density output
			if ( !$coord_hash1{$start} || ( $coord_hash1{$start} && $end > $coord_hash1{$start} ) ) {
				$coord_hash1{$start} = $end; # for assembly
				# print "$hash_field_array[15]\t$hash_field_array[16]\t";
			}
		}
	}
	# print "coord_hash1 done\n"; # for testing only
	
	my @density_array = ();
	foreach my $key0 ( sort {$a <=> $b} keys %coord_array_hash) { 
		foreach my $val0 ( sort {$a <=> $b} @{$coord_array_hash{$key0}} ) {
			for ( my $i = $key0; $i <= $val0; $i ++ ) {
				$density_array[$i] ++;
			}
		}
	}
	
	# output chromosome-size density data
	$density =~ s/[\|\:\;\,\>\*\^\&\?\#\$\@\!\~\s]/\_/g;
	$density =~ s/\_+/\_/g;
	unlink $density if (-e $density);
	open (DENSITY, ">$density") or die "$!: $density\n";
	foreach (my $j =1; $j <= $chr_size; $j ++)  {
		if ( !$density_array[$j] ) { 
			# print DENSITY "$j\t0\n";
			print DENSITY "0\n";
			# print "$j\t0\n";
		} else {
			print DENSITY "$density_array[$j]\n";
			# print DENSITY "$j\t$density_array[$j]\n";
			# print "$j\t$density_array[$j]\n";
		}
	}
	@density_array = ();
	close DENSITY;
	
	my ($min, $max, $key1, $val1);
	my %coord_hash2 = ();
	foreach my $key1 ( sort {$a <=> $b} keys %coord_hash1 ) { #numeric sort keys
		
		$val1 = $coord_hash1{$key1};
	
		# print "key_val\t$key1\t$val1\n"; # for testing only
		if ( !$min ) { # for 1st set of coordinates
			$min = $key1;
			$max = $val1;
		} else {
			if ($key1 > $max + 1) {
			$coord_hash2{$min} = $max; # store larger contig coordinates in a new coord_hash2
			$min = $key1;
			$max = $val1;
			} elsif ($val1 > $max) {
				$max = $val1;
				$coord_hash2{$min} = $max;
			}
		}
		# print "min_max\t$min\t$max\n"; # for testing only
	}

	%coord_hash1 = ();
	my $sum = 0; 
	my $len = 0;
	my $cnt = 0;
	
	# output stored larger contig coordinates and calculate mroe info
	foreach my $key2 ( sort {$a <=> $b} keys %coord_hash2 ) {
		my $val2 = $coord_hash2{$key2};
		$len = $val2 - $key2 + 1;
		$sum += $len;
		# print "$key2\t$val2\t$len\n";
		print COVR "$key2\t$val2\t$len\n";
		$cnt ++;
	}
	
	print "\t# of contigs= " . $cnt . "\n\tSUM= $sum\n\tChromosome size= ". $chr_size."\n\tCoverage= ". $sum/$chr_size."\n\n";
	print COVR "\t# of contigs= " . $cnt . "\n\tSUM= $sum\n\tChromosome size= ". $chr_size."\n\tCoverage= ". $sum/$chr_size."\n\n";
	
	if ($opt_p) {
		my $pltFile = $density.'.'.$pid.'temp.plt';
		unlink $pltFile if (-e $pltFile);
		open (PNG_OUT, ">$pltFile") or die "$!: $pltFile\n";
		print PNG_OUT "set output \"$density.png\"\n";
		if ( $chr_size >= 10000 ) {
			print PNG_OUT "set term png size 10000, 1000 \n";
		} else {
			print PNG_OUT "set term png size ". $chr_size .", 1000 \n";
		}
		print PNG_OUT "set xlabel \"Base Position (bp)\"\n";
		print PNG_OUT "set ylabel \"Base Density\"\n";

		# print PNG_OUT "set title \"$density\"\n";
		print PNG_OUT "plot \"$density\" w i\n";
		close PNG_OUT;

		`gnuplot422 $pltFile`;
		unlink $pltFile;
		
		open (PNG_OUT, ">$pltFile") or die "$!: exiting\n";
		if ( $chr_size >= 10000 ) {
			print PNG_OUT "set term png size 10000, 1000 \n";
		} else {
			print PNG_OUT "set term png size ". $chr_size .", 1000 \n";
		}
		print PNG_OUT "set output \"$density.y100.png\"\n";
		print PNG_OUT "set term png \n";

		print PNG_OUT "set xlabel \"Base Position (bp)\"\n";
		print PNG_OUT "set ylabel \"Base Density\"\n";
		print PNG_OUT "set yrange [0:100]\n";

		# print PNG_OUT "set title \"$density\"\n";
		print PNG_OUT "plot \"$density\" w i\n";
		close PNG_OUT;

		`gnuplot422 $pltFile`;
		#### only run on new machines
		unlink $pltFile;
	}

}
close COVR;

print "RGA job $pid finished!\n";
`rm -rf RGA_job.$pid.Err$error_rate.finished`;
`touch  RGA_job.$pid.Err$error_rate.finished`;

1;

__END__

=head1 NAME

RGA - Reference Guided Assembler

=head1 SYNOPSIS

RGA.pl Reference Microread_file other_options

=head1 Description

RGA aligns input sequences to their best match in reference sequences and creates guided consensus sequences from overlapping alignments.

=head2 Features

1. Calcuate the coverage of Solexa microreads to the ref genome

2. Report the coordinates of all the maximum contigs and length

3. Output density graphes of chromosome size


BLAT result file psl format

match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts  tStarts
        match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count

=head1 Some issues

1. BLAT T_start + 1 = BLAST T_start 
(kind of bug in BLAT!)
For BLAT data:  $start should be = T_start + 1;
but BLAT T-end is fine.

2. We should use $val > $max + 1


error rate control:

(Q_size - match + mismatch)/Q_size

for 25mers

match		mismatch	Q_size	error_rate

25		0			25		0.00 (25-25)/25

cut-----------------

24		0			25		0.04 (25-24)/25

24		1			25		0.04 (25-24)/25

cut-----------------

23		0			25		0.08 (25-23)/25

23		1			25		0.08 (25-23)/25

23		2			25		0.08 (25-23)/25

cut-----------------


(Q_size - match + mismatch)/Q_size

for 32mers

match		mismatch	Q_size	error_rate

32		0			32		0.00000 (32-32)/32

cut-----------------

31		0			32		0.03125 (32-31)/32 #mismatch at ends

31		1			32		0.03125 (32-31)/32 #mismatch at ends and in the middle

cut-----------------

30		0			32		0.06250 (32-30)/32

30		1			32		0.06250 (32-30)/32

30		2			32		0.06250 (32-30)/32

cut-----------------

29		0			32		0.09375 (32-29)/32

29		1			32		0.09375 (32-29)/32

29		2			32		0.09375 (32-29)/32

29		3			32		0.09375 (32-29)/32

cut-----------------


03:28 PM 10-31-2007

found a bug in using the coordinates about the following cases since BLAT just reports T_start and T_end of part of 32mer

fixed 03:41 PM 10-31-2007

31      0       0       0       0       0       0       0       +       37|FC20|1_1|164_853|1_32        32      0       31      super_4 26695973        13194666        13194697        1       31,     0,      13194666,       gttatgacttagatacttgtggtagcatccc,        gttatgacttagatacttgtggtagcatccc,

30      0       0       0       0       0       0       0       +       28|FC20|1_1|62_704|1_32 32      2       32      super_4 26695973        25845676        25845706        1       30,     2,      25845676,       gggattcaagaatgtgaagaacctggaagg, gggattcaagaatgtgaagaacctggaagg,

30      1       0       0       0       0       0       0       -       248|FC20|1_1|993_611|1_32       32      1       32      super_7 17724951        4778553 4778584 1       31,     0,      4778553,        ttgccatctccatatacttgaactttctctt,        ttcccatctccatatacttgaactttctctt,


but no problem with the following cases

31      1       0       0       0       0       0       0       -       53|FC20|1_1|506_938|1_32        32      0       32      super_5 20333241        3538045 3538077 1       32,     0,      3538045,        tttgtaacttgttacagcttcattgcctatta,       tttgtaacttggtacagcttcattgcctatta,

30      2       0       0       0       0       0       0       -       129|FC20|1_1|574_792|1_32       32      0       32      super_6 19186327        7977727 7977759 1       32,     0,      7977727,        caggcacggttccaatgccaccaatcttgtac,       caggcacagttccaatgccaccgatcttgtac,


=head1 AUTHOR

Rongkun Shen (rongkun.shen@gmail.com)

=head1 ACKNOWLEDGEMENTS

This software was developed in the Department of Botany and Plant Pathology at 
Oregon State Univeristy, Corvallis, OR.

=head1 COPYRIGHT

Copyright (C) 2007-2008 Rongkun Shen. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut
