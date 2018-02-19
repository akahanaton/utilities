#!/usr/bin/perl
use strict;
# use warnings;
use Carp;

###########################################################
# Rongkun Shen
# Mockler lab
# Oregon State University
# Corvallis, OR
###########################################################

# Usage
my $usage = "
RGA_eland - Reference Guided Assember using ELAND results. It does:
  1. calcuate the real coverage of Solexa microreads to the reference sequence(s)
  2. report the coordinates of all the maximum contigs and length
  3. output density graphes of chromosome size (optional)
  
  Note: only takes microreads with the unique hits [U0, U1, U2] from ELAND results
	
Usage: perl $0 options
 required:
  -r	Reference size file [only the SIZE of each chromosome, e.g.,
		chr1[TAB]29231457  
		[NOTE: [TAB] means tab-delimited; chromosome names should be the same as that in ELAND file]
  -l	ELAND result file
 optional:
  -m	Mismatches: a string of U0,U1,U2without space, e.g., U0, U0U1, U1U2, U0U1U2 [default: U0]
  -p	Switch of Plotting or not [default: not plot; need gnuplot version 4.2.2 or above]
  
";

use Getopt::Std;
use vars qw( $opt_r $opt_l $opt_m $opt_p);
# command line processing #
getopts('r:l:m:p');
die $usage unless ($opt_r and $opt_l);
my ($ref_file, $eland_file, $mismatch );

my $pid = getpgrp(0);

$ref_file = $opt_r if $opt_r;
open (REF,   "$ref_file") or die "$!, $ref_file not exists: exiting\n";
$eland_file = $opt_l if $opt_l;
open (ELAND, "$eland_file") or die "$!, $eland_file not exists: exiting\n";
$mismatch	= $opt_m ? $opt_m : 'U0'; 

my %chr_hash = ();
my $chr;

# print "Load reference sequence sizes ...\n";
while (<REF>) {
	chomp;
	my ($chr, $size) = split (/\t/, $_);
	$chr_hash{$chr} = $size;
	# print "$chr\t$size\n";
}
close REF;
# print "\tDone\n";


my %eland_chr_hash = ();  # for categorize each chromosome - memory intensive
my %assembly_coord_array_hash = ();
my $chr_size = 0;
my $first_read_flag = 1;

print "Load ELAND result file...\n";
while (<ELAND>) {
	chomp;
	# s/\>//g;
	my $line = $_;
	my @item_array = split(/\t/, $line);
	
	#
	# column definitions from ELAND result:
	#
	# $item_array[0]: read_head	
	# $item_array[1]: sequence	
	# $item_array[2]: align_type	
	# $item_array[3]: #_of_0-MisMatch	
	# $item_array[4]: #_of_1-MM	
	# $item_array[5]: #_of_2-MM	
	# $item_array[6]: chr	
	# $item_array[7]: coord	
	# $item_array[8]: strand	
	# $item_array[9]: base substitutions
	
	# push all the selected lines into hash
	if ( $mismatch =~ /U0/ and ("$item_array[2]" eq 'U0') ) { 		
		push ( @{ $eland_chr_hash{$item_array[6]} }, $line);
	}
	
	if ( $mismatch =~ /U1/ and ("$item_array[2]" eq 'U1') ) { 		
		push ( @{ $eland_chr_hash{$item_array[6]} }, $line);
	}
	if ( $mismatch =~ /U2/ and ( "$item_array[2]" eq 'U2') ) { 		
		push ( @{ $eland_chr_hash{$item_array[6]} }, $line);
	}
}
close ELAND;
print "\tDone\n";

# coverage output file
my $coverage_file = $eland_file.'_'.$mismatch.'.RGA.txt';
unlink $coverage_file if (-e $coverage_file);
open (COVR, ">$coverage_file") or die "$!: exiting\n";

# process each reference sequence
foreach my $key ( sort {$a cmp $b} keys %eland_chr_hash ) {  # $key here is chromosome
	$chr_size = $chr_hash{$key};
	print ">$key\t$chr_size\n";
	print COVR ">$key\t$chr_size\n";
	my $density  = $eland_file.'_'.$mismatch.'_'.$key."_density.txt";
	my $read_length = 0;
	my @array_line = @{ $eland_chr_hash{$key} }; 
	my %density_coord_array_hash = ();
	foreach my $line (@array_line) {
		# print "$line\n";
		my @input = split(/\t/, $line);
		my ($start, $end);
		
		$read_length = length($input[1]);
		$start = $input[7]; # start position
		$end = $start + $read_length - 1;
		push ( @{$density_coord_array_hash{$start}}, $end ); # for density output
		if ( !$assembly_coord_array_hash{$start} or ( $assembly_coord_array_hash{$start} and $end > $assembly_coord_array_hash{$start} ) ) {
			$assembly_coord_array_hash{$start} = $end; # for assembly
		}
	}
	
	my @density_array = ();
	foreach my $key3 ( sort {$a <=> $b} keys %density_coord_array_hash) { #numeric sort keys
		foreach my $val3 ( sort {$a <=> $b} @{$density_coord_array_hash{$key3}} ) {
			for ( my $i = $key3; $i <= $val3; $i ++ ) {
				$density_array[$i] ++;
			}
		}
	}
	
	# output chromosome-size density data
	$density =~ s/[\|\:\;\,\>\*\^\&\?\#\$\@\!\~\s]/\_/g;
	$density =~ s/\_+/\_/g;
	my $block = $density;
	$block =~ s/density/block/g;
	unlink $density if (-e $density);
	unlink $block if (-e $block);
	open (DENSITY, ">$density") or die "$!: exiting\n";
	open (BLOCK,   ">$block")   or die "$!: exiting\n";

	my $block_start = 1;
	foreach (my $j =1; $j <= $chr_size; $j ++)  {
		if ( ($density_array[$j] != $density_array[$j+1]) or ($j == $chr_size) ) {
			if (!$density_array[$j]) {
				$density_array[$j] = 0;
			}
			# print "$block_start\t$j\t$density_array[$j]\n"; # for testing
			print BLOCK "$block_start\t$j\t$density_array[$j]\n";
			$block_start = $j + 1;
		} 
		if ( !$density_array[$j] ) { 
			print DENSITY "0\n";
			# print "0\n"; # for testing
		} else {
			print DENSITY "$density_array[$j]\n";
			# print "$density_array[$j]\n"; # for testing
		}
	}
	@density_array = ();
	close DENSITY;
	close BLOCK;
	
	my ($min, $max, $key1, $val1);
	my %assembly_coord_array_hash2 = ();
	foreach my $key1 ( sort {$a <=> $b} keys %assembly_coord_array_hash ) { #numeric sort keys
		$val1 = $assembly_coord_array_hash{$key1};
	
		# print "key_val\t$key1\t$val1\n"; # for testing only
		if ( !$min ) { # for 1st set of coordinates
			$min = $key1;
			$max = $val1;
		} else {
			if ($key1 > $max + 1) {
			$assembly_coord_array_hash2{$min} = $max; # store larger contig coordinates in a new hash2
			$min = $key1;
			$max = $val1;
			} elsif ($val1 > $max) {
				$max = $val1;
				$assembly_coord_array_hash2{$min} = $max;
			}
		}
		# print "min_max\t$min\t$max\n"; # for testing only
	}

	%assembly_coord_array_hash = ();
	my $sum = 0; 
	my $len = 0;
	my $cnt = 0;
	
	# output stored larger contig coordinates and calculate mroe info
	foreach my $key2 ( sort {$a <=> $b} keys %assembly_coord_array_hash2 ) {
		my $val2 = $assembly_coord_array_hash2{$key2};
		$len = $val2 - $key2 + 1;
		$sum += $len;
		# print "$key2\t$val2\t$len\n";
		print COVR "$key2\t$val2\t$len\n";
		$cnt ++;
	}
	
	print "\t# of contigs= " . $cnt . "\n\tSUM= $sum\n\tChromosome size= ". $chr_size."\n\tCoverage= ". $sum/$chr_size."\n\n";
	print COVR "\t# of contigs= " . $cnt . "\n\tSUM= $sum\n\tChromosome size= ". $chr_size."\n\tCoverage= ". $sum/$chr_size."\n\n";
	
	# plot or not using gnuplot
	if ($opt_p) { 
		my $pltFile = $density.'.'.$pid.'temp.plt';
		unlink $pltFile if (-e $pltFile);
		open (PNG_OUT, ">$pltFile") or die "$!: exiting\n";
		print PNG_OUT "set output \"$density.png\"\n";
		if ( $chr_size >= 50000 ) {
			print PNG_OUT "set term png size 50000, 1000 \n";
		} else {
			print PNG_OUT "set term png size ". $chr_size .", 1000 \n";
		}
		print PNG_OUT "set xlabel \"Seq_base\"\n";
		print PNG_OUT "set ylabel \"Count\"\n";

		# print PNG_OUT "set title \"$density\"\n";
		print PNG_OUT "plot \"$density\" notitle w i\n";
		close PNG_OUT;

		`gnuplot $pltFile`;
		unlink $pltFile;
	}
}
close COVR;

# `rm -f temp.R`;
print "Job finished!\n";
`rm -rf RGA_density_job_finished`;
`touch RGA_density_job_finished`;
exit;


1;

__END__

=head1 NAME

RGA - Reference Guided Assembler

=head1 SYNOPSIS

RGA_eland.pl options

=head1 Description

=head2 Features

1. Calcuate the coverage of Solexa microreads to the ref genome

2. Report the coordinates of all the maximum contigs and length

3. Reference(scaffold)-guided assembly

4. Density output with an array of chromosome size and a block file

=head2 ELAND result format:

read_head sequence    align_type  #_of_0-MisMatch #_of_1-MM   #_of_2-MM   chr coord   strand  substitutions

>HWI-EAS153_3_2008FAAXX_3_1_639_891     GTTTGGACGCAGATTTGTCAGGCTCTCAGATT        U1      0       1       0       chr1 151717413       R       ..      2C

>HWI-EAS153_3_2008FAAXX_3_1_110_477     TGATTCGAGGCTTCCTTGAATTGTAGCACTTA        NM      0       0       0

>HWI-EAS153_3_2008FAAXX_3_1_101_555     TTCTTCTTTGACAATTTGTTATATTCACCCAA        U1      0       1       0       chr9 4822003 F       ..      19G

>HWI-EAS153_3_2008FAAXX_3_1_633_916     GCCCTTCTTAAAGATCCCTGTCTTCATCATGA        U0      1       0       0       chr15        44037561        F       ..


=head1 AUTHOR

Rongkun Shen (rongkun.shen@gmail.com)

=head1 ACKNOWLEDGEMENTS

This software was developed in

Mockler lab

Department of Botany and Plant Pathology

Oregon State Univeristy, Corvallis, OR.

=head1 COPYRIGHT

Copyright (C) 2007-2008 Rongkun Shen. All Rights Reserved.

=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=cut

