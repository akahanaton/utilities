#!/usr/bin/perl -w
#
# Author: Ruan Jue
#
use warnings;
use strict;

my $scaff_file = shift or die;
my $blat_file  = shift or die;

my %contigs = ();

print STDERR "Loading alignment file ...\n";

open(IN, $blat_file) or die($!);
while(<IN>){
	my @tabs = split;
	next unless(@tabs == 21);
	next if($tabs[0] < 0.95 * $tabs[10]);
	next if($tabs[5] + $tabs[7] > 20);
	my ($ref, $strand, $c_start, $c_end);
	if($tabs[8] eq '+'){
		$c_start = $tabs[15] - $tabs[11];
		$c_end   = $tabs[16] + ($tabs[10] - $tabs[12]);
	} else {
		$c_start = $tabs[15] - ($tabs[10] - $tabs[12]);
		$c_end   = $tabs[16] + $tabs[11];
	}
	$strand = $tabs[8];
	$ref    = $tabs[13];
	push(@{$contigs{$tabs[9]}}, [$ref, $strand, $c_start, $c_end]);
}
close IN;

=pod

print STDERR "Analyse contigs' positions on reference ...\n";

foreach my $key (keys %contigs){
	my $ctg = $contigs{$key};
	next if(@$ctg == 1);
	my @loci = sort {$a->[2] <=> $b->[2]} @$ctg; 
	$contigs{$key} =\@loci;
}

=cut

print STDERR "Checking scaffolds ...\n";

my ($out_pass, $out_fail);

open(IN, $scaff_file) or die($!);
open($out_pass, ">$scaff_file.pass") or die;
open($out_fail, ">$scaff_file.fail") or die;

my $scf_name;
my @ctgs;
my @regions = ();
while(<IN>){
	if(/^>(\S+)/){
		&evaluate_scaff($scf_name, \@ctgs, \%contigs, \@regions, $out_pass, $out_fail) if(@ctgs);
		$scf_name = $1;
		@ctgs = ();
	} else {
		my @tabs = split;
		push(@ctgs, \@tabs);
	}
}
&evaluate_scaff($scf_name, \@ctgs, \%contigs, \@regions, $out_pass, $out_fail) if(@ctgs);
close IN;
close $out_pass;
close $out_fail;

print STDERR "Program terminated normally.\n";

1;

sub evaluate_scaff {
	my ($scf_name, $ctgs, $loci_hash, $regions, $out_pass, $out_fail) = @_;
	my $scf_len = $ctgs->[-1]->[1] + $ctgs->[-1]->[3];
	my $matrix = [];
	my $dists  = [];
	my $pos    = 0;
	my $ctgs2  = [];
	foreach my $c (@$ctgs){
		my $ctg = $loci_hash->{$c->[0]};
		next unless(defined $ctg);
		push(@$matrix, $ctg);
		push(@$ctgs2, $c);
		push(@$dists, $c->[1] - $pos - 1);
		$pos = $c->[3] + $c->[1];
	}
	# Dynamic programming
	my $best_hit   = 0;
	my $best_score = 0;
	for(my $i=0;$i<@$matrix;$i++){
		my $col = $matrix->[$i];
		if($i == 0){
			foreach my $cell (@$col) {
				$cell->[4] = 0;
				$cell->[5] = -1;
			}
			next;
		}
		my $last = $matrix->[$i - 1];
		for(my $j=0;$j<@$col;$j++){
			my $cell = $col->[$j];
			my ($max_score, $idx);
			$max_score = -9999999;
			$idx       = 0;
			for(my $k=0;$k<@$last;$k++){
				my $lc = $last->[$k];
				my $score = 0;
				if($cell->[0] eq $lc->[0]){
					$score = - abs((($ctgs2->[$i]->[2] eq $cell->[1])? ($cell->[2] - $lc->[3] - 1):($lc->[2] - $cell->[3] - 1)) - $dists->[$i]);
					$score = -1000 if($score < -1000);
				} else {
					$score = -1000;
				}
				$score += $lc->[4];
				if($score > $max_score){
					$max_score = $score;
					if($score == -1000){
						$idx = - $k - 1;
					} else {
						$idx = $k;
					}
				}
				$cell->[4] = $max_score;
				$cell->[5] = $idx;
			}
		}
		if($i == @$matrix - 1){
			my $max_score = -999999;
			my $idx       = 0;
			for(my $j=0;$j<@$col;$j++){
				my $cell = $col->[$j];
				if($cell->[4] > $max_score){
					$max_score = $cell->[4];
					$idx = $j;
				}
			}
			$best_hit   = $idx;
			$best_score = $max_score;
		}
	}
	# Track back
	my %hash = ();
	my $mis = 0;
	my $n_mis = 0;
	for(my $i=@$matrix-1;$i>=0;$i--){
		if($best_hit < 0){
			$best_hit = -1 - $best_hit;
			$n_mis ++;
			$mis += $ctgs->[$i]->[3];
		}
		$hash{$ctgs2->[$i]->[0]} = $matrix->[$i]->[$best_hit];
		$best_hit = $matrix->[$i]->[$best_hit]->[5];
	}
	my $output;
	$best_score = -$best_score - $n_mis * 1000;
	#if($mis >= $scf_len * 0.05 || $best_score/$scf_len>0.1*scalar(@$ctgs)){
	if($mis >= $scf_len * 0.05){
		$output = $out_fail;
	} else {
		$output = $out_pass;
	}
	my $sum = 0;
	my $num = 0;
	my $ref = '';
	foreach my $c (@$ctgs){
		my $corr = $hash{$c->[0]};
		if($corr){
			if($num == 0){
				$ref = $corr->[0];
				$sum = $corr->[2];
			} else {
				$sum = $corr->[2] if($sum > $corr->[2]);
			}
			$num ++;
		}
	}
	if($num<1){
		return;
	}
	# Print
	print $output ">$scf_name\tlength=$scf_len\tn_ctg=", scalar(@$ctgs), "\tdeviate=$best_score\tn_mis=$n_mis\tmis_len=$mis\n";
	$sum = 0;
	$num = 0;
	$ref = '';
	foreach my $c (@$ctgs){
		my $corr = $hash{$c->[0]};
		if($corr){
			if($num == 0){
				$ref = $corr->[0];
				$sum = $corr->[2];
			} else {
				$sum = $corr->[2] if($sum > $corr->[2]);
			}
			$num ++;
			print $output "$corr->[0] $corr->[1] $corr->[2] ~ $corr->[3]";
		} else {
			print $output "No hit";
		}
		print $output "\t#";
		print $output join("\t", @$c);
		print $output "\n";
	}
	if($num > 1 and $scf_len > 200){ 
		push(@$regions, [$ref, $scf_name, $sum, $sum + $scf_len, $scf_len]);
	}
}

