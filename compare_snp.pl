#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/basename dirname/;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <snp2>  <snp3>
};


getopts('', \%opts);
die($usage) if (@ARGV != 2);
my ($file2,$file3) = @ARGV;

my(%snp2, %snp3);

open(F2, $file2) || die("** Fail to open $file2.\n");
while(my $line = <F2>)
{
	chomp $line;
	my @arr = split(/\s+/, $line);
	if(! defined($snp2{$arr[0]})){
		$snp2{$arr[0]}{str} = $arr[2];
		$snp2{$arr[0]}{reffile} = $arr[1];
		$snp2{$arr[0]}{allel} = $arr[3];
	}
}
close F2;

open(F3, $file3) || die("** Fail to open $file3.\n");
while(my $line = <F3>)
{
	chomp $line;
	my @arr = split(/\s+/, $line);
	if(! defined($snp3{$arr[0]})){
		$snp3{$arr[0]}{str} = $arr[3];
		$snp3{$arr[0]}{allel} = $arr[4];
	}
}
close F3;


my (@prob_snps_miss2, @prob_snps_miss3, @prob_snps_err);
foreach my $snp ( keys %snp2 ){
	if( defined($snp2{$snp}{str}) ) {
		if (defined($snp3{$snp}{str}) ){
			if ( $snp2{$snp}{str} eq '+' ){
				my @ale2 = split(/\//,$snp2{$snp}{allel});
				#--------------------------------------------------
				# if ($snp2{$snp}{allel} eq 'A/T' || $snp2{$snp}{allel} eq 'T/A' || $snp2{$snp}{allel} eq 'G/C' || $snp2{$snp}{allel} eq 'C/G'){
				# 	print "$snp\t\t$snp2{$snp}{reffile}\t$snp2{$snp}{str}\t$snp2{$snp}{allel}\t$file3\t$snp3{$snp}{str}\t$snp3{$snp}{allel}\n";
				# }
				#-------------------------------------------------- 
				foreach my $ale ( @ale2 ){
					if (! $snp3{$snp}{allel} =~ /$ale/){
						push(@prob_snps_err, $snp);
					}
				}
			}else{
				my @ale2 = split(/\//,$snp2{$snp}{allel});
				if ($snp2{$snp}{allel} eq 'A/T' || $snp2{$snp}{allel} eq 'T/A' || $snp2{$snp}{allel} eq 'G/C' || $snp2{$snp}{allel} eq 'C/G'){
					print "$snp\t\t$snp2{$snp}{reffile}\t$snp2{$snp}{str}\t$snp2{$snp}{allel}\t$file3\t$snp3{$snp}{str}\t$snp3{$snp}{allel}\n";
				}
				#--------------------------------------------------
				# print "$snp\t\t$snp2{$snp}{reffile}\t$snp2{$snp}{str}\t$snp2{$snp}{allel}\t$file3\t$snp3{$snp}{str}\t$snp3{$snp}{allel}\n";
				#-------------------------------------------------- 
				foreach my $ale ( @ale2 ){
					$ale =~ tr/ATCG/TAGC/ ;
					if (! $snp3{$snp}{allel} =~ /$ale/){
						push(@prob_snps_err, $snp);
					}
				}
			}
		}else{
			push(@prob_snps_miss3, $snp);
		}
	}else{
		push(@prob_snps_miss2, $snp);
	}
}
#--------------------------------------------------
# print Dumper(%mirnaLen);
#-------------------------------------------------- 

foreach my $snp( @prob_snps_miss2){
	print "none in $file2\t$snp\n";
}

foreach my $snp( @prob_snps_miss3){
	print "only in $file2\t$snp\t$snp2{$snp}{str}\t$snp2{$snp}{allel}\n";
}

foreach my $snp( @prob_snps_err){
	print "error\t$snp\t$snp2{$snp}{reffile}\t$snp2{$snp}{str}\t$snp2{$snp}{allel}\t$file3\t$snp3{$snp}{str}\t$snp3{$snp}{allel}\n";
}

