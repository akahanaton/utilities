#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <mgf.file> <fdr.file1> <fdr.file2> <fdr.file3> .......
    multiple fdr files supported
};

sub complement {
	$_[0] =~ y/CGAT/GCTA/;
	return $_[0];
}

getopts('', \%opts);
die($usage) if (@ARGV < 2);

my $mgf_file = $ARGV[0];

my @headLines;
my %scanInfo;
my $scan = undef;
my $begin = 0;
open(M, $mgf_file) || die("** Fail to open $mgf_file.\n");
while(my $line = <M>)
{
    if($line =~ m/^###MSMS:/)
    {
        $begin = 1;
        chomp $line;
        if($line =~ m/\//){
            ($scan) = ($line =~ m/\s(\d+)\//);
        }else{
            ($scan) = ($line =~ m/\s(\d+)/);
        }
        $scanInfo{$scan} = ();
        while( $line = <M>){
            push(@{$scanInfo{$scan}}, $line);
            if( $line =~ /^END IONS/){
                $line = <M>;
                push(@{$scanInfo{$scan}}, $line);
                last;
            }
        }
    }elsif( !  $begin ){
        push(@headLines, $line);
    }
}
close M;

#--------------------------------------------------
# 5.08E-025       0       16246   R.HQGVMVGMGQKDSYVGDEAQSKR.G
#-------------------------------------------------- 
my %scanMatched;
for(my $i = 1; $i <=$#ARGV; $i++){
    my $fdr_file = $ARGV[$i];
    open(F, $fdr_file) || die("** Fail to open $fdr_file.\n");
    while(my $line = <F>)
    {
        my @arr = split(/\t/, $line);
        if($arr[1] <= 0.01){
            $scanMatched{$arr[2]} = 1;
        }
    }
    close F;
}
#--------------------------------------------------
# my @matched =  keys %scanMatched;
# print $#matched,"\t", $matched[0], "\n";
#-------------------------------------------------- 

foreach my $line (@headLines)
    {print $line;}

my @scans =  sort {$a <=> $b} (keys %scanInfo);
foreach my $scan (@scans){
    if(!defined($scanMatched{$scan})){
        print "###MSMS: $scan/$scan\n";
        foreach my $line (@{$scanInfo{$scan}}){
            print $line;
        }
    }
}
#--------------------------------------------------
# print $#scans,"\n";
#-------------------------------------------------- 
