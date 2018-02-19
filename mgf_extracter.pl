#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;

my %opts;
my $usage = qq{
    Usage: $0 <mgf.file> <ID or ID file>
    multiple ID supported;
};

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
if(-e $ARGV[1]){
    my $scanID_file = $ARGV[1];
    open(F, $scanID_file) || die("** Fail to open $scanID_file.\n");
    while(my $line = <F>)
    {
        chomp $line;
        $scanMatched{$line}=1;
    }
    close F;
}else{
    for(my $i=1; $i<=$#ARGV; $i++){
        $scanMatched{$ARGV[$i]}=1;
    }
}
#--------------------------------------------------
# my @matched =  keys %scanMatched;
# print $#matched,"\t", $matched[0], "\n";
#-------------------------------------------------- 

#--------------------------------------------------
# foreach my $line (@headLines)
#     {print $line;}
#-------------------------------------------------- 

my @scans =  sort {$a <=> $b} (keys %scanInfo);
foreach my $scan (@scans){
    if(defined($scanMatched{$scan})){
        print "###MSMS: $scan/$scan\n";
        foreach my $line (@{$scanInfo{$scan}}){
            print $line;
        }
    }
}
#--------------------------------------------------
# print $#scans,"\n";
#-------------------------------------------------- 
