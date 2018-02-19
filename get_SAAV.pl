#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <snp.fasta> <raw.fasta> <snp.code>
};

getopts('', \%opts);
die($usage) if (@ARGV != 3);
my ($snp_fas, $raw_fas, $snp_code_file) = @ARGV;

my %refSeq;
my $seqin=Bio::SeqIO->new(-file=>$snp_fas,-format=>"fasta");
while(my $seq=$seqin->next_seq()){
    $refSeq{$seq->id} = $seq->seq();
    #--------------------------------------------------
    # $refSeq{$seq->id} =~ s/I/L/g;
    #-------------------------------------------------- 
}

my %mutationInfo;
my $seqin2=Bio::SeqIO->new(-file=>$raw_fas,-format=>"fasta");
while(my $seq=$seqin2->next_seq()){
    if($seq->description =~ m/Variant=/){
        my ($varintInfo) = ($seq->description =~ m/Variant=(.*)\s+\\Process?/);
        my ($seqID) = ($seq->id =~ m/NX_(.*)/);
        foreach my $var( $varintInfo =~ /\|(\d+\|\w)\)/g){
            #--------------------------------------------------
            # print $var,"\n";
            #-------------------------------------------------- 
            my @arr = split('\|', $var);
            my $acid = substr($seq->seq, $arr[0]-1, 1);
            $mutationInfo{$seqID}{$arr[0]} = $acid."->".$arr[1];
        }
    }
}

#--------------------------------------------------
# print Dumper(%refSeq);
#-------------------------------------------------- 
my %snpInfo;
open(S, $snp_code_file) || die("** Fail to open $snp_code_file.\n");
while( my $line = <S>) {
    if( $line=~ /^Pep/){
        chomp $line;
        my @arr = split('\t', $line);
        next if $arr[2] < 2;
        next if($arr[1] =~ m/:.*:/);
        #--------------------------------------------------
        # print Dumper(@arr);
        #-------------------------------------------------- 
        $snpInfo{$arr[8]}{$arr[1]}{code} = $arr[2];
        $snpInfo{$arr[8]}{$arr[1]}{count} = $arr[7];
        $snpInfo{$arr[8]}{$arr[1]}{Skewness} = $arr[3];
        $snpInfo{$arr[8]}{$arr[1]}{Kurtosis} = $arr[4];
        $snpInfo{$arr[8]}{$arr[1]}{weighted} = $arr[5];
        $snpInfo{$arr[8]}{$arr[1]}{seq} = $arr[8];
        $snpInfo{$arr[8]}{$arr[1]}{mod} = $arr[9];
        $snpInfo{$arr[8]}{$arr[1]}{line} = $line;
    }
}
close S;

print "Protein ID\tVariation Site at Protein\tVariation Type\tPeptide Identified\tVariation Site at Peptide\tSpectrum Count\tSkewness\tKurtosis\tWeighted Mean\tPeptide Evidence Code\n"; 
foreach my $snpSeq (keys %snpInfo){
    my @snpIDs = keys %{$snpInfo{$snpSeq}};
    my $beginPos= undef;
    if(defined($refSeq{$snpIDs[0]})){
        my ($rawSeqID, $snpRawPos, $snpPepPos) = split(/:|\./, $snpIDs[0]); 
        $beginPos  =  index($refSeq{$snpIDs[0]}, $snpInfo{$snpSeq}{$snpIDs[0]}{seq});
        my $snpPepPosNew = $snpPepPos - $beginPos;
        my $refSeqLen = length $refSeq{$snpIDs[0]};
        my $pepLen = length $snpInfo{$snpSeq}{$snpIDs[0]}{seq};
        if($snpPepPosNew < $pepLen && $beginPos >= 0){
            print "$rawSeqID\t$snpRawPos\t$mutationInfo{$rawSeqID}{$snpRawPos}\t$snpInfo{$snpSeq}{$snpIDs[0]}{seq}\t$snpPepPosNew\t$snpInfo{$snpSeq}{$snpIDs[0]}{count}\t$snpInfo{$snpSeq}{$snpIDs[0]}{Skewness}\t$snpInfo{$snpSeq}{$snpIDs[0]}{Kurtosis}\t$snpInfo{$snpSeq}{$snpIDs[0]}{weighted}\t$snpInfo{$snpSeq}{$snpIDs[0]}{code}\n";
        }
        if($beginPos<0)
        {
            print STDERR  $snpInfo{$snpSeq}{$snpIDs[0]}{line}, "\n";
        }
    }
}
