#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <database.file>  <proteomicstool.fdr>
};

getopts('', \%opts);
die($usage) if (@ARGV != 2);
my ($fas, $pt_file) = @ARGV;

my @replaceCode = ("L", "X"); 

my %refSeq;
my $seqin=Bio::SeqIO->new(-file=>$fas,-format=>"fasta");
while(my $seq=$seqin->next_seq()){
    $refSeq{$seq->id} = $seq->seq();
    $refSeq{$seq->id} =~ s/I/L/g;
}

my %replaceSeq = %refSeq;

open(PT, $pt_file) || die("** Fail to open $pt_file.\n");
my  $line = <PT>;
while( $line = <PT>)
{
    my @arr = split('\t', $line);
    last if $#arr < 3;
    my $curPep = "";
    my $sucess = undef;
    if($arr[1] =~ / ! /){
        my @tmpArr = split(/!/, $arr[1]);
        foreach my $curPep (@tmpArr){
            my $beginIndex = index($curPep,'.');
            my $endIndex = rindex($curPep,'.');
            $curPep =  substr $curPep, $beginIndex+1, $endIndex-$beginIndex-1;
            #--------------------------------------------------
            # print "here1\t",$curPep,"\n";
            #-------------------------------------------------- 
            $sucess = &replace_seq($curPep,$arr[2]);
        }
    } else{
        my $beginIndex = index($arr[1],'.');
        my $endIndex = rindex($arr[1],'.');
        $curPep =  substr $arr[1], $beginIndex+1, $endIndex-$beginIndex-1;
        #--------------------------------------------------
        # print "here2\t",$curPep,"\n";
        #-------------------------------------------------- 
        $sucess = &replace_seq($curPep,$arr[2]);
    }
    if(!defined($sucess))
    {
        print STDERR "pt\t", $curPep,"\t", $line, "\t", $sucess,"\n";
    }

}
close PT;

my %geneCov;
foreach my $key (keys %replaceSeq)
{
    my (undef, $proteinName, undef) = split /\|/, $key;
    my ($geneName, undef) = split /\-/, $proteinName;
    my $covered_num = ($replaceSeq{$key} =~ tr/#//);
    if(defined($geneCov{$geneName})){
        push(@{$geneCov{$geneName}}, $covered_num/length($replaceSeq{$key})*100)
    }else{
        $geneCov{$geneName} = ();
        push(@{$geneCov{$geneName}}, $covered_num/length($replaceSeq{$key})*100)
    }
}
foreach my $key (keys %geneCov)
{
    if($#{$geneCov{$key}} > 1){
        my $sumCov = 0;
        foreach my $cov (@{$geneCov{$key}}){
            $sumCov = $sumCov + $cov;
        }
        print $key,"\t", $sumCov / $#{$geneCov{$key}},"\n";
    }else
    {
        print $key,"\t", $geneCov{$key}[0],"\n";
    }
}

sub replace_seq {
    my($peptide, $proIDs) = @_;
    $peptide =~ s/[\*\^%#&@]//g;
    $peptide =~ s/I/L/g;
    #--------------------------------------------------
    # print "here3\t",$peptide,"\n";
    #-------------------------------------------------- 
    my $sucess = undef;
    my $replaceStr = "#" x length($peptide);
    my @proteinIDs = split /[\/!]/, $proIDs;
    #--------------------------------------------------
    # print Dumper(@proteinIDs);
    #-------------------------------------------------- 
    foreach my $id (@proteinIDs){
        $id =~ s/\s//g;
        if(defined($refSeq{$id})){
            my $beginPos  =  index($refSeq{$id}, $peptide);
            if($beginPos >= 0){
                $replaceSeq{$id} =~ s/\Q$peptide\E/$replaceStr/g;
                $sucess = 1;
            }
            else 
            {
                my $sucess = undef;
                # process "X" replacement
                for(my $i=0; $i<= length $peptide; $i++)
                {
                    my $tmpPep = $peptide;
                    substr($tmpPep, $i, 1, "X");
                    my $beginPos  =  index($refSeq{$id}, $tmpPep);
                    if($beginPos >= 0){
                        $replaceSeq{$id} =~ s/\Q$tmpPep\E/$replaceStr/g;
                        #--------------------------------------------------
                        # print STDERR "good pt\t", $peptide,"\t", $tmpPep, "\t", $id, "\t",  $line;
                        #-------------------------------------------------- 
                        $sucess = 1;
                        last;
                    }
                }
                # process "KQ" replacement
                if(!defined($sucess))
                {
                    for(my $i=0; $i<= length $peptide; $i++)
                    {
                        my $tmpPep = $peptide;
                        if(substr($tmpPep, $i, 1) eq "K")
                        {
                            substr($tmpPep, $i, 1, "Q");
                        }elsif(substr($tmpPep, $i, 1) eq "Q"){
                            substr($tmpPep, $i, 1, "K");
                        }elsif(substr($tmpPep, $i, 1) eq "N"){
                            substr($tmpPep, $i, 1, "B");
                        }elsif(substr($tmpPep, $i, 1) eq "B"){
                            substr($tmpPep, $i, 1, "N");
                        }
                        my $beginPos  =  index($refSeq{$id}, $tmpPep);
                        if($beginPos >= 0){
                            $replaceSeq{$id} =~ s/\Q$tmpPep\E/$replaceStr/g;
                            #--------------------------------------------------
                            # print STDERR "good pt\t", $peptide,"\t", $tmpPep, "\t", $id, "\t",  $line;
                            #-------------------------------------------------- 
                            $sucess = 1;
                            last;
                        }
                    }
                }
                if(!defined($sucess))
                {
                    for(my $i=0; $i<= length $peptide; $i++)
                    {
                        my $tmpPep = $peptide;
                        if(substr($tmpPep, $i, 1) eq "N"){
                            substr($tmpPep, $i, 1, "D");
                        }
                        if(substr($tmpPep, $i, 1) eq "N"){
                            substr($tmpPep, $i, 1, "B");
                        }
                        my $beginPos  =  index($refSeq{$id}, $tmpPep);
                        if($beginPos >= 0){
                            $replaceSeq{$id} =~ s/\Q$tmpPep\E/$replaceStr/g;
                            #--------------------------------------------------
                            # print STDERR "good pt\t", $peptide,"\t", $tmpPep, "\t", $id, "\t",  $line;
                            #-------------------------------------------------- 
                            $sucess = 1;
                            last;
                        }
                    }
                }
                if(!defined($sucess))
                {
                    my $tmpPep = $peptide;
                    for(my $i=0; $i<= length $peptide; $i++)
                    {
                        if(substr($tmpPep, $i, 1) eq "K")
                        {
                            substr($tmpPep, $i, 1, "Q");
                        }elsif(substr($tmpPep, $i, 1) eq "Q"){
                            substr($tmpPep, $i, 1, "K");
                        }elsif(substr($tmpPep, $i, 1) eq "N"){
                            substr($tmpPep, $i, 1, "B");
                        }elsif(substr($tmpPep, $i, 1) eq "B"){
                            substr($tmpPep, $i, 1, "N");
                        }
                        my $beginPos  =  index($refSeq{$id}, $tmpPep);
                        if($beginPos >= 0){
                            $replaceSeq{$id} =~ s/\Q$tmpPep\E/$replaceStr/g;
                            $sucess = 1;
                            last;
                        }
                    }
                }
                if(!defined($sucess))
                {
                    my $tmpPep = $peptide;
                    for(my $i=0; $i<= length $peptide; $i++)
                    {
                        if(substr($tmpPep, $i, 1) eq "N")
                        {
                            substr($tmpPep, $i, 1, "B");
                        }elsif(substr($tmpPep, $i, 1) eq "Q"){
                            substr($tmpPep, $i, 1, "Z");
                        }
                        my $beginPos  =  index($refSeq{$id}, $tmpPep);
                        if($beginPos >= 0){
                            $replaceSeq{$id} =~ s/\Q$tmpPep\E/$replaceStr/g;
                            $sucess = 1;
                            print STDERR "here4\t",$tmpPep, "\t", $peptide,"\t",$sucess,"\n";
                            last;
                        }
                    }
                }
            }
        }
    }
    return($sucess)
}
