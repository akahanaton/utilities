#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;

my %opts;
my $usage = qq{
    Usage: $0 <Filename> <LIST.tag> <REVERSE.tag> <ExpectValue or Score>
    set  LIST.tag = 1 to treat  Filename as list
};

sub complement {
	$_[0] =~ y/CGAT/GCTA/;
	return $_[0];
}

getopts('', \%opts);
die($usage) if (@ARGV != 4);
my ( $file_list, $list_tag, $reverse_tag, $value_tag) = @ARGV;

my @peptidesFiles;
if($list_tag){
    open(F, "$file_list") || die("** Fail to open $file_list.\n");
    while(my $file = <F>)
    {
        chomp $file;

        # Use a regular expression to ignore files beginning with a period
        next if ($file !~ m/.peptides$/);

        push(@peptidesFiles, $file);
    }
}
else{
    push(@peptidesFiles, $file_list);
}
my %fdr_data;
foreach my $file(@peptidesFiles){
    open(M, $file) || die("** Fail to open $file.\n");
    my $line = <M>;
    my $ValueIndex = undef;  
    my $SequenceIndex = undef;  
    my $ReferenceIndex = undef;
    my $ScanIndex = undef;
    my @arr = split(/\t/, $line);
    for(my $i=0; $i <= $#arr; $i++){
       $ValueIndex = $i if($arr[$i] eq $value_tag); 
       $SequenceIndex = $i if($arr[$i] eq "Sequence"); 
       $ReferenceIndex = $i if($arr[$i] eq "Reference"); 
       $ScanIndex = $i if($arr[$i] eq "FileScan"); 
    }
    while(my $line = <M>)
    {
        last if $line =~ m/^\s+$/;
        my @arr = split(/\t/, $line);
        my (undef, $scanNum) = split(',', $arr[$ScanIndex]);
        if(defined($fdr_data{$arr[$ValueIndex]})){
            push (@{$fdr_data{$arr[$ValueIndex]}{Sequence}}, $arr[$SequenceIndex]);
            push (@{$fdr_data{$arr[$ValueIndex]}{Reference}}, $arr[$ReferenceIndex]);
            push (@{$fdr_data{$arr[$ValueIndex]}{SCANID}}, $scanNum);
            if($arr[$ReferenceIndex] =~ m/$reverse_tag/)
            {
                push (@{$fdr_data{$arr[$ValueIndex]}{Reverse}}, 1);
            }else{
                push (@{$fdr_data{$arr[$ValueIndex]}{Reverse}}, 0);
            }
        }else{
            $fdr_data{$arr[$ValueIndex]}{Sequence} = ();
            $fdr_data{$arr[$ValueIndex]}{Reference} = ();
            $fdr_data{$arr[$ValueIndex]}{Reverse} = ();
            $fdr_data{$arr[$ValueIndex]}{SCANID} = ();
            push (@{$fdr_data{$arr[$ValueIndex]}{Sequence}}, $arr[$SequenceIndex]);
            push (@{$fdr_data{$arr[$ValueIndex]}{Reference}}, $arr[$ReferenceIndex]);
            push (@{$fdr_data{$arr[$ValueIndex]}{SCANID}}, $scanNum);
            if($arr[$ReferenceIndex] =~ m/$reverse_tag/)
            {
                push (@{$fdr_data{$arr[$ValueIndex]}{Reverse}}, 1);
            }else{
                push (@{$fdr_data{$arr[$ValueIndex]}{Reverse}}, 0);
            }
        }
    }
    close M;

}
close F;


my @Scores;
if($value_tag  eq  "Score")
{
    @Scores =  sort {$b <=> $a} (keys %fdr_data);
}else{
    @Scores =  sort {$a <=> $b} (keys %fdr_data);
}

my $targetNum = 0;
my $reverseNum = 0;
my $totalNum = 0;
for(my $i=0; $i<=$#Scores; $i++)
{
    my @reverseArr = @{$fdr_data{$Scores[$i]}{Reverse}};
    # method 1
    foreach my $reversed (@reverseArr){
        if($reversed){
            $reverseNum++;
        }
        $totalNum++;
        my $TPFP = $totalNum - $reverseNum;
        if($TPFP > 0){
            $fdr_data{$Scores[$i]}{FDR} = $reverseNum / $TPFP;
        } else{
            $fdr_data{$Scores[$i]}{FDR} = 0;
        }
    }

    # method 2
    #--------------------------------------------------
    # foreach my $reversed (@reverseArr){
    #     if($reversed){
    #         $reverseNum++;
    #     }else{
    #         $targetNum++;
    #     }
    # }
    # my $totalNum = $targetNum;
    # if($totalNum>0){
    #     $fdr_data{$Scores[$i]}{FDR} = $reverseNum / $totalNum;
    # }else{
    #     $fdr_data{$Scores[$i]}{FDR} = 0;
    # }
    #-------------------------------------------------- 
}
#--------------------------------------------------
# print Dumper(%fdr_data);
#-------------------------------------------------- 

for(my $i=0; $i<=$#Scores; $i++){
    #--------------------------------------------------
    # if($fdr_data{$Scores[$i]}{FDR} <= 0.01)
    # {
    #-------------------------------------------------- 
        my @sequenceArr = @{$fdr_data{$Scores[$i]}{Sequence}}; 
        my @reverseArr = @{$fdr_data{$Scores[$i]}{Reverse}}; 
        my @referenceArr = @{$fdr_data{$Scores[$i]}{Reference}}; 
        my @scanArr = @{$fdr_data{$Scores[$i]}{SCANID}}; 
        for(my $j=0; $j<=$#sequenceArr; $j++){
            if(!$reverseArr[$j]){
                print "$Scores[$i]\t$fdr_data{$Scores[$i]}{FDR}\t$scanArr[$j]\t$sequenceArr[$j]\t$referenceArr[$j]\n";
            }
        }
    #--------------------------------------------------
    # }
    #-------------------------------------------------- 
}
exit 0;
