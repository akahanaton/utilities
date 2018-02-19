#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <allgene.gff-sort> <chr.size> <up_dist> <down_dist>
	allgene.gff-sort: only gene information is needed
};

getopts('', \%opts);
die($usage) if (@ARGV != 4);
my ($gene_gff, $chr_file, $up_d, $down_d) = @ARGV;

my (%chrInfo, %chrSize);
#Chr1    MSU_osa1r6      gene    10218   11435   .       +       .       ID=13101.t00002;Name=expressed%20protein;Alias=LOC_Os01g01019
open(H, $gene_gff) || die("** Fail to open $gene_gff\n");
while(my $line = <H>)
{
	chomp $line;
	my @tmpArr = split(/\t+/,$line);
	$tmpArr[8] =~ m/(Os\d+g\d+)/;
	my $geneID =  $1;
	if( defined $chrInfo{$tmpArr[0]}{$tmpArr[6]}){
		push(@{$chrInfo{$tmpArr[0]}{$tmpArr[6]}{"start"}}, $tmpArr[3]);
		push(@{$chrInfo{$tmpArr[0]}{$tmpArr[6]}{"end"}}, $tmpArr[4]);
		push(@{$chrInfo{$tmpArr[0]}{$tmpArr[6]}{"id"}}, $geneID);
	}else{
		$chrInfo{$tmpArr[0]}{$tmpArr[6]}{"start"} = ();
		$chrInfo{$tmpArr[0]}{$tmpArr[6]}{"end"} = ();
		$chrInfo{$tmpArr[0]}{$tmpArr[6]}{"id"} = ();
		push(@{$chrInfo{$tmpArr[0]}{$tmpArr[6]}{"start"}}, $tmpArr[3]);
		push(@{$chrInfo{$tmpArr[0]}{$tmpArr[6]}{"end"}}, $tmpArr[4]);
		push(@{$chrInfo{$tmpArr[0]}{$tmpArr[6]}{"id"}}, $geneID);
	}
}
close H;

open(C, $chr_file) || die("** Fail to open $chr_file\n");
while(my $line = <C>)
{
    chomp $line;
    my @tmpArr = split(/\s+/,$line);
    $chrSize{$tmpArr[0]} = $tmpArr[1];
}
close C;

#chr01   .       miRNA_primary_transcript        672278  672374  .       -       .       ID=MI0008244_1;accession_number=MI0008244;Name=osa-MIR1859
foreach my $curChr (keys %chrInfo){
    foreach my $curStr (keys %{$chrInfo{$curChr}}){
		my @curStart = @{$chrInfo{$curChr}{$curStr}{"start"}};
		my @curEnd   = @{$chrInfo{$curChr}{$curStr}{"end"}};
		my @curIDs   = @{$chrInfo{$curChr}{$curStr}{"id"}};
		if( $curStr eq "+"){
			my($downEndLen, $upBeginLen) = (0,0);
			printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[0],  $curChr, $curStr, "up", $curStart[0]>$up_d?$curStart[0]-$up_d:1, $curStart[0]-1;
			$downEndLen = $curStart[1]-$curEnd[0];
			printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[0],  $curChr, $curStr, "down", $curEnd[0]+1,  $downEndLen > $down_d? $curEnd[0]+$down_d:$curStart[1]-1;

			for(my $i=1; $i <= $#curStart-1; $i++){
				$upBeginLen = $curStart[$i]-1 - $curEnd[$i-1];
				printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[$i], $curChr, $curStr, "up", $upBeginLen>$up_d?$curStart[$i]-1-$up_d:$curEnd[$i-1], $curStart[$i]-1;
				$downEndLen = $curStart[$i+1]-1 - $curEnd[$i];
				printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[$i], $curChr, $curStr, "down",$curEnd[$i]+1, $downEndLen>$down_d?$curEnd[$i]+$down_d:$curStart[$i+1]-1;
			}

			$upBeginLen = $curStart[$#curStart]-1 - $curEnd[$#curStart-1];
			printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[$#curStart], $curChr, $curStr, "up",$upBeginLen>$up_d?$curStart[$#curStart]-1-$up_d:$curEnd[$#curStart-1],$curStart[$#curStart]-1;
			$downEndLen = $chrSize{$curChr} - $curEnd[$#curStart];
			printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[$#curStart], $curChr, $curStr, "down",$curEnd[$#curStart]+1, $downEndLen>$down_d?$curEnd[$#curStart]+$down_d:$chrSize{$curChr};
		}else{
			my($downEndLen, $upBeginLen) = (0,0);
			printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[0],  $curChr, $curStr, "down", $curStart[0]>$up_d?$curStart[0]-$up_d:1, $curStart[0]-1;
			$downEndLen = $curStart[1]-$curEnd[0];
			printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[0],  $curChr, $curStr, "up", $curEnd[0]+1,  $downEndLen > $down_d? $curEnd[0]+$down_d:$curStart[1]-1;

			for(my $i=1; $i <= $#curStart-1; $i++){
				$upBeginLen = $curStart[$i]-1 - $curEnd[$i-1];
				printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[$i], $curChr, $curStr, "down", $upBeginLen>$up_d?$curStart[$i]-1-$up_d:$curEnd[$i-1], $curStart[$i]-1;
				$downEndLen = $curStart[$i+1]-1 - $curEnd[$i];
				printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[$i], $curChr, $curStr, "up",$curEnd[$i]+1, $downEndLen>$down_d?$curEnd[$i]+$down_d:$curStart[$i+1]-1;
			}

			$upBeginLen = $curStart[$#curStart]-1 - $curEnd[$#curStart-1];
			printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[$#curStart], $curChr, $curStr, "down",$upBeginLen>$up_d?$curStart[$#curStart]-1-$up_d:$curEnd[$#curStart-1],$curStart[$#curStart]-1;
			$downEndLen = $chrSize{$curChr} - $curEnd[$#curStart];
			printf "%s\t%s\t%s\t%s\t%d\t%d\n", $curIDs[$#curStart], $curChr, $curStr, "up",$curEnd[$#curStart]+1, $downEndLen>$down_d?$curEnd[$#curStart]+$down_d:$chrSize{$curChr};
		}
	}
}
