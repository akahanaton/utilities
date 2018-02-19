#! /usr/bin/perl -w
#===============================================================================
#         FILE:  dnaDistanceK2P.pl
#
#  DESCRIPTION:  Calculate DNA Distance using Kimura 2-Parameters Model
#        NOTES:  ---
#       AUTHOR:  Lv Yang
#      VERSION:  1.0
#      CREATED:  01/22/2010
#     REVISION:  11/06/2011
#===============================================================================
my $version = 1.00;

use strict;
use Getopt::Long;

#necessary arguments
my %opts;
GetOptions(\%opts, "i=s", "o=s", "n1=s", "n2=s", "h");
&usage if (!(defined $opts{i} and defined $opts{o} and defined $opts{n1} and defined $opts{n2}) or defined $opts{h});

my $input   =   $opts{i};
my $output  =   $opts{o};
my $tax_id1 =   $opts{n1};
my $tax_id2 =   $opts{n2};

#create file handle
open my $in, "<", "$input" or die "Cannot open file $input: $!.\n";
open my $out, ">", "$output" or die "Cannot create file $output: $!.\n";

my %aln_seq = ( );

print STDERR "Generating hash...\n";
while(<$in>) {
    chomp;
    my ($trid, $taxon, $seq) = split;
    $aln_seq{$trid}{$taxon} = $seq;
}

print STDERR "Calculating the divergence matrix...\n";

for my $trid (sort keys %aln_seq) {
    print STDERR "Processing $trid...\n";
    print $out "$trid\t";
    next unless exists $aln_seq{$trid}{$tax_id1} and $aln_seq{$trid}{$tax_id2};
    my($div, $p, $q) = &div_k2p($aln_seq{$trid}{$tax_id1}, $aln_seq{$trid}{$tax_id2});
    if(defined($div) && $div < 5){
        print $out "$div\t$p\t$q";
    }
    else {
        print $out "nan\tnan\tnan";
    }
    print $out "\n";
}

#close file handle
close $in;
close $out;
	
#========================================SUB==============================================
sub div_k2p {
    my($seq1, $seq2) = @_;
    my ($transition, $transversion, $n, $k, $p, $q) = (0, 0, 0, 0, 0, 0);
    for my $i (0..length($seq1)-1) {
        my $base1 = substr($seq1, $i, 1);
        my $base2 = substr($seq2, $i, 1);
        next if($base1 eq "-" or $base2 eq "-" or $base1 eq "N" or $base2 eq "N");
        $n++;
        $base1 =~ tr/atcgu/ATCGT/;
        $base1 =~ tr/atcgu/ATCGT/;
        next if($base1 eq $base2); 
        if($base1 eq "A" and $base2 eq "G") {
            $transition++;
        } 
        elsif($base1 eq "T" and $base2 eq "C") {
            $transition++;
        } 
        elsif($base1 eq "C" and $base2 eq "T") {
            $transition++;
        } 
        elsif($base1 eq "G" and $base2 eq "A") {
            $transition++;
        }
        else {
            $transversion++;
        }
    }

    #1) Two sequences are not available
    if(!$n) {
        return (0, 0, 0);
    }

    #2) K2P distance
    $p = $transition / $n;
    $q = $transversion / $n;
    if(1-2*$p-$q > 0 and 1-2*$q >0 ) {
        $k = (-1/2)*log((1-2*$p-$q)*sqrt(1-2*$q));
        return (abs($k), $p, $q);
    } 

    #3) Two sequences are far diverged and cannot calculate K2P distance
	else {
        return (5, -1, -1);
    }
}
#========================================SUB==============================================
sub usage{
    print <<"USAGE";
Version $version
Usage:
    $0 -i <input file> -o <output file>
options:
    -i  input file (sequence alighment tab-delimated:Name TaxID Aligned_sequence )
    -o  output file
    -n1 species id 1 (e.g. dm3)
    -n2 species id 2 (e.g. droSim1)
    -h  help
USAGE
    exit(1);
}
#========================================END==============================================
