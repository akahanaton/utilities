#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <dat.file> 
};

getopts('', \%opts);
die($usage) if (@ARGV != 1);
my ($dat_file) = @ARGV;

open(D, $dat_file) || die("** Fail to open $dat_file.\n");
my $fileID = undef;
while(my $line = <D>)
{
	if($line =~ m/^FILE=/)
    {
        if ( $line =~ m/\\(.*)\./){
            $fileID = $1;
            last;
        }
    }
}
print $fileID;
close D;

#--------------------------------------------------
# my %seen = ();
# my @unique = grep { ! $seen{ $_ }++ } @array;
#-------------------------------------------------- 

#--------------------------------------------------
# my @sorted_exps= sort by_number @exps;
# sub by_number {
#     if ($a < $b){ -1 } elsif ($a > $b) { 1 } else { 0 }
# }
#-------------------------------------------------- 
