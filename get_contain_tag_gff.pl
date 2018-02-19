#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <gff.file> <id_file> <Y/N>
	field count begin with 1, 
	Y: print tag-contained line
	N: print non-tag line
};

getopts('', \%opts);
die($usage) if (@ARGV != 3);
my ($tab_file, $id_file, $YN) = @ARGV;

my (@tags, @lib1_info);

my %tagInfo;
open(F1, $id_file) || die("** Fail to open '$id_file'.\n");
while( my $line = <F1>)
{
	chomp $line;
	$tagInfo{$line} = 1;
}
close F1;

open(COMP, $tab_file) || die("** Fail to open $tab_file.\n");
while(my $line = <COMP>)
{
  chomp $line;
  my @arr = split(/\s+/, $line);
  if($arr[2] eq 'gene' || $arr[2] eq 'chromosome'){
	print "$line\n" ;
	next;
  }
  my ($ID) = ($arr[8] =~ m/ID=(.*?)[:;]/);
  if($YN eq 'Y'){
    if (defined($tagInfo{$ID})){
	print "$line\n";
    }
  }else{
    if (!defined($tagInfo{$ID})){
	print "$line\n";
    }
  }
}
close COMP;
