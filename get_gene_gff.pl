#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <gene.id.file> <gff.file>
};

getopts('', \%opts);
die($usage) if (@ARGV != 2);
my ($id_file, $gff_file) = @ARGV;

my %id;
open(ID, $id_file) || die("** Fail to open '$id_file'.\n");
while ( my $line = <ID>){
	chomp $line;
	my @arr = split(/\s+/, $line);
	$id{$arr[0]} = $line;
}
close ID;

open(GFF, $gff_file) || die("** Fail to open '$gff_file'.\n");
while( my $line = <GFF>)
{
    my @arr = split(/\s+/, $line);
    print $line if ( $arr[2] eq "gene");
    if($arr[8] =~ m/ID=(.*?)[;:]/ ){
	#--------------------------------------------------
	#   print $1,"\n";
	#-------------------------------------------------- 
	  if(defined($id{$1})){
		print $line;
	  }
    }
}
close GFF;
