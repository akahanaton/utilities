#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/basename dirname/;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <lib1.name> <lib2.name>
};

getopts('', \%opts);
die($usage) if (@ARGV != 2);
my ($libName1,$libName2) = @ARGV;

my(%tagCounter);

open(LIB1, $libName1) || die("** Fail to open '$libName1'.\n");
while(my $line = <LIB1>)
{
	chomp $line;
	my @tmpArr = split(/\s+/,$line);
	$tagCounter{uc($tmpArr[0])} = $tmpArr[1];
}
close LIB1;

print "Tags present in both file\n";

open(LIB2, $libName2) || die("** Fail to open '$libName2'.\n");
my %tagCounter_lib2_ongly;
while(my $line = <LIB2>)
{
	chomp $line;
	my @tmpArr = split(/\s+/,$line);
	if (defined($tagCounter{uc($tmpArr[0])}) ){
		printf "%s\t%d\t\%d\n", $tmpArr[0], $tagCounter{uc($tmpArr[0])}, $tmpArr[1];
		delete $tagCounter{uc($tmpArr[0])};
	}else{
		$tagCounter_lib2_ongly{uc($tmpArr[0])} = $tmpArr[1];
	}
}
close LIB2;

print "\n\nTags present only in $libName1\n";
for my $key (keys %tagCounter){
	printf "%s\t%d\t%d\n", $key, $tagCounter{$key}, 0;
}

print "\n\nTags present only in $libName2\n";
for my $key (keys %tagCounter_lib2_ongly){
	printf "%s\t%d\t%d\n", $key, 0,$tagCounter_lib2_ongly{$key};
}
#--------------------------------------------------
# print Dumper(%mirnaLen);
#-------------------------------------------------- 
