#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

my %opts;
my $usage = qq{
Usage: $0 <tag.id.file>  <@ pattern> <id.field> <target.field> <loop_num>
	\@:  for wild match in the filename
	field: 0-based;
};

getopts('', \%opts);
die($usage) if (@ARGV != 5);
my ($tag_id, $pat, $id_field, $tgt_field, $loop_num) = @ARGV;
my (@tags, @lib1_info);



open(COMP, $tag_id) || die("** Fail to open $tag_id.\n");
while(my $line = <COMP>)
{
	chomp $line;
	push(@tags,$line);
}
close COMP;


for (my $i =1; $i <= $loop_num; $i++) 
{
	my %libi;
	my $file1_name = $pat;
	$file1_name =~ s/@/$i/g;
	open(F1, $file1_name) || die("** Fail to open '$file1_name'.\n");
	while( my $line = <F1>)
	{
		next if ( $line =~ /^$/);
		next if ( $line =~ /^#/);
		chomp $line;
		my @arr = split(/\s/, $line);
		$libi{$arr[$id_field]} = $arr[$tgt_field];
	}
	push (@lib1_info, \%libi);
	close F1;
}

for my $eachkey (@tags)
{
	printf "%s\t", $eachkey;
	for (my $i = 0; $i < $loop_num; $i++) 
	{
		if (defined(${lib1_info[$i]}{$eachkey})){
			print ${lib1_info[$i]}{$eachkey},"\t";
		}else{
			print "0\t";
		}
	}
	print "\n";
}
