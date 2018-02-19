#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <tag.id.file> <id_contain_file> <pos_beg> <string_len>
	pos_beg count begin with 1
};

getopts('', \%opts);
die($usage) if (@ARGV != 4);
my ($tag_id, $id2id_file,$beg, $len) = @ARGV;

my (@tags, @lib1_info);

open(COMP, $tag_id) || die("** Fail to open $tag_id.\n");
while(my $line = <COMP>)
{
	chomp $line;
	push(@tags,$line);
}
close COMP;

my %tagInfo;
open(F1, $id2id_file) || die("** Fail to open '$id2id_file'.\n");
while( my $line = <F1>)
{
	if($line =~ />/){
		my $id = substr($line,$beg-1,$len);
		$tagInfo{$id} = $line;
	}
}
close F1;

for my $eachkey (@tags)
{
	if (defined($tagInfo{$eachkey})){
		printf "%s", $tagInfo{$eachkey};
	}else{
		printf "nomatch\t %s\n", $eachkey;
	}
}
