#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::Perl;
use Bio::AlignIO;
use Data::Dumper;
use Getopt::Std;
use List::Util qw[min max];

my %opts;
my $usage = qq{
  Usage: $0 <aln.file>  <aln.format> <begin> <end>
};

getopts('', \%opts);
die($usage) if (@ARGV != 4);

exit(main(@ARGV));

sub main{
    my($lcbFile,$format, $begin, $end) = @_;
    my $in=Bio::AlignIO->new(-file=>$lcbFile,-format=>$format);
    my $aln=$in->next_aln;
    my $slicedAln=$aln->slice($begin,$end);
	$slicedAln->set_displayname_flat();

    my $outFile="$lcbFile.$begin.$end";
    my $out=Bio::AlignIO->new(-format=>$format,-file=>">$outFile");
    $out->write_aln($slicedAln);

    foreach my $seq ( $slicedAln->each_seq() ) {
	#--------------------------------------------------
	#   print $seq->seq(),"\t",$seq->display_id(),"\n";
	#-------------------------------------------------- 
	  printf "%-10s\t%s\n",$seq->display_id(),$seq->seq();
    }
    my $mat_line = $slicedAln->match_line();
    printf "%-10s\t%s\n","Consensus",$mat_line;
}

sub is_col_identical{
    my($matchLine,$col)=@_;
    my $t=substr($matchLine,$col,1);
    if($t eq '*'){
	  return 1;
    }
	  return 0;
}
