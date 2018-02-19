#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <fa.file> 
};

sub complement {
	$_[0] =~ y/CGAT/GCTA/;
	return $_[0];
}

getopts('', \%opts);
die($usage) if (@ARGV != 1);
my ($fas) = @ARGV;

open(G, $gmap_result) || die("** Fail to open $gmap_result.\n");
while(my $line = <G>)
{
	print $line if $line =~ m/^#/;
	if ( $line =~ m/131(\d+?)\..*Chr(\d+)/){
		print $line if $1 eq $2;
	}elsif ( $line =~ m/Chr(\d+)\s+.*131(\d+?)\./){
		print $line if $1 eq $2;
	}
}
close G;

my @sequences;
my @ids;
my $seqin=Bio::SeqIO->new(-file=>$fas,-format=>"fasta");
while(my $seq=$seqin->next_seq()){
	push @sequences, $seq->seq();
	push @ids, $seq->id();
}

my $seqOut = Bio::SeqIO->new(-file => ">$dir/$id\_merged.fa", -format => 'fasta');
    my $out_seq_str = join("N" x 60, values %sequences)
    my $out_seq = Bio::Seq->new(
    -seq => $out_seq_str,
    -display_id => $id,
    -alphabet => "dna" );
$seqOut->write_seq($out_seq);

my $payoff = { match      => 4,  # $match,
			mismatch   => -3, # $mismatch,
			gap_open   => -2, # $gap_open,
			gap_extend => -1 # $gap_extend
		};


for(my $i=1; $i <= $#sequences; $i++){
	my $nw = new Align::NW $sequences[0], $sequences[$i], $payoff;
	$nw->score;
	$nw->align;
	print "$ids[0]\t$ids[$i]\n";
	#--------------------------------------------------
	# print "$sequences[0]\t$sequences[$i]\n";
	#-------------------------------------------------- 
	my $score = $nw->get_score;
	my $align = $nw->get_align;
	$nw->print_align;
	print "$score\n";
	print "\n";

	#--------------------------------------------------
	# my $tmpSeq = reverse $sequences[$i];
	# $tmpSeq = complement $tmpSeq;
	#-------------------------------------------------- 

	#--------------------------------------------------
	# $nw = new Align::NW $sequences[0], $tmpSeq, $payoff;
	# $nw->score;
	# $nw->align;
	# my $score = $nw->get_score;
	# my $align = $nw->get_align;
	# $nw->print_align;
	# print "$score\n";
	# print "$align\n";
	# print "\n";
	#-------------------------------------------------- 

}
#--------------------------------------------------
# $nw->dump_score;
#-------------------------------------------------- 

my %work = map { $_ => 1 } 0 .. $#TE_list;
my @works = keys %work;
my %pids;
while (%work) {
    while (@works and keys %pids < $processes) {
	  my $cur_work = shift @works;
	  die "could not fork" unless defined(my $pid = fork);
	  # parent
	  if ($pid) {
		$pids{$pid} = 1;
		next;
	  }
	  system ("maq pileup -s -P $aligned\/$TE_list[$cur_work]\_LEFT\.bfa $aligned\/\..\/mapping\/out.map \> $aligned\/$TE_list[$cur_work]\_LEFT.pileup");
	  exit $cur_work;
    }
    my $pid = waitpid -1, WNOHANG;
    if ($pid > 0) {
	  delete $pids{$pid};
	  my $rc = $? >> 8; #get the exit status
	  delete $work{$rc};
    }
    select undef, undef, undef, .10;
}

my %seen = ();
my @unique = grep { ! $seen{ $_ }++ } @array;

my @sorted_exps= sort by_number @exps;
sub by_number {
    if ($a < $b){ -1 } elsif ($a > $b) { 1 } else { 0 }
}
