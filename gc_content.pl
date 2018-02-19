#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

=head1 Description

    This script will first calculate, for each nucleotide position, the %GC of the surrounding nucleotid sequence according to the window size.
    For the first nucleotid, half of the window size will be taken from the end of the sequence, because this script has been developped for cicular genoms.
    Then, i will print a value of % each "step" nucleotides
	
=head1 Usage

    perl gc_content.pl --fasta genome.fasta [--window 1000] [--step 100] [--log] [--help] 
    
    
    arguments details :
    --fasta   a one sequence fasta file, it can be multiline. Program abort if multiple sequences are found one sequence only window(even number ONLY) step(even number ONLY)
    --window  optional. Window size, even number ONLY. Sets the number of nucleotides used to calculate the %GC value of each position. Default 1000
    --step  optional. step size. The output will contain a sliding GC% value every "step" nucleotides. The numbers of values you get is therefore (length genome)/(step). Default 100
    --log   optional. for debugging purposes only
    --help   optional. Shows this help
   
=head1 Output files

    genome.fasta.gc_content
    
    column 1: "chr1" for convenience when the program is used to create an input for Circos (http://www.circos.ca/) 
    column 2: start position of interval
    column 3: end position of interval 
    column 4: GC% of the interval 
    
    example:
    chr1 1 100 0.499
    chr1 101 200 0.5
    chr1 201 300 0.5
    chr1 301 400 0.5
    chr1 401 500 0.5
    chr1 501 600 0.5
    chr1 601 700 0.499
    chr1 701 800 0.5

    genome.fasta.gc_deviation

    column 1: "chr1" for convenience when the program is used to create an input for Circos (http://www.circos.ca/) 
    column 2: start position of interval
    column 3: end position of interval 
    column 4: GC deviation of the interval 
    
    example:
    chr1 1 100 0.162
    chr1 101 200 0.16
    chr1 201 300 0.16
    chr1 301 400 0.16
    chr1 401 500 0.16
    chr1 501 600 0.16
    chr1 601 700 0.162

=cut


my ($fasta,$help,$log);
GetOptions(
	"fasta=s" => \$fasta,
	"window=s" => \(my $window = 1000),
	"step=s" => \(my $step = 100),
	"log" => \$log,
	"help" => \$help
  );

die `pod2text $0` unless $fasta;
die `pod2text $0` if $help;

$fasta =~ s/\r?\n//g;
$window =~ s/\r?\n//g;
$step =~ s/\r?\n//g;


if ($window % 2 == 0)
{print "Okay, $window is an even number.\n\n"}
else
{print "\n!!! Warning : window must be an even number... !!!\n\n"; die `pod2text $0`}


my $output;
$output = $fasta . ".gc_content";
#if ($fasta =~ /\//){ $output = (split(/\//,$fasta))[-1] . "gc_content" }else{ $output = $fasta . "gc_content"}

print "output file is $output\n";

open (OUT, ">", $output ) or die "can't open $!";
open (LOG, ">", $output . ".log") or die "can't open $!";
open (IN, "<", $fasta) or die "can't open  $!";

my (%h_pcent, $size_fasta, $sequence, $not_sequence);
my %h_gcdev; #GC_deviation
$not_sequence = 0;

while (my $line = <IN>) {
	$line =~ s/\r?\n//g;
	$line =~ s/\s//g;
	$line !~ m/^>/ ? $sequence .= $line : $not_sequence++;
}

#### Running few tests ... ####
if ($not_sequence > 1) {print "Your fasta file contains more than one sequence (i found the > symbol $not_sequence times, i'll die now ...\n"; die}
if ($sequence =~ m/[^ATGC]/)
{print "Your fasta sequence contains other caracters than [ATGC], i'll keep working normally, it's just so you know ...\n ";}

$size_fasta = length($sequence);
print "fasta file is $size_fasta nucleotid long\n\n";

#### Calculating mean GC content ####
my $count1 = 0;
while ($sequence =~ /[GC]/g) { $count1++ }
my $pcent =  $count1 / $size_fasta * 100 ;
print "It's computationaly-free to display it, so while i'm at it ...\nthe mean GC content of the sequence is $pcent % \n";
my $line_mod = (substr $sequence, -$window/2 ) . $sequence . (substr $sequence, 0, ($window / 2));
print LOG $line_mod . "\n" if $log;

my $start = 0;

until ($start > (length($sequence)))
{
	my $count = 0;
	my $countG = 0; my $countC = 0;  #GC_deviation
	my $tmp = substr $line_mod, $start, $window;
	my $position = $start +1 ;
	while ($tmp =~ /[GC]/g) {$count++}
	while ($tmp =~ /[C]/g) {$countC++} #GC_deviation
	while ($tmp =~ /[G]/g) {$countG++} #GC_deviation
	my $gcdev = ($countG-$countC)/($countG+$countC) ; #GC_deviation
	#print LOG $count . " " . $tmp . "\n";
	my $pcent =  $count / $window ;
	$h_pcent{$position} = $pcent;
	$h_gcdev{$position} = $gcdev; #GC_deviation
	print LOG $position . " " . $tmp . " $pcent\n" if $log;
	$start ++;  # $start += $window; if i choose to calculate the mean on the windw instead of making it for each nucleotid
}

my $k = 0;
my $acc = 0;



my $last_end;
foreach my $key2 ( sort {$a<=>$b} keys %h_pcent)
{
	#if ($key2 % $step == 0 || $key2  == ($size_fasta + 1)) {
	if ($key2 % $step == 0 || $key2  == ($size_fasta + 1)) {
		my $start_window2 = $key2 - ($step - 1);
		my $end_window2 = $key2; if($end_window2 >= $size_fasta){$start_window2 = $last_end + 1; $end_window2 = $size_fasta }
		print OUT "chr1 " . $start_window2 . " " . $end_window2 . " " . $h_pcent{$key2} . "\n";
		$last_end = $end_window2;
	}
}
close OUT;

open (OUT3, ">", $fasta . ".GC_deviation") or die "can't open $!";
my $last_end2;
foreach my $key3 ( sort {$a<=>$b} keys %h_gcdev)
{
	#if ($key2 % $step == 0 || $key2  == ($size_fasta + 1)) {
	if ($key3 % $step == 0 || $key3  == ($size_fasta + 1)) {
		my $start_window2 = $key3 - ($step - 1);
		my $end_window2 = $key3; if($end_window2 >= $size_fasta){$start_window2 = $last_end2 + 1; $end_window2 = $size_fasta }
		print OUT3 "chr1 " . $start_window2 . " " . $end_window2 . " " . $h_gcdev{$key3} . "\n";
		$last_end2 = $end_window2;
	}
}
close OUT3;

close IN;
close OUT;
close LOG;

unlink $output . ".log" if !$log;
