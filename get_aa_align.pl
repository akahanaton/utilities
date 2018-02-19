#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::Tools::Run::Alignment::Clustalw;
use File::Basename qw/basename dirname/;
use Data::Dumper;
use MyMod::Bio::Tools::MolGen;
use Bio::Tools::CodonTable;

my %opts;
my $usage = qq{
	Usage: $0 <orth.list> <species_num> <out_aln>
};
getopts('', \%opts);
die($usage) if (@ARGV != 3);
my ($list_file, $species_num, $out_aln) = @ARGV;

my @lib_name = ("ENSDARG","FBgn","ENSGALG","ENSG0","AgaP_","Sjc_","Smp_","csin");
my %speciesInfo = (
	"ENSDARG" => "D._rerio",
	"FBgn" => "D._melanogaster",
	"ENSGALG" => "G._gallus",
	"ENSG0" => "H._sapiens",
	"AgaP_" => "A._gambiae",
	"Sjc_" => "S._japonicum",
	"Smp_" => "S._mansoni",
	"csin" => "C._sinensis",
	"c_elegans" => "C._elegans"
);

open(LIST, "<", $list_file)  || die("** Fail to open $list_file.\n");
my %all4FoldSite;
while(my $line = <LIST>){
	chomp $line;
	my @infile = split(/\s+/,$line);
	my($mapInfo);
	my $clustalw_file = $infile[0];
	$clustalw_file =~ s/id/pro.aln/;
	next if( ! -e $clustalw_file);
	my $alnIn = Bio::AlignIO->new(-format=> 'clustalw', -file=>$clustalw_file);
	my $alnobj = $alnIn->next_aln;
	my $alnobj2 = $alnobj->remove_gaps('-');
    my $alnobj_len = $alnobj->length;
    my $alnobj2_len = $alnobj2->length;
	next if (100 * $alnobj2_len / $alnobj_len) < 50;
	print $clustalw_file,"\n";
	foreach my $seq ($alnobj2 ->each_seq){
	#--------------------------------------------------
	# foreach my $seq ($alnobj ->each_seq){
	#-------------------------------------------------- 
		#--------------------------------------------------
		# print $seq->id(),"\t",$seq->seq(),"\n";
		#-------------------------------------------------- 
		foreach my $lib (@lib_name){
			if ($seq->id() =~ m/$lib/){
				$all4FoldSite{$lib} .= $seq->seq();
				$alnobj2->remove_seq($seq);
				#--------------------------------------------------
				# $alnobj->remove_seq($seq);
				#-------------------------------------------------- 
				last;
			}
		}
	}
	my @tmpSeq = $alnobj2->each_seq();
	#--------------------------------------------------
	# my @tmpSeq = $alnobj->each_seq();
	#-------------------------------------------------- 
	$all4FoldSite{"c_elegans"} .= $tmpSeq[0]->seq();
}
close LIST;

#--------------------------------------------------
# print Dumper(%all4FoldSite);
#-------------------------------------------------- 
open (FOLD4, ">", $out_aln);
my @ids = keys %all4FoldSite;
my $seq_len = length $all4FoldSite{"c_elegans"};
print FOLD4 " ".($#ids+1)." ".$seq_len."\n";
foreach my $id (@ids){
	print FOLD4 $speciesInfo{$id}."\t".substr($all4FoldSite{$id},0)."\n";
	#--------------------------------------------------
	# print FOLD4 $speciesInfo{$id}."\t".$all4FoldSite{$id}."\n";
	#-------------------------------------------------- 
}
close FOLD4;
