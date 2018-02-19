#!/usr/bin/perl

use strict;
use Bio::SeqIO;

my $infile = $ARGV[0];
#cds|$protein_id|$gene_name|$locus_tag|$location|$product|$len\n$sequence\n
$infile=~/(\w+)\./;
my $dir=$1;
#system("mkdir $dir");
#open(OUT, ">$dir.cds_info");
#open(OUT2, ">./$dir/cdstable.txt");
#open(OUT3, ">./$dir.cds_location");
#open(OUT4, ">./$dir/cds.fa")||die "can not open dir";
open(OUT5, ">./$dir/spliced.fa")||die "can not open dir";
my $seq_io = Bio::SeqIO->new(-file => "$infile", -format => "genbank");
my $count=0;
while( my $seq_obj = $seq_io->next_seq )
{
my $accession = $seq_obj->accession_number;
my $location;
my $genoseq = $seq_obj->seq;
foreach my $feat_object ( $seq_obj->get_SeqFeatures )
{
	my ($gene_name) = $feat_object->get_tag_values('gene') if $feat_object->has_tag('gene');
	my ($locus_tag) = $feat_object->get_tag_values('locus_tag') if $feat_object->has_tag('locus_tag');
	my ($codon_start) = $feat_object->get_tag_values('codon_start');
	my ($product) = $feat_object->get_tag_values('product') if $feat_object->has_tag('product');
	my ($protein_id) = $feat_object->get_tag_values('protein_id');
	my $strand = $feat_object->strand;
	my $location;
	my $sequence;
	my $translation;

	if ($feat_object->primary_tag eq "CDS")
	{
		next if ($feat_object->has_tag('pseudo'));
		$count++;
		my $loc_strand = $feat_object->strand;
	if ( $feat_object->location->isa('Bio::Location::SplitLocationI') )
				{
					my $endf=0;
					my $startf=0;
					foreach my $loc ( $feat_object->location->sub_Location )
					{
						if($endf)
						{
							if(($endf-$startf+$loc->end-$loc->start)<35)
							{
								print $loc->start,"\t",$loc->end,"\t";
							}
							elsif(($endf-$startf)<35)
							{
								my $str=substr($genoseq,$startf-1,$endf-$startf+1).substr($genoseq,$loc->start-1,35);
							print OUT5 $endf-34,"\t",$endf,"\t",$loc->start,"\t",$loc->start+34,"\n$str\n";
							}
							elsif(($loc->end-$loc->start)<35)
							{
								my $str=substr($genoseq,$endf-35,35).substr($genoseq,$loc->start-1,$loc->end-$loc->start+1);
							print OUT5 $endf-34,"\t",$endf,"\t",$loc->start,"\t",$loc->start+34,"\n$str\n";
							}
							else{
							my $str=substr($genoseq,$endf-35,35).substr($genoseq,$loc->start-1,35);
							print OUT5 $endf-34,"\t",$endf,"\t",$loc->start,"\t",$loc->start+34,"\n$str\n";
						}
						}
						$location = $location . "(" . $loc->start . ".." . $loc->end . ")";
						$endf=$loc->end;
						$startf=$loc->start;
					}
					$sequence = $feat_object->spliced_seq->seq;
					my $len=length($sequence);
					
					#$translation = $feat_object->spliced_seq->translate->seq;
					#chop $translation;	$count++;
					if($strand eq "-1")
					{
						$location=changeloc($location);
					}
				print OUT ">cds$count|$accession|$strand|$protein_id|$gene_name|$locus_tag|$location|$product|$len\n$sequence\n"; 
				print OUT2 ">cds$count\t$accession\t$strand\t$protein_id\t$gene_name\t$locus_tag\t$location\t$product\t$len\t$sequence\n"; 
				print OUT3 ">cds$count|$accession|$strand|$protein_id|$gene_name|$locus_tag|$location|$product|$len\n"; 	
				print OUT4 ">cds$count\n$sequence\n"; 	
				
				}			
				else
				{
					$location = "(" . $feat_object->location->start . ".." . $feat_object->location->end . ")";
					$sequence = $feat_object->seq->seq;
					my $len=length($sequence);
					#$translation = $feat_object->seq->translate->seq;
					#chop $translation;
					if($strand eq "-1")
					{
						$location=changeloc($location);
					}
					print OUT ">cds$count|$accession|$strand|$protein_id|$gene_name|$locus_tag|$location|$product|$len\n$sequence\n";
					print OUT2 ">cds$count\t$accession\t$strand\t$protein_id\t$gene_name\t$locus_tag\t$location\t$product\t$len\t$sequence\n"; 
					print OUT4 ">cds$count\n$sequence\n"; 
					print OUT3 ">cds$count|$accession|$strand|$protein_id|$gene_name|$locus_tag|$location|$product|$len\n"; 
				}
			}
		}
		
	}print $count,"\n";

close(OUT);

sub changeloc()
{
	(my $loc)=@_;
	my	@a=split(/\)/,$loc);
	my $ttt=@a;
	my $loct="";
	for(my $i=$ttt-1;$i>=0;$i--)
	{
		$a[$i]=~/(\d+)..(\d+)/;
		my $h=$1;my $y=$2;
		$a[$i]="(".$y."..".$h.")";
		$loct=$loct.$a[$i];
	}
	return $loct;
}
