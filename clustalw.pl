#!/usr/bin/perl -w
use MyMod::DB::GenericDB;
use MyMod::Bio::Seq::mirna;
use strict;
use warnings;
use MyMod::Bio::Tools::miRNA;
use Bio::SeqIO;
use Bio::AlignIO;
use Data::Dumper;


my $seqana=MyMod::Bio::Tools::SeqAna->new();
my $mirna_fact=MyMod::Bio::Tools::miRNA->new();


my $seqin=Bio::SeqIO->new(-file=>$ARGV[2],-format=>"fasta");
my $seed;
my %mirna_seq;
my $seed_info;
while(my $seq=$seqin->next_seq()){
	my ($tmp_acc,$index,$mature,$grade)=split(/_/,$seq->display_id());
	my $acc = $tmp_acc.'_'.$index."_".$grade;
	$seed->{$acc} = substr $mature,1,7;
	$mirna_seq{$acc} = $seq->seq();
}

my $anno_info;
if(defined($ARGV[3])){
	open(ANNO,$ARGV[3]);
	while(my $line=<ANNO>){
		chomp($line);
		my @ary=split(/\t/,$line);
		#my ($acc,undef)=split(/_/,$ary[0]);
		my $acc=$ary[0];
		$anno_info->{$acc}=$ary[1];
	}
}

#read tags information
my $mapping_result=$ARGV[1]; #soap mapping results
my $tags_info;
open(TAG,$mapping_result);
my $total_libs;
while(my $line=<TAG>){
	chomp($line);
	my @ary=split(/\t/,$line);
	my ($lib,$mir,$copy)=split(/_/,$ary[0]);
	#next if $lib=~ m/dsi/;
	next if !($lib=~ m/dm/ && $ary[8]=~ m/dm3/) && !($lib=~ m/ds/ && $ary[8]=~ m/droSim1/) ;
	#next if !($lib=~ m/dm/ && $ary[8]=~ m/dm3/) && !($lib=~ m/ds/ && $ary[8]=~ m/droSim1/) 
	#&& !($lib=~ m/ds/ && $ary[8]=~ m/droSec1/) && !($lib=~ m/DPS/ && $ary[8]=~ m/dp4/)
	#&& !($lib=~ m/DME/ && $ary[8]=~ m/dm3/) && !($lib=~ m/DSI/ && $ary[8]=~ m/droSim1/)
	#&& !($lib=~ m/DPS/ && $ary[8]=~ m/dm3/) && !($lib=~ m/DPS/ && $ary[8]=~ m/droSim1/);
	my $elem={"copy"=>$copy,"seq"=>$ary[1],"acc"=>$lib,"strand"=>$ary[6]};
	push(@{$tags_info->{$ary[8]}},$elem);
	$total_libs->{$lib}=1;
}
my @all_libs=sort(keys(%$total_libs));

my $align_file=$ARGV[0];
my $alnin=Bio::AlignIO->new(-file=>$align_file,-format=>"clustalw", -displayname_flat=>1);

while(my $aln=$alnin->next_aln()){	
	my $major_id;
	my $acc;
	my @seqs=$aln->each_seq();
	($acc,undef)=($seqs[0]->display_id()=~ m/(.*)_(.*)/sg);
	my ($sst_aln)=$mirna_fact->get_sst_aln($aln,1);

	my @aln_taxons;
	foreach my $seq(@seqs){
		my @tmpArr = split(/_/, $seq->display_id());
		push(@aln_taxons, $tmpArr[-1]);
	}
	my %seen = ();
	@aln_taxons = grep { ! $seen{$_}++ } @aln_taxons;

	my $used_taxon;
	my @all_tags;
	foreach my $tmp_taxon(@aln_taxons){
		my $hit_acc=$acc."_".$tmp_taxon;
		next if !defined($tags_info->{$hit_acc});
		$used_taxon->{$tmp_taxon}=1;
		my @tags=@{$tags_info->{$hit_acc}};
		#combined tags
		my $all_copy;
		foreach my $tmptag (@tags){
			$all_copy->{$tmptag->{"seq"}."_".$tmptag->{"strand"}}->{$tmptag->{"acc"}}=$tmptag->{"copy"};
		}
		my @new_tags;
		foreach my $tmpseqid(keys %{$all_copy}){
			my ($tmptag,$tmpstrand)=split(/_/,$tmpseqid);
			my @lib_copys;
			foreach my $lib(@all_libs){
				my $tmpcopy=0;
				$tmpcopy=$all_copy->{$tmpseqid}->{$lib} if defined($all_copy->{$tmpseqid}->{$lib});
				push(@lib_copys,$tmpcopy);
			}
			my $new_copy=join("\t",@lib_copys);
			push(@new_tags,{"acc"=>"","copy"=>$new_copy,"seq"=>$tmptag,"strand"=>$tmpstrand});
		}

		#--------------------------------------------------
		# foreach my $seq_tmp($aln->each_seq()) {
		# 	print $seq_tmp->display_id()."\n";
		# 	print $seq_tmp->seq()."\n";
		# }
		#-------------------------------------------------- 
		
		
		my ($tmpseq)=$aln->each_seq_with_id($hit_acc);
		next if !defined($tmpseq);
		#--------------------------------------------------
		# $tmpseq->end(length($tmpseq->seq()));
		#-------------------------------------------------- 
		$tmpseq->seq(uc($tmpseq->seq()));
		my @tmp_all_tags;
		foreach my $elem(@new_tags){
			my $tmptag=Bio::LocatableSeq->new(
				-display_id=>"tag",
				-seq=>uc($elem->{"seq"}),
				-start=>1,
				-end=>length($elem->{"seq"}),
				-strand=>1
			);
			my ($aln_seq,$tagstart)=$seqana->aln_identity_seq($tmpseq,$tmptag);
			next if !defined($aln_seq);
			my $aln_seq_str=$aln_seq->seq();
			$aln_seq_str=~ s/T/U/sg;
			$aln_seq_str=lc($aln_seq_str) if $elem->{"strand"} eq "-";
			my $tag_line=$aln_seq_str." ".$elem->{"acc"}."_".$tmp_taxon." ".$elem->{"copy"}." ".$elem->{"strand"}."\n";
			push(@tmp_all_tags,{"tag"=>$tag_line,"pos"=>$tagstart});
		}
		my @sort_tags=sort{$a->{"pos"}<=>$b->{"pos"}}@tmp_all_tags;
		grep{push(@all_tags,$_)}@sort_tags;

		if(defined($seed->{$acc})){
			my $tmpseed=Bio::LocatableSeq->new(
				-display_id=>"seed",
				-seq=>uc($seed->{$acc}),
				-start=>1,
				-end=>length($seed->{$acc}),
				-strand=>1
			);
			my $seed_rc = reverse $seed->{$acc};
			$seed_rc =~ tr/ACGTacgt/TGCAtgca/;
			my $tmpseed_rc=Bio::LocatableSeq->new(
				-display_id=>"seed",
				-seq=>uc($seed_rc),
				-start=>1,
				-end=>length($seed_rc),
				-strand=>1
			);
			my ($aln_seq,$tagstart);
			($aln_seq,$tagstart)=$seqana->aln_identity_seq($tmpseq,$tmpseed);
			if (! defined($aln_seq)){
				($aln_seq,$tagstart)=$seqana->aln_identity_seq($tmpseq,$tmpseed_rc);
			}
			next if !defined($aln_seq);
			my $aln_seq_str=$aln_seq->seq();
			$aln_seq_str=~ s/T/U/sg;
			$seed_info->{$acc}= $aln_seq_str;
		}
	}
	
	print ">$acc\n";
	print "$mirna_seq{$acc}\n";
	
	if(defined($anno_info->{$acc})){
		#my @rep_taxons=split(/__/,$rep_info->{$acc});
		$anno_info->{$acc}=~ s/__/; /sg;
		print "Another annotations: ".$anno_info->{$acc}."\n\n";
	}
	
	print "Sequences:\n\n";
	foreach my $other_taxon(@aln_taxons){
		my ($tmpseq)=$aln->each_seq_with_id($acc."_".$other_taxon);
		next if !defined($tmpseq);
		print $tmpseq->seq()."\t\t".$other_taxon."\n";
	}
	print "\n\n";
	
	print "Seed:\n\n";
	if(defined($seed_info->{$acc})){
		print $seed_info->{$acc}
	}
	print "\n\n";

	print "Structures:\n\n";
	foreach my $other_taxon(@aln_taxons){
		#--------------------------------------------------
		# my ($tmpseq)=$sst_aln->each_seq_with_id($acc."_".$other_taxon."_"."sst");
		# next if !defined($tmpseq);
		# print $tmpseq->seq()."\t\t".$other_taxon."\n";
		#-------------------------------------------------- 
		my ($tmpseq)=$sst_aln->{$acc."_".$other_taxon};
		next if !defined($tmpseq);
		print $tmpseq."\t\t".$other_taxon."\n";
	}
	print "\n\n";
	
	print "Reads:\n\n";
	if( $#all_tags >= 0){
		foreach my $tag_elem(@all_tags){
			print $tag_elem->{"tag"};
		}
	}else{
		print "No reads mapped to this miRNA.\n";
	}
	print "\n\n";
	
}
