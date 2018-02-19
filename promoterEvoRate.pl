#!/usr/bin/perl -w
#===============================================================================
#         FILE:  dnaDistanceK2P.pl
#
#  DESCRIPTION:  Calculate DNA Distance using Kimura 2-Parameters Model
#        NOTES:  ---
#       AUTHOR:  Lv Yang
#      VERSION:  1.0
#      CREATED:  01/22/2010
#     REVISION:  11/06/2011
#===============================================================================
use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Data::Dumper;

my $version = 1.00;
#necessary arguments
my %opts;
GetOptions(\%opts, "i=s", "o=s", "s=s", "m=s","f=s","p=s","h");
&usage if (!(defined $opts{i} and defined $opts{o} and defined $opts{s} or defined $opts{h} ));

my $payoff = { match   => 3,  # $match,
			mismatch   => -3, # $mismatch,
			gap_open   => 1, # $gap_open,
			gap_extend => 1 # $gap_extend
};

my $aln_input   =   $opts{i};
#--------------------------------------------------
# my $ma_file = $opts{m};
# my $hp_file = $opts{f};
#-------------------------------------------------- 
my $major_taxon =   $opts{s};
#--------------------------------------------------
# my $prj_name =   $opts{p};
#-------------------------------------------------- 
my $output  =   $opts{o};

#create file handle
open OUT, ">", "$output" or die "Cannot create file $output: $!.\n";

my %aln_seq = ();
my @aln_taxons;
my (%identity_info, %kmir_info);
my $alnin=Bio::AlignIO->new(-file=>$aln_input,-format=>"maf",-verbose => -1, -displayname_flat=>1 );
while(my $aln=$alnin->next_aln(-verbose => -1)){
	my @seqs=$aln->each_seq(-verbose => -1);
	# get reference sequence
	my($ref_seq,$refAcc);
	foreach my $seq(@seqs){
		if( $seq->display_id() =~ /$major_taxon/){
			$ref_seq = $seq;
			$refAcc = $seq->display_id();
			last;
		}
	}

	# hairpin conservation
	my $tmp_aln = Bio::SimpleAlign->new();
	$tmp_aln ->add_seq($ref_seq);
	foreach my $seq(@seqs){
		my $curTaxon = $seq->display_id();
		if ($curTaxon ne $major_taxon) {
			push(@aln_taxons, $curTaxon);

			# kmir_info
			my($div, $p, $q) = &div_k2p($ref_seq->seq(), $seq->seq());
			if(defined($div) && $div < 5){
				$kmir_info{$refAcc}{$curTaxon} = $div;
			} else {
				$kmir_info{$refAcc}{$curTaxon} = 'undef';
			}

			$tmp_aln->add_seq($seq);
			my $identity = $tmp_aln->percentage_identity();
			$tmp_aln->remove_seq($seq);
			$identity_info{$refAcc}{$curTaxon} = $identity;
		}
	}
}

my %seen = ();
@aln_taxons = grep { ! $seen{$_}++ } @aln_taxons;

#--------------------------------------------------
# print Dumper(%identity_info);
# print Dumper(%kmir_info);
#-------------------------------------------------- 

foreach my $acc (keys %kmir_info){
	print OUT "$acc:\t";
        foreach my $taxon ( @aln_taxons){
            printf OUT "%s\t", $taxon;
        }
    print OUT "\n";
    print OUT "Kmir of:\t";
        foreach my $taxon ( @aln_taxons){
            if (defined($kmir_info{$acc}{$taxon})){
				printf OUT "%2.2f\t", $kmir_info{$acc}{$taxon};
			}else{
				print OUT "undef\t";
			}
		}
	print OUT "\n";
    print OUT "Identity of:\t";
        foreach my $taxon ( @aln_taxons){
            if (defined($identity_info{$acc}{$taxon})){
                printf OUT "%2.2f\t", $identity_info{$acc}{$taxon};
            }else{
                print OUT "undef\t";
            }
        }
	print OUT "\n";
}


#close file handle
close OUT;
	
#========================================SUB==============================================
sub div_k2p {
    my($seq1, $seq2) = @_;
    my ($transition, $transversion, $n, $k, $p, $q) = (0, 0, 0, 0, 0, 0);
    for my $i (0..length($seq1)-1) {
        my $base1 = substr($seq1, $i, 1);
        my $base2 = substr($seq2, $i, 1);
        next if($base1 eq "-" or $base2 eq "-" or $base1 eq "N" or $base2 eq "N");
        $n++;
        $base1 =~ tr/atcgu/ATCGT/;
        $base1 =~ tr/atcgu/ATCGT/;
        next if($base1 eq $base2); 
        if($base1 eq "A" and $base2 eq "G") {
            $transition++;
        } 
        elsif($base1 eq "T" and $base2 eq "C") {
            $transition++;
        } 
        elsif($base1 eq "C" and $base2 eq "T") {
            $transition++;
        } 
        elsif($base1 eq "G" and $base2 eq "A") {
            $transition++;
        }
        else {
            $transversion++;
        }
    }

    #1) Two sequences are not available
    if(!$n) {
        return (0, 0, 0);
    }

    #2) K2P distance
    $p = $transition / $n;
    $q = $transversion / $n;
    if(1-2*$p-$q > 0 and 1-2*$q >0 ) {
        $k = (-1/2)*log((1-2*$p-$q)*sqrt(1-2*$q));
        return (abs($k), $p, $q);
    } else {

    #3) Two sequences are far diverged and cannot calculate K2P distance
        return (5, -1, -1);
    }
}

#========================================SUB==============================================

sub usage{
    print <<"USAGE";
Usage:
	$0 -i <input file> -s <reference species name> -o <output file>
options:
	-i  input file (maf alighment file)
	-s  reference species name in your MAF file (e.g. dm3)
	-o  output file
	-h  help
USAGE
exit(1);
}

sub rev_com_seq {
	my $seq = shift ;
	$seq =~ tr/atcgnATCGN/tagcnTAGCN/ ;
	$seq =reverse($seq);
	return($seq);
}

sub get_residue_number{
	my ($raw_hairpin_seq, $mat_seq) = @_;
	my $tag;
	my $hairpin_seq = lc($raw_hairpin_seq->seq());
	$hairpin_seq =~ s/\-//g;
	my $hairpin_lng = length $hairpin_seq;
	my ($mature_beg, $mature_end);
	my $tmp_pos;
	$tmp_pos = index($hairpin_seq, $mat_seq);
	if ($tmp_pos < 0){
		my $mat_seq_rec = rev_com_seq($mat_seq);
		$tmp_pos = index($hairpin_seq, $mat_seq);
		if ($tmp_pos < 0){
			print STDERR "here\t$mat_seq\t$hairpin_seq\n";
			return (undef, undef, undef);
		}
	}
	if( $tmp_pos> ($hairpin_lng / 2)){
		$tag = '3p';
	}else{
		$tag = '5p';
	}

	$mature_beg = $tmp_pos + $raw_hairpin_seq->start;
	$mature_end = $mature_beg + length $mat_seq;
	return ($mature_beg, $mature_end, $tag);
}

#========================================END==============================================
