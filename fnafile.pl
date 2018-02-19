#!/usr/bin/perl

=head1 Name
	输出含有指定 repeat 的序列的 id
	get_repeat_id -i id.file -s seq.fna out.file
=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use Bio::Perl;
use Bio::SeqIO;
use DBI;


my ($idfile,$seqfile,$outfile);
my ($Help,%seqs);

GetOptions(
	"i:s"=>\$idfile,
	"s:s"=>\$seqfile,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV==0 || $Help);
print $seqfile."\n";
my $seq_in = Bio::SeqIO->new(-file => "$seqfile",-format => 'Fasta') or die " can't readin $seqfile \n";

while(my $seq = $seq_in->next_seq() )
{
	$seqs{$seq->id} = $seq->seq();
}


my $outfile = shift;
if ( !open OUT,"> $outfile")
{
	die "failed to open repeats file ($!)";
}


if ( !open ID,$idfile)
{
	die "failed to open repeats file ($!)";
}

while (<ID>)
{
	chomp;
	print OUT ">",$_,"\n";
	print OUT $seqs{$_},"\n";
}
