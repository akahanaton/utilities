#!/usr/bin/perl 
$format1 = shift;
$format2 = shift || die
	"Usage: $0 format1 format2 input_file";
use Bio::SeqIO;
$in  = Bio::SeqIO->newFh(-format => $format1, -fh => \*ARGV );
$out = Bio::SeqIO->newFh(-format => $format2 );
print $out $_ while <$in>;
