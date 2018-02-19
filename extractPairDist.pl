#!/usr/bin/perl

# This program extract the pairwise distance outputted from Paup showdist. 
# assumes that the distance is in the lower triangle.

my $usage = "Usage: $0 [-h] [-i] [-s seqNames] [infile]\n" .
            "  -h  help\n" .
            "  -i  infinite distance is replaced with \"NA\"\n" .
            "  -s  seq names contained in the seqNames file is used\n" .
            " infile should be an output of PAUP showdist, " .
            "STDIN is used if no infile is given\n";

use Getopt::Std;
getopts('his:') || die "$usage\n";
die "$usage\n" if (defined ($opt_h));

my $state = 0;
my $colNumExpected=1;
my @seqName=();
my $numOTUs = 0;

my @colNum = ();
while(<>) {
    if ($colNumExpected == 1 && /^(\s+\d+)+$/) {
	$state ++;
	s/^\s+//; s/\s+$//;
	@colNum = split;
	$colNumExpected = 0;
	next;
    }
    
    next if ($state == 0);

    if (/\(continued\)\s*$/) {
	$colNumExpected = 1;
	next;
    }

    # When reached here, it is data
    chomp; s/^\s+//; s/\s+$//;
    next if (/^$/);

    s/\s+-\s*$//;  # get rid of the '-' for the diagonal
    if (/\s+\*\d+/) { # Paup attaches "*" for infinite distances
	warn "WARN: Undefined distance found\n$_\n";
	if (defined($opt_i)) {  # substitute with NA
	    s/(\s+)\*(\d+\.?\d*|\.\d+)/$1NA/g;
	} else  {
	    s/(\s+)\*(\d+\.?\d*|\.\d+)/$1$2/g;
	}
# note these pattern match was using (\d), so decimals got screwed up
# fixed 1/30/2004
    }

    my @line = split(/\s+/);
    $rowNum = shift @line;  # get rid of the first column (row number)
    $numOTUs = $rowNum if ($rowNum > $numOTUs);

    my $tmpName = shift(@line);
    if ($state == 1) {       # saving the sequence names
	$seqIndex{$tmpName} = $rowNum - 1;
	push @seqName, $tmpName;
    }

    if (@line > 0) {
	if ($state == 1)  { # first block
	    push @dataArray, [ @line ];
	} else {                         # attach to the end
	    push @{$dataArray[$rowNum - 2]}, @line;
	}
    }
} # end of reading in the data

if ($state == 0) {
    die "End of file reached without data\n";
}

# deal with -s
if (defined($opt_s)) {
    @seqName = ReadSeqNameFile($opt_s);
    @selectedSeqIdx =();
    foreach my $i (@seqName) {
	if (exists ($seqIndex{$i})) {
	    push @selectedSeqIdx, $seqIndex{$i};
	} else {
	    die "The name ($i) in the file $opt_S doesn't exists\n";
	}
    }
} else {
    @selectedSeqIdx = 0..($numOTUs -1);
}

print "seq1\tseqName1\tseq2\tseqName2\tdist\n";
for ($i = 1 ; $i < @selectedSeqIdx; $i++) {
    for ($j = 0; $j < $i; $j ++) {
	my $curRowNum = $i + 1;
	my $curColNum = $j + 1; 
	print "$curRowNum\t$seqName[$i]\t";
	print "$curColNum\t$seqName[$j]\t";
	
	if ($selectedSeqIdx[$i] > $selectedSeqIdx[$j]) {
	    print "$dataArray[$selectedSeqIdx[$i]-1][$selectedSeqIdx[$j]]\t\n";
	} else {
	    print "$dataArray[$selectedSeqIdx[$j]-1][$selectedSeqIdx[$i]]\t\n";
	}
    }
}

sub ReadSeqNameFile {
    my $file = shift;
    my @result = ();
    open(INFILE, "<$file") || die "Can't open $file\n";

    while (<INFILE>) {
        chomp;
        s/#.*$//;    # remove comments (#)
        s/^\s+//; s/\s+$//;
        next if (/^$/);
	push @result, $_;
    }
    close(INFILE);
    return @result;
}
