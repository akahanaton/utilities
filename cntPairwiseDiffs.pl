#!/usr/bin/perl
# Gives FASTA seq as an argument or STDIN.
# Calculate the pairwise differences

my $usage="Usage: $0 [-t geneticCodeTbl] [-d distanceCorrection] [-f output_format] [alignedDNASeq]\n" .
 "  -c: eliminate codons which include one or more gaps\n" .
 "      This option can be combined with -e\n" .
 "  -e: eliminate base sites where gaps are observed in at least one seq\n".
 "      Gap sites are eliminated pair-wise without this option\n" .
 "  -s: use a file which list names of sequences to be used in the analysis\n".
 "  -a: average differece and average number of sites are printed\n" .
 "  -d: \"kimuraProt\", \"poissonProt\", \"p\".  (default is p)\n" .
 "  -n: do not print seq names\n" .
 "  -f: output format, currently only phylip, PAML style if not specfied";

use Getopt::Std;
getopts('af:d:hces:n') || die "$usage\n";

if (defined($opt_h)) {
    die "$usage\n";
}

my %synoSitesPerCodon = ();
my (%synoDiff, %nonSynoDiff);
my %aaSeq;

die "$usage\n" if (@ARGV > 1);

@ARGV = ('-') unless @ARGV;
my $dnaFile = shift;

# $mode is set "nucleotide" or "protein"
# read in the data and store them in @seqName and @seqDat
ReadInFASTA($dnaFile);
#foreach $i (0..$#seqName) {
#    print ">$seqName[$i]\n$seqDat[$i]\n";
#}

if (defined($opt_c) && $mode eq "protein") {
    die "ERROR: infile was Protein sequences, -c (codon) doesn't make sense\n";
}
# select the sequences (and reorder if necessary)
if (defined($opt_s)) {
    ReorderNameAndDat($opt_s);
}

if (defined($opt_e)) {
    @seqDat = CleanSeqs(@seqDat);
}

# for debugging
# while (($k, $v) = each %aaSeq) { print "**$k** => ++$v++\n"; };
# foreach $k (keys %synoDiff)
#   {print "**$k** => $synoDiff{$k} $nonSynoDiff{$k}\n";};

PrintAllPairDiff();

exit (0);

# takes an arg; name of a file from which data are read
# Then read in the data and store them in @seqName and @seqDat
sub ReadInFASTA {
    my $infile = shift;
    my @line;
    my $i = -1;
    my @result = ();

    open (INFILE, "<$infile") || die "Can't open $infile\n";

    $mode = "nucleotide";
    while (<INFILE>) {
	chomp;
	if (/^>/) {  # name line in fasta format
	    $i++;
	    s/^>\s*//; s/^\s+//; s/\s+$//;
	    #@line = split (/\|/);     # note it takes only the name before |
	    #$line[0] =~ s/\s+$//;
	    $seqName[$i] = $_;
	    $seqDat[$i] = "";
	} else {
	    s/^\s+//; s/\s+$//;
	    s/\s+//g;
	    next if (/^$/);            # skip empty line
	    if (! /^[ATGCUatgcu\-\.\?]+$/) {
		$mode="protein";
	    }
	    $seqDat[$i] = $seqDat[$i] . uc($_);
	}
    }
    close(INFILE);

    if ($mode eq "nucleotide") {
	foreach (@seqDat) {  # CHECK THIS
	    s/U/T/g;   # change U to T
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
    return @result;
}

# Reorder @seqName and @seqDat
sub ReorderNameAndDat {
    my $filename = shift;
    my $i;

    my @selectedSeq = ReadSeqNameFile($filename);

    my %nameIdx=();
    foreach $i (0..$#seqName) {
	if (exists ($nameIdx{$seqName[$i]})) {
	    die "There are multiple seqs with a same name: $seqName[$i]\n";
	} else {
	    $nameIdx{$seqName[$i]} = $i;
	}
    }

    my @orderedName = ();
    my @orderedDat =();
    foreach $i (@selectedSeq) {
	if (exists ($nameIdx{$i})) {
	    push @orderedName, $seqName[$nameIdx{$i}];
	    push @orderedDat, $seqDat[$nameIdx{$i}];
	} else {
	    die "ERROR: $i in seqName files doesn't exist in the data file\n";
	}
    }

    @seqName = @orderedName;
    @seqDat = @orderedDat;
}

sub PrintAllPairDiff {

    my ($i, $j, $seq1, $seq2, $seqName1, $seqName2);

    # calc pair diff
    my @nameVect = ();
    my @cntVect = ();
    my @pDist = ();
    for ($i = 2; $i <= @seqDat; $i++) {
        $seq1 = $seqDat[$i-1]; $seqName1 = $seqName[$i-1];
        for ($j = 1; $j < $i; $j++) {
            $seq2 = $seqDat[$j-1]; $seqName2 = $seqName[$j-1];
	    my @diff = CntPairWiseDiff($seq1, $seq2);
	    push @pDist, $diff[0]/$diff[1];  # proportion of different sites
	    push (@cntVect, join("\t", @diff));
            push @nameVect, "$i\t$seqName1\t$j\t$seqName2";
        }
    }

    if (defined($opt_a)) { # print ave diff and ave number of sites
	my $sumDiff = 0;
	my $sumN = 0;
	foreach $i (@cntVect) {
	    my @record = split(/\t/, $i);
	    $sumDiff += $record[0];
	    $sumN += $record[1];
	}
	my $aveDiff = $sumDiff / @cntVect;
	my $aveN = $sumN / @cntVect;
	print "$aveDiff\t$aveN\n";
    } elsif (defined($opt_d)) {
	if ($opt_d eq "kimuraProt") { # kimura's amino acid
	    @pDist=KimuraProteinDist(@pDist);
	} elsif ($opt_d eq "poissonProt") { # poisson correction for amino acid
	    @pDist=PoissonProteinDist(@pDist);
	}
	if (defined ($opt_f)) { # print some other format
	    if ($opt_f eq "phylip") {
		PrintPhylipDist(@pDist);
	    } else  {
		die "ERROR: $opt_f as the argument of -f is not implemented\n";
	    }
	} elsif (defined ($opt_n)) {
	    map {print "$_\n"} @pDist;
	} else {
	    print "seq1\tseqName1\tseq2\tseqName2\tdist\n";
	    foreach $i (0..$#nameVect) {
		print "$nameVect[$i]\t$pDist[$i]\n";
	    }
	}
    } else { # print pairwise diff's
	if (defined ($opt_n)) {
	    map {print "$_\n"} @cntVect;
	} else {
	    print "seq1\tseqName1\tseq2\tseqName2\tdiff\tn\n";
	    foreach $i (0..$#nameVect) {
		print "$nameVect[$i]\t$cntVect[$i]\n";
	    }
	}
    }
}

sub CntPairWiseDiff {
    my ($seq1, $seq2) = @_;
    my $i;

    # get rid of codons with gaps    
    my @cleaned = CleanSeqs($seq1, $seq2);

    $cleaned[0] =~ s/\t//g;
    $cleaned[1] =~ s/\t//g;
    my @noGapSeq1 = split(//, $cleaned[0]);
    my @noGapSeq2 = split(//, $cleaned[1]);
#    print "$cleaned[0]\t$seq1\n$cleaned[1]\t$seq2\n"; #debug

    # counting number of differences
    my $numBases = @noGapSeq1;
    my $numDiff = 0;
    for ($i = 0; $i < $numBases; $i++) {
	$numDiff++ if ($noGapSeq1[$i] ne $noGapSeq2[$i])
    }

    return ($numDiff, $numBases)
}

sub KimuraProteinDist {
    my $maxP = 0.854101966249684544;
    return (map { ($_ >= $maxP) ? "Inf" : 
		      (($_ == 0) ? 0 : - log(1 - $_ - $_ ** 2 / 5)) } @_);
}

sub PoissonProteinDist {
    return (map { ($_ >= 1) ? "Inf" : (($_ == 0) ? 0 :- log(1 - $_)) } @_);
}

sub PrintPhylipDist {
    my @pVect = @_;
    my $numSeq = @seqName;

    if ($numSeq * ($numSeq-1) != 2 * @pVect) {
	die ("ERROR: number of pVect is weird in PrintPhylip\n" .
	     "$numSeq sequences and pVect has ", $#pVect +1, " elements\n");
    }
    my @distTbl = ();
    my @initVal = Repeat(0, $numSeq);
    for my $i (0 .. ($numSeq-1)) {
	push @distTbl, [ @initVal ];
    }

    my $i = 1; my $j=0;
    for my $d (@pVect) {
	$distTbl[$i][$j] = $distTbl[$j][$i] = $d;

	if ($i - 1 == $j) {
	    $i++; $j = 0;
	} else {
	    $j++;
	}
    }
    if ($i != $numSeq || $j != 0) {
	die "ERROR: conversion to Phylip format didn't work (\$i=$i, \$j=$j)\n";
    }

    print "$numSeq\n";
    for my $i (0 .. ($numSeq-1)) {
	print "$seqName[$i]";
	for my $j (0..$numSeq-1) {
	    print "\t$distTbl[$i][$j]";
	}
	print "\n";
    }
}

# returns an array with the 1st argument repeated n times (2nd arg)
sub Repeat {
    my ($val, $n) = @_;
    my @result = ();
    foreach my $i (0..($n-1)) {
	push @result, $val;
    }
    return @result;
}

# takes multiple string; each string represent a sequence
# Get rid of the sites where at least one seq has a gap (or ambiguious base)
# It returns an array of sequences.  Each seq is tab-delimited
# e.g, returned array:  ("AAA\tGGG\tATA", "AAA\tGGA\tATA", "AAG\tGGG\tATA")
sub CleanSeqs {
    my $minLength = CharLen($_[0]);
    my ($i, $seqNum);
    my @noGapSites = ();
    my @result = ();
    my @codonMat = ();

    # get rid of codons with gaps    
    foreach $i (@_) {
	my @tmpArray = (defined($opt_c)) ? MkTripletArray($i) : split(//, $i);
	push @codonMat, [ @tmpArray ];
	$minLength = @tmpArray if (@tmpArray < $minLength); # find min length
    }

    # identify the sites without any gaps.
    for $i (0 .. $minLength - 1) {
	my $gap = 0;
	if (defined($opt_c)) {
	    for $seqNum (0 .. $#_) {
		if ($codonMat[$seqNum][$i] !~ /^[ACGT]{3}$/) {
		    $gap = 1;
		    last;
		}
	    }
	} elsif ($mode eq "nucleotide") {
	    for $seqNum (0 .. $#_) {
		if ($codonMat[$seqNum][$i] !~ /^[ACGT]{1}$/) {
		    $gap = 1;
		    last;
		}
	    }
	} else { # protein
	    for $seqNum (0 .. $#_) {
		if ($codonMat[$seqNum][$i] =~ /^[\-\?]{1}$/) {
		    $gap = 1;
		    last;
		}
	    }
	}
	push (@noGapSites, $i) if ($gap == 0);
    }

    # select the sites without any gaps
    for $seq (0 .. $#_) {
	my @oneSeq = ();
	foreach $i (@noGapSites) {
	    push @oneSeq, $codonMat[$seq][$i];
	}
	my $tabbedSeq = join "\t", @oneSeq;
	push @result, $tabbedSeq;
    }
    return @result;
}

# this function take two scalars and return the larger value
sub Smaller {
    my ($a, $b) = @_;
    return (($a < $b) ? $a : $b);
}

sub Min {
    my $min = shift;
    foreach my $i (@_) {
	$min = ($i < $min) ? $i : $min;
    }
    return $min;
}

sub Sum {
    my $sum = 0;
    foreach my $i (@_) { 
	$sum += $i;
    }
    return $sum;
}

# count the number of characters in a string
sub CharLen {
    my $string = shift;
    my @charString = split (//, $string);
    return scalar(@charString);
}

# for a given string, it separates into triplets, and return
# the resulting array.
# if the last element is less than a triplet, it will be removed
sub MkTripletArray {
    my $seq = shift;
    $seq =~ s/\s+//g;
    $seq =~ s/(.{3})/$1 /g;
    $seq =~ s/\s+$//;
    my @result = split(/ /, $seq);
    pop @result unless ($result[$#result] =~ /.{3}/);
    return @result;
}

# When one of the codons corresponds to termination, (-1, -1) is returned.
# If two codons are identical, (0,0) is returned even if the codon 
# correspond to termination.

# take a list as the argument and extract the unique elements.
# The order of elements will not be preserved.
sub ExtractUnique {
    my %seen=();
    my @unique = ();

    foreach my $item (@_) {
        push (@unique, $item) unless $seen{$item}++;
    }
    return @unique;
}
