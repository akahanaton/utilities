#!/usr/bin/perl

#
# sharcgs.pl
#
# This program assembles very short sequence reads into genomic contigs.
# Further documentation is found at the end of this file, call the 
# script with no command-line parameters to read it.
#
# written by Juliane Dohm and Claudio Lottaz
# Max Planck Institute for Molecular Genetics (Berlin, Germany)
# September 2006 - February 2007
#

use strict;
use POSIX qw(strftime);

my $version = "1.2.11";
$| = 1; # force flush after any output
my $debug = 0; # setting this to 1 causes additional file output

########## read command line

# just help if no command-line arguments are specified
if (scalar(@ARGV) == 0) { 
    system("perldoc $0"); 
    exit(0); 
}

my $usage = "Usage: $0 [options] <read-file>\n" .
    "Options: -d <max-dangling-end>\n" .
    "         -e <expected-error-rate> [%]\n" .
    "         -g <gap-span>\n" .
    "         -l <min-contig-length>\n" .
    "         -n <target-length> [kb]\n" .
    "         -o <output-file>\n" .
    "         -p <prefix-trie-depth>\n" .
    "         -q <quality-file>\n" .
    "         -r <read-repetition-requirement>[,<read-repetition-requirement>]+\n" .
    "         -t <end-trimming>\n";

sub checkInteger {
    my($arg) = @_;
    if ($arg =~ m/^-?\d+$/) { return($arg); }
    die("Integer expected instead of '$arg'\n$usage");
}

sub checkIntegerVector {
    my($arg) = @_;
    if ($arg =~ m/^\d+(,\d+)*$/) { return(split(/,/, $arg)); }
    die("Set of integers expected instead of '$arg'\n$usage");
}

sub checkFloat {
    my($arg) = @_;
    if ($arg =~ m/^\d+\.?\d*$/) { return($arg); }
    die("Integer expected instead of '$arg'\n$usage");
}

# defaults
my $readFile          = "";
my $qualityFile       = "none";
my @minConfirm        = ();      # if not provided, adapted to input data
my $prefixTrieDepth   = 0;       # if not provided, adapted to input data
my $maxGapSpan        = -1;      # if not provided, adapted to input data
my $minGapSpan        = 3;
my $minContigLen      = 50;
my $trimming          = 0;
my $errRate           = "estimated"; # if not provided, guessed from read frequencies
my $targetLen         = "estimated"; # if not provided, guessed from read frequencies
my $force             = 0;
my $errMin            = 0.1;
my $errMax            = 2.0;
my $maxDangling       = 2;
my $contigFile        = "contigs.seq";
my $skipAssembly      = 0; # undocumented feature

my $arg;
while ($arg = shift()) {
    if ($arg =~ m/^-/) {
	if ($arg eq "-d") { $maxDangling       = checkInteger(shift());       next; }
	if ($arg eq "-e") { $errRate           = checkFloat(shift());       next; }
	if ($arg eq "-f") { $force             = 1;;                          next; }
	if ($arg eq "-g") { $maxGapSpan        = checkInteger(shift());       next; }
	if ($arg eq "-l") { $minContigLen      = checkInteger(shift());       next; }
	if ($arg eq "-n") { $targetLen         = checkInteger(shift());       next; }
	if ($arg eq "-o") { $contigFile        = shift();                     next; }
	if ($arg eq "-p") { $prefixTrieDepth   = checkInteger(shift());       next; }
	if ($arg eq "-q") { $qualityFile       = shift();                     next; }
	if ($arg eq "-r") { @minConfirm        = checkIntegerVector(shift()); next; }
	if ($arg eq "-s") { $skipAssembly      = 1;;                          next; }
	if ($arg eq "-t") { $trimming          = checkInteger(shift());       next; }
	die("Unknown switch '$arg'\n$usage");
    } 
    if ($readFile ne "") { 
	print("Unexpected token on command-line: '$arg', ignored.\n"); 
	next; 
    }
    $readFile = $arg;
}

my $minConfirm = "automatic";
if (@minConfirm) { $minConfirm = join(", ", @minConfirm); }
if ($minConfirm eq "0") { $minConfirm = "not required"; }
my $gapSpan = "automatic";
if ($maxGapSpan > -1) { $gapSpan = $maxGapSpan; }
my $prefixSetting = "automatic";
if ($prefixTrieDepth > 0) { $prefixSetting = $prefixTrieDepth;}

print("sharcgs version $version\n",
      "Using this software implies acceptance of the terms of the GNU General\n",
      "Public License (see http://www.gnu.org/licenses/gpl-3.0.txt)\n",
      "Cite Dohm et al., Genome Res. 2007, 17:1697-1706, PMID 17908823\n",
      "----------------------------------------------------------------------\n",
      "Read sequences                : $readFile\n",
      "Quality measures              : $qualityFile\n",
      "Required read confirmation    : $minConfirm\n",
      "Error rate                    : $errRate\n",
      "Target sequence lenggh        : $targetLen [kb]\n",
      "Gap span                      : $gapSpan\n",
      "Prefix trie depth             : $prefixSetting\n",
      "Required contig length        : $minContigLen\n",
      "End trimming                  : $trimming\n",
      "Output file                   : $contigFile\n",
      "----------------------------------------------------------------------\n");
if ($qualityFile eq "none") { $qualityFile = ""; }

########## read input reads

print(strftime("%b %e %H:%M:%S", localtime), 
      " - Read input data................");

if (!(-s $readFile)) { die("\nRead sequences in '$readFile' not found.\n"); }
if ($qualityFile) {
    if (!(-s $qualityFile)) { die("\nQuality measures in '$qualityFile' not found.\n"); }
}

my($q, @qs, $i, $j);
my $nrawReads = 0;
my $readLen = 0;
my %rawReads = ();
open(INFILE, "<$readFile") or die "$!\n";
if ($qualityFile) {
    open(QFILE, "<$qualityFile") or die "$!\n";
}
while(<INFILE>) {
    chomp();
    if ((m/^>/) | (m/\./)) { 
	if ($qualityFile) { $q = <QFILE>; }
	next; 
    }
    if (m/@/) { die("'@' found in read '$_'.\n"); next; }
    if (m/\t/) { 
	@_ = split("\t", $_);
	$_ = @_[4]; 
    }
    if ($readLen == 0) { $readLen = length($_) - $trimming;
    } else {
	if (length($_) != $readLen + $trimming) { die("Read '$_' from '$readFile' is not " . $readLen . " nts long"); }
    }
    $_ = substr($_, 0, $readLen);
    if ($qualityFile) {
	$q = <QFILE>;
	chomp($q);
	if ($q =~ m/\t/) {
	    @qs = split(/\t/, $q);
	    for ($i = 0; $i < $readLen; $i++) {
		my $bestq = -30;
		foreach $q (split(" ", $qs[$i])) {
		    if ($bestq < $q) { $bestq = $q; }
		}
		$qs[$i] = $bestq;
	    }
	} else {
	    @qs = split(/ /, $q);
	    while(scalar(@qs) > $readLen) { pop(@qs); }
	}
	if (scalar(@qs) < $readLen) { die("Found" . scalar(@qs) . "quality values for read '$_'.\n"); }
	if (exists($rawReads{$_})) {
	    $rawReads{$_}[0]++;
	    # choose max quality value, confirmation on same strand not independent
	    for ($i = 1; $i <= $readLen; $i++) {
		if ($rawReads{$_}[$i] < $qs[$i-1]) {
		    $rawReads{$_}[$i] = $qs[$i-1];
		}
	    }
	} else {
	    unshift(@qs, 1);
	    for ($i = 0; $i <= $readLen; $i++) { $rawReads{$_}[$i] = $qs[$i]; }
	}
    } else { # no quality file
	$rawReads{$_}[0]++;
    }
    $nrawReads++;
}
close(INFILE);
if ($qualityFile) { close(QFILE); }
print(sprintf("%8d", $nrawReads), " reads (", $readLen, "nt)\n");

########## generate reverse strand reads

sub revComp {
    my($x) = @_;
    $x = reverse($x);
    $x =~ tr/ATCG/TAGC/;
    return($x);
}

print(strftime("%b %e %H:%M:%S", localtime),
      " - Add reverse complements........");
my($mer,$revmer);
my %reads = ();
foreach $mer (keys(%rawReads)) {
    @qs = @{$rawReads{$mer}};
    $revmer = revComp($mer);
    if(exists($reads{$mer})) {
	$reads{$mer}[0] += @qs[0]; 
	for ($i = 1; $i < scalar(@qs); $i++) { $reads{$mer}[$i] += @qs[$readLen - $i + 1]; }
    } else { 
	# the same confirmation level is used for a read and its complement
	for ($i = 0; $i < scalar(@qs); $i++) { $reads{$mer}[$i] = $qs[$i]; }
	$reads{$revmer} = $reads{$mer};
    }
}
print(sprintf("%8d", scalar(keys(%reads))), " unique reads\n");

if ($debug) {
    my $uniqueFile = $contigFile;
    $uniqueFile =~ s/[^\.]+$/unique/;
    open(OUTFILE, ">$uniqueFile");
    foreach $mer (keys(%reads)) {
	$q = $reads{$mer};
	print(OUTFILE "$mer\t" . join(" ", @$q) . "\n");
    }
    close(OUTFILE);
}

########## Estimate sequence length and error rate
                

if ($targetLen eq "estimated") {
    print(strftime("%b %e %H:%M:%S", localtime),
	  " - Estimate sequence length/error rate...");
}

sub logBinomDensity {
    my($n, $size, $p) = @_;
    my($res, $i);
    my $pi = 3.141592653589793116;

    # choose($size, $n)
    $res = 0;
    if ($n < 9) { # exact
	for ($i = 0; $i < $n; $i++) {
	    $res += log($size - $i) - log($i+1)
	    }
    } else { # stirling approximation
	$res = log(sqrt(2*$pi*$size)) + $size*(log($size) - 1);
	$res -= log(sqrt(2*$pi*$n)) + $n*(log($n) - 1);
	$res -= log(sqrt(2*$pi*($size-$n))) + ($size-$n)*(log($size-$n) - 1);
    }

    # p^n (1-p)^(size-n)
    if (($p != 0) & ($p != 1)) { $res += $n * log($p) + ($size-$n)*log(1-$p); }
    return($res);
}

# find the number of reads observed n times
my @readConfirm = ();
foreach (keys(%reads)) { $readConfirm[$reads{$_}[0]]++; }

# see how often the most frequent reads are generated
my $i;
my $bestI = 8;
my $bestMulti = 0;
my($er, $erhat, $nhat, $n, $dist, $binSearch, @thisP, $thisBestP, $thisBestN, $tmpP, $tmpN, $p);

# find the typical number of read confirmation
for ($i = 8; $i < scalar(@readConfirm); $i++) {
    if ($bestMulti < $readConfirm[$i]) {
	$bestMulti = $readConfirm[$i];
	$bestI = $i;
    }
}

# find the minimum sequence length (number of unique reads)
my $tmpMin = $bestI - 5;
my $tmpMax = $bestI + 5;
my $nmin = 0;
for ($i = $tmpMin; $i <= scalar(@readConfirm); $i++) { $nmin += $readConfirm[$i]/2; }

if ($targetLen ne "estimated") {
    $nhat = $targetLen * 1000;
} else {
    if ($nmin < 1000) { # we don't expect estimation of sequence length feasible
	print(". undetermined\n");
	$nhat = 1000;
	$erhat = 0.1;
    } else { # fit the length/error model to the most frequent reads
	my $bestP = -1e307;
	my $oldP  = -1e308;
	if ($errRate ne "estimated") { 
	    $errMin = $errRate;
	    $errMax = $errRate;
	}
	for ($er = $errMin; $er <= $errMax; $er += 0.01) {
	    # find the first maximum in the lieklihood
	    $p = (1 - $er/100) ** $readLen;
	    $dist  = 2**18;
	    $binSearch = 0;
	    $n = $nmin + $dist;
	    while(($n > 0) & ($n <= 5000000) & (abs($dist) > 1)) {
		for ($tmpN = 0; $tmpN < 2; $tmpN++) {
		    $thisP[$tmpN] = 0;
		    for ($i = $tmpMin; $i <= $tmpMax; $i++) {
			$tmpP = logBinomDensity($i, $nrawReads, $p/($n+$tmpN));
			$thisP[$tmpN] += logBinomDensity($readConfirm[$i]/2, $n+1, exp($tmpP));
		    }
		}
		if ($thisP[1] > $thisP[0]) { # ascending
		    if ($binSearch) { $dist = abs($dist)/2; }
		} else { # descending
		    $binSearch = 1;
		    $dist = -abs($dist)/2;
		}
		$n += $dist;
	    }
	    if ($thisP[0] > $bestP) {
		$bestP = $thisP[0];
		$nhat = $n;
		$erhat = $er;
	    }
	}
	print(sprintf("%6.1fkb/%4.2f%%\n", $nhat/1000, $erhat));
    }
}

########## fine-tune unspecified parameters

# configure gap span
if ($maxGapSpan > -1) { # -g option provided
    $prefixTrieDepth = $maxGapSpan + 1;
    $minGapSpan = $maxGapSpan;
}

# configure prefix trie according to sequence length
my $maxMissingRatio = 0.5;
if ($prefixTrieDepth) { # -p or -g option was given
    $maxMissingRatio = exp(-log($nhat)/$prefixTrieDepth); 
} else { 
    $prefixTrieDepth = int(-log($nhat)/log(0.5)) + 1; 
}
print(strftime("%b %e %H:%M:%S", localtime),
      sprintf(" - Prefix trie depth to handle %2d%% missing reads.....%2d\n", 
	      100*$maxMissingRatio, $prefixTrieDepth));
my $prefixLen = $readLen - $prefixTrieDepth - 1;
$maxGapSpan = $prefixTrieDepth - 1;

# configure read-confirmation-level (min repeat or min worst quality value)
my $provision = 1000; # extra reads allowed in first filter step
my $minQuality = 1; # used when no quality file is provided
if ($qualityFile) {
    $minQuality <- 15;
    @readConfirm = ();
    foreach $mer (keys(%reads)) {
	# store worst quality value in $reads{$mer}[0] (instead of #repeats)
	$reads{$mer}[0] = $reads{$mer}[1];
	for ($i = 2; $i < scalar(@qs); $i++) {
	    if ($reads{$mer}[0] > $reads{$mer}[$i]) { $reads{$mer}[0] = $reads{$mer}[$i]; }
	}
	if ($reads{$mer}[0] < 1) {$reads{$mer}[0] = 1; }
	$readConfirm[$reads{$mer}[0]]++;
    }
}
my $cumsum = 0;
if (!@minConfirm) { # -r option not provided 
    for ($i = scalar(@readConfirm); ($i > $minQuality) & ($cumsum/$nhat < 1-$maxMissingRatio); $i--) {
	$cumsum += $readConfirm[$i]/2; # /2 (both strands)
    }
    push(@minConfirm, $i); $i--;
    while($i > $minQuality) {
	$cumsum += $readConfirm[$i]/2; # /2 (both strands)
	if ($cumsum < $nhat+$provision) { push(@minConfirm, $i); }
	$i--;
    }
    print(strftime("%b %e %H:%M:%S", localtime),
	  sprintf(" - Confirmation required for %2d%% missing reads....", 100*$maxMissingRatio)); 
    print(sprintf("%2d-%2d\n", $minConfirm[$#minConfirm], $minConfirm[0])); 
    if (scalar(@minConfirm) > 3) {
	@minConfirm = ($minConfirm[0], $minConfirm[int(scalar(@minConfirm)/2)], $minConfirm[$#minConfirm]);
    }
}

########## fill and access prefix trie

my(%prefixTrie, $ncells, $prefix, @p, $p, $r, $rold, $k);
sub addRead {
    my($read) = @_;
    $prefix = substr($read, 0, $prefixLen);
    @p = split(//, substr($read, $prefixLen));
    $p = pop(@p);
    if (!exists($prefixTrie{$prefix})) { $prefixTrie{$prefix} = {}; $ncells++; }
    $r = $prefixTrie{$prefix};
    while(@p) {
	$k = shift(@p);
	$rold = $r;
	if (!exists($$r{$k})) { $$r{$k} = {}; $ncells++; }
	$r = $$r{$k};
    }
    if (ref($$rold{$k})) { $$rold{$k} = $p; }
    else { 
	if ($$rold{$k} !~ m/$p/) {
	    $$rold{$k} .= $p; 
	}
    }
}


sub removeSuffix {
    my($suffix, $ref) = @_;
    if (length($suffix) == 1) {
	$ref =~ s/$suffix//;
	return($ref);
    }
    my $key = substr($suffix, 0, 1);
    if (!exists($$ref{$key})) { return($ref); }
    $$ref{$key} = removeSuffix(substr($suffix,1), $$ref{$key});
    if (!$$ref{$key}) { delete($$ref{$key}); }
    if (scalar(keys(%$ref)) == 0) { return(""); }
    return($ref);
}
sub removeRead {
    my($read) = @_;
    $prefix = substr($read, 0, $prefixLen);
    if (exists($prefixTrie{$prefix})) { 
	my $suffix = substr($read, $prefixLen);
	my $keep = removeSuffix($suffix, $prefixTrie{$prefix});
	if (!$keep) { delete($prefixTrie{$prefix}); }
    }
}


# The following function matches its only parameter to the prefix tree 
# stored in the global variable %prefixTrie. It also looks up $readLen 
# and $prefixLen in the global environment. If the parameter is longer 
# than $minOverlap, this function does nothing useful. Otherwise it
# returns the extension needed to complete the prefix given as parameter
# to yield a complete read stored in the prefix trie. The extension 
# contains a "@" if the extension is ambiguous and may be emtpy if 
# no read is matched.
sub matchPrefix {
    my($prefix) = @_;
    my $prePrefix = substr($prefix, 0, $prefixLen);
    if (!exists($prefixTrie{$prePrefix})) { return(""); }    
    my $p = $prefixTrie{$prePrefix};
    $prefix = substr($prefix, $prefixLen);
    while($prefix) {
	$_ = substr($prefix, 0, 1);
	if (!exists($$p{$_})) { return(""); }
	$p = $$p{$_};
	$prefix = substr($prefix, 1);
    }
    my $extension = "";
    while(ref($p)) {
	if (scalar(keys(%$p)) > 1) { return($extension . "@"); }
	$extension .= (%$p)[0];
	$p = $$p{(%$p)[0]};
    }
    if (length($p) > 1) { return($extension . "@"); }
    return($extension . $p);
}

########## start long loop

my(%remainingReads, %filteredReads, %contigs, %mergedContigs, %shortContigs);
foreach $minConfirm (@minConfirm) {
    if ($minConfirm == 0 ) {
	%filteredReads = %reads;
    } else {
	%filteredReads = ();

	########## filter reads with too little confirmation

	if ($qualityFile) {
	    print(strftime("%b %e %H:%M:%S", localtime),
		  sprintf(" - Remove low quality reads (Q<%2d)..", $minConfirm));
	} else {
	    print(strftime("%b %e %H:%M:%S", localtime),
		  " - Remove rare reads (<$minConfirm times).....");
	}
	if ($debug) {
	    my $dropFile = $contigFile;
	    $dropFile =~ s/[^\.]+$/$minConfirm\_rare_dropped/;
	    open(OUTFILE, ">$dropFile");
	}
	foreach (keys(%reads)) {
	    if ($reads{$_}[0] < $minConfirm) { 
		if ($debug) { print(OUTFILE "$_\n"); }
	    }  else { $filteredReads{$_} = 1; }
	}
	close(OUTFILE);
	print(" kept ", sprintf("%7d", scalar(keys(%filteredReads))), " reads\n");
    }

    ########## filter reads with no matching partners

    my $read;
    print(strftime("%b %e %H:%M:%S", localtime),
	  " - Store quality reads to prefix trie....");
    $ncells = 0;
    %prefixTrie = ();
    foreach $read (keys(%filteredReads)) { addRead($read); }
    print(sprintf("%8d", $ncells), " nodes\n");

    my($minOverlap, $guessGapSpan);
    my $nreads = scalar(keys(%filteredReads));
    if ($minConfirm > 0) {
	print(strftime("%b %e %H:%M:%S", localtime),
	      " - Remove reads with no partners....");

	if ($targetLen eq "estimated") {
	    $minOverlap = int($readLen/2) + 1;
	} else {
	    # the following relies on the sequence length. This is estimated, it is 
	    # unlikely to work well on real data.
	    if ($nreads/2 < $nhat) { # otherwise we under estimated nhat/included error reads
		$p = $nreads/2 / $nhat; # /2 since %reads is from both strands
		$guessGapSpan = int(-(log($nhat) + log($p)) / log(1-$p));
	    } else { $guessGapSpan = $minGapSpan; }
	    $minOverlap = $readLen - $guessGapSpan - 1;
	    if ($minOverlap < int($readLen/2)+1) { $minOverlap = int($readLen/2)+1; }
	}
	my @singlets = ();
	foreach $read (keys(%filteredReads)) {
	    my @extensions = ();
	    for (my $strand = 0; $strand < 2; $strand++) {
		my $overlap = $read;
		while((length($overlap) > $minOverlap) & !$extensions[$strand]) {
		    $overlap = substr($overlap, 1);
		    $extensions[$strand] = matchPrefix($overlap);
		}
		$read = revComp($read);
	    }
	    if (!$extensions[0] | !$extensions[1]) { push(@singlets, $read); }
	}
	foreach $read (@singlets) {	
	    delete $filteredReads{$read};
	    removeRead($read);
	}
	print(" kept ", sprintf("%7d", scalar(keys(%filteredReads))), " reads\n");

	if ($debug) {
	    # reads removed by partner check
	    my $dropFile = $contigFile;
	    $dropFile =~ s/[^\.]+$/$minConfirm\_singlet_dropped/;
	    open(OUTFILE, ">$dropFile");
	    foreach $read (@singlets) { print(OUTFILE "$read\n"); }
	    close(OUTFILE);

	    # reads used for assembly
	    my $filterFile = $contigFile;
	    $filterFile =~ s/[^\.]+$/$minConfirm\_filtered/;
	    open(OUTFILE, ">$filterFile");
	    my $ctr = 1;
	    foreach(keys(%filteredReads)) { print(OUTFILE ">read$ctr\n$_\n"); $ctr++; }
	    close(OUTFILE);
	}
    }
    if ($skipAssembly) { next; }

    ########## Estimate gap start

    # choose initial gap span such that we expect 1 or less larger
    # holes according to Poisson-distribution with intensity nhat*gap^gapSpan

    print(strftime("%b %e %H:%M:%S", localtime),
	  " - Guess read fraction/best gap span.......... ");
    $nreads = scalar(keys(%filteredReads));
    if ($nreads/2 < $nhat) { # otherwise we under estimated nhat/included error reads
	$p = $nreads/2 / $nhat; # /2 since %reads is from both strands
	$guessGapSpan = int(-(log($nhat) + log($p)) / log(1-$p));
    } else {
	$p = 1;
	$guessGapSpan = 0;
    }
    print(sprintf("%4.1f%%/%2d\n", 100*$p, $guessGapSpan));
    if ($guessGapSpan < $minGapSpan) { $guessGapSpan = $minGapSpan; }
    if (($guessGapSpan > $maxGapSpan) & (!$force)){ 
	print("                  Too many gaps expected, the assembly is skipped...\n",
	      "----------------------------------------------------------------------\n");
	next;
    }

    ########## assemble contigs

    print(strftime("%b %e %H:%M:%S", localtime),
	  " - Assenble contigs... \n");

    sub progressBar {
	my(($progress, $nreads, $readsRef)) = @_;
	my $remaining = scalar(keys(%$readsRef));
	my $newProgress = int(39*($nreads - $remaining)/$nreads);
	for (my $i = $progress; $i < $newProgress; $i++) { print("-"); }
	return($newProgress);
    }

    %contigs = ();
    print("                Span|0%.......:.........:........:......100%| \#Contigs\n");

    my $contig;
    my $minOver = $readLen - $maxGapSpan - 1;
    my $maxOver = $readLen - $minGapSpan - 1;
    $minOverlap = $readLen - $guessGapSpan - 1;
    if ($minOverlap < $minOver) { $minOverlap = $minOver }
    my $dangling = $maxDangling + 1; # make sure to kick off
    if ($maxDangling < 0) { 
	$minOver = $minOverlap; 
	$maxOver = $minOverlap; 
    }
    while(($minOverlap >= $minOver) & ($minOverlap <= $maxOver) & ($dangling > $maxDangling)) {
	%contigs = ();
	%remainingReads = %filteredReads;
	print(strftime("%b %e %H:%M:%S", localtime),
	      sprintf(": %2d |", $readLen - $minOverlap - 1));
	$dangling = 0;
	my $progress = 0;
	while($contig = (%remainingReads)[0]) {
	    delete($remainingReads{$contig});
	    my $thisDangling = 0;	    $progress = progressBar($progress, $nreads, \%remainingReads); 

	    for (my $strand = 0; $strand < 2; $strand++) {
		my $extension = "start";
		while($extension) {
		    # find shortest possible extension
		    $extension = "";
		    my $overlap = substr($contig, -$readLen);
		    while((length($overlap) > $minOverlap) & !$extension) {
			$overlap = substr($overlap, 1);
			$extension = matchPrefix($overlap);
		    }
		    if (!$extension) { $thisDangling++; }
		    $extension =~ s/\@$//;

		    # avoid reusing reads
		    my $usedRead = $overlap . $extension;
		    if (!exists($remainingReads{$usedRead}) & (length($usedRead)==$readLen)) { 
			$extension = ""; 
		    }

		    # check for local ambiguities
		    my($tmpSnip, $tmpOverlap, $tmpEnd);
		    if ($extension) { $tmpSnip = substr($contig, -$readLen).$extension; }
		    for (my $tmpStrand = 0; $tmpStrand<2 & $extension ne ""; $tmpStrand++) {
			$tmpOverlap = substr($tmpSnip, 0, $minOverlap);
			$tmpEnd = substr($tmpSnip, $minOverlap);
			while($tmpEnd ne "" & $extension ne "") {
			    my $tmpExt = matchPrefix($tmpOverlap);
			    if (($tmpExt !~ m/$tmpEnd/) & ($tmpEnd !~ m/$tmpExt/)) { 
				$extension = ""; 
			    }
			    $tmpOverlap = substr($tmpOverlap, 1) . substr($tmpEnd, 0, 1);
			    $tmpEnd = substr($tmpEnd, 1);
			}
			$tmpSnip = revComp($tmpSnip);
		    }

		    # extend the current contig
		    if ($extension) { 
			$contig .= $extension;
			if (length($usedRead) == $readLen) {
			    delete $remainingReads{$usedRead};
			    delete $remainingReads{revComp($usedRead)};
			}
		    }
		}
		# 5'-extension on direct strand is equivalent to 3'-extension on reverse strand
		$contig = revComp($contig);
	    }
	    # store new contig
	    if (length($contig) >= $minContigLen) { 
		$dangling += $thisDangling;
		$contigs{$contig}++;
	    } else { $shortContigs{$contig}++; }
	}
	print(sprintf("|%5d/%3d\n", scalar(keys(%contigs)), $dangling));
	$minOverlap--;
    }

    # merge into former contigs, sort for decreasing length first
    my @contigs = sort({length($b) <=> length($a)} keys(%contigs));
    my @mergedContigs = sort({length($b) <=> length($a)} keys(%mergedContigs));

    # replace merged contigs by new ones which contain them
    my $j = 0;
    for(my $i = 0; $i < scalar(@contigs); $i++) {
	my $len = length($contigs[$i]);
	while((length($mergedContigs[$j])>$len) & ($j<scalar(@mergedContigs))) {$j++;}
	for (my $jj = $j; $jj < scalar(@mergedContigs); $jj++) {
	    if ((index($contigs[$i], $mergedContigs[$jj]) >-1) |
		(index($contigs[$i], revComp($mergedContigs[$jj])) >-1)) {
		delete $mergedContigs{$mergedContigs[$jj]};
	    }
	}
    }
    
    # merge new contigs to merged contigs
    my $pos;
    while($contig = pop(@contigs)) {
	my $startProbe = substr($contig, 0, $readLen);
	my $startProbeRev = revComp($startProbe);
	my $endProbe = substr($contig, -$readLen);
	my $endProbeRev = revComp($endProbe);
	$pos = -1; 
	foreach my $merged (keys(%mergedContigs)) {
	    my $extLen = 0;
	    my $newMerged = "";
	    if (($pos = index($merged, $startProbe)) > -1) {
		$extLen = length($contig) + $pos - length($merged);
		if ($extLen > 0) { $newMerged = $merged . substr($contig, -$extLen); }
	    } elsif (($pos = index($merged, $startProbeRev)) > -1) {
		$extLen = length($contig) - $readLen - $pos;
		if ($extLen > 0) {$newMerged = substr(revComp($contig),0,$extLen).$merged;}
	    } elsif (($pos = index($merged, $endProbe)) > -1) {
		$extLen = length($contig) - $readLen - $pos;
		if ($extLen > 0) { $newMerged = substr($contig, 0, $extLen) . $merged; }
	    } elsif (($pos = index($merged, $endProbeRev)) > -1) {
		$extLen = length($contig) + $pos - length($merged);
		if ($extLen > 0) { $newMerged = $merged.substr(revComp($contig),-$extLen);}
	    }
	    if ($newMerged) {
		delete $mergedContigs{$merged};
		push(@contigs, $newMerged);
	    }
	    if ($pos >= 0) { last; }
	}
	if ($pos < 0) { $mergedContigs{$contig} = 1; }
    }
    print(strftime("%b %e %H:%M:%S", localtime),
	  " - Merged contigs generated...................... ", 
	  sprintf("%5d\n", scalar(keys(%mergedContigs))),
	  "----------------------------------------------------------------------\n");
} # next configuration $minConfirm/$neighborThreshold
if ($skipAssembly) { exit(0); }

########## output file containing long contigs

print(strftime("%b %e %H:%M:%S", localtime),
      " - Generate sequence output.............. ");

# write resulting contigs
my $contigCtr = 0;
my @lens = ();
open(OUTFILE, ">$contigFile");
foreach my $contig (keys(%mergedContigs)) {
    $contigCtr++;
    push(@lens, length($contig));
    print(OUTFILE ">contig$contigCtr\n$contig\n");
}
close(OUTFILE);
print(sprintf("%5d contigs\n", $contigCtr));

# generate quality output file for contigs
if ($qualityFile) {
    print(strftime("%b %e %H:%M:%S", localtime),
	  " - Generate quality output......... ");
    $contigCtr = 0;
    my @qsrev = ();
    my $contigQFile = $contigFile;
    $contigQFile =~ s/.seq$//; $contigQFile .= ".qual";
    open(OUTFILE, ">$contigQFile");
    foreach my $contig (keys(%mergedContigs)) {
	$contigCtr++;
	@qs = ();
	@qsrev = ();
	for ($i = 0; $i <= length($contig)-$readLen; $i++){
	    $mer = substr($contig, $i, $readLen);
	    if (exists($rawReads{$mer})) {
		for ($j = 0; $j < $readLen; $j++) {
		    if ($qs[$i+$j] < $rawReads{$mer}[$j+1]) { $qs[$i+$j] = $rawReads{$mer}[$j+1]; }
		}
	    }
	    $revmer = revComp($mer);
	    if (exists($rawReads{$revmer})) {
		for ($j = 0; $j < $readLen; $j++) {
		    if ($qsrev[$i+$j] < $rawReads{$revmer}[$readLen-$j]) { $qsrev[$i+$j] = $rawReads{$revmer}[$readLen-$j]; }
		}
	    }
	}
	for ($i = 0; $i < length($contig); $i++) { $qs[$i] += $qsrev[$i]; }
	for ($i = 0; $i < length($contig); $i++) { if ($qs[$i] > 97) { $qs[$i] = 97; } }
	print(OUTFILE ">contig$contigCtr\n" . join(" ", @qs) . "\n");
    }
    close(OUTFILE);
}

# report summary on contigs
my $ncontigs = scalar(@lens);
@lens = sort {$b <=> $a} @lens;
$cumsum = 0;
my $meanlen = 0;
my $maxlen = 0;
my $n50 = 0;
my $n80 = 0;
foreach my $len (@lens) {
    if ($maxlen < $len) { $maxlen = $len; }
    $meanlen += $len;
    if ($cumsum < 0.5 * $nhat) { $n50 = $len; }
    if ($cumsum < 0.8 * $nhat) { $n80 = $len; }
    $cumsum += $len;
}
$meanlen = sprintf("%.2f", $meanlen/scalar(@lens));
my $coverage = sprintf("%.2f", 100*$cumsum/$nhat);
my $nhatf = sprintf("%6.1f", $nhat/1000);
$nhatf =~ s/^ *//;
print("\nAssembly statistics:\n",
    "----------------------------------------------------------------------\n",
    "Minimal contig length         : $minContigLen\n",
    "Number of contigs             : $ncontigs\n",
    "Length of largest contig      : $maxlen\n",
    "Average length of contigs     : $meanlen\n",
    "----------------------------------------------------------------------\n");

# write short contigs
my $shortsFile = $contigFile;
$shortsFile =~ s/[^\.]+$/shorts/;
open(OUTFILE, ">$shortsFile");
foreach my $contig (keys(%shortContigs)) {
    $contigCtr++;
    print(OUTFILE ">contig$contigCtr\n$contig\n");
}
close(OUTFILE);


########## Documentation

=head1 NAME

    SHARCGS - SHort-read Assembler based on Robust Contig-extension
    for Genomic Sequencing.

=head1 SYNOPSIS

    sharcgs.pl [options] <input-file>

=head1 LICENSE

      Using this software implies acceptance of the terms of the GNU General
      Public License (see http://sharcgs.molgen.mpg.de/docs/license.txt)

=head1 DESCRIPTION

    This program generates long contigs of genomic sequence based on very
    short 25-40mer, error prone reads as produced by 2nd generation sequencing
    machines. A large number of equal-sized reads is read from the input
    file and concatenated to generate contigs of several 1000 bases
    in length. The reads may contain errors in as many as 2% of all
    base calls.

    The algorithm filters the input for reads which fulfill the
    following criteria: they are read several times and their ends are
    matched perfectly by other, overlapping reads. When quality scores
    are provided for each call, the first step of the filtering is
    modified as follows: If a read is read several times, it is stored
    only once. For each call the highest observed quality score is
    retained. If two reads are reverse complements, their
    corresponding quality scores are added up. The criterion for the
    first filter step is then the smallest cumulated quality score for
    a read.

    The assembly algorithm then assumes that only correct reads pass
    this filter despite the unreliability of the raw data. The
    algorithm then only joins reads into contigs when this is possible
    unambiguously. In order to assemble contigs covering almost the
    complete target sequence, reads from a clear majority of positions
    must pass the mentioned filters.  When there is no read available
    from a position, this is called gap in this
    documentation. Actually, a gap can consist of several positions in
    a row, for which no read passes the filter.

    SHARCGS can perform several runs of the core assembly algorithm
    for different filter settings. When running the assembly for
    sloppy filtering, contig breaks are expected to occur frequently
    due to reads containing sequencing errors. With stringent
    filtering on the other hand, contig breaks occur because larger
    gaps must be spanned. The contig breaks, however, are not expected
    at the same positions. Therefore, SHARCGS gives the possibility to
    merge contigs from different core assembly runs.

    The final assembly is stored in a file containing one entry for
    each generated contig in FASTA format. If quality scores are
    provided as input, this output file is complemented by a
    FASTA-format file containing entries in the same order as the
    contig file. Instead of base calls, however, quality scores
    separated by spaces are reported. Quality scores are thereby
    computed by taking the maximum when several reads confirm a base
    call on the same strand. For each contig, quality scores
    corresponding to its reverse complement are also computed. The
    final quality scores are computed by adding up the corresponding
    quality scores from forward and reverse contig.

    The only mandatory argument is the input file, which must contain
    one read per line. Lines starting with a '>' sign are skipped from
    the input. All of these reads must have exactly the same
    length. Since we normally do not know, which strand has been read,
    all reads are converted to their reverse complement and added to
    the input data.

=head1 OPTIONS

    The following options can be used to fine-tune the assembly
    process by overriding the corresponding defaults given in
    parentheses:

=head2 -r <read-repetition-requirement> (default according to input data)

    Reads which are occuring in the input data less than this value
    are discarded from the assembly. Setting this value to 1 most
    probably integrates too many faulty reads into the assembly, since
    as many as 20% of all 30-nt reads contain sequencing errors, when
    the error rate per base call is 0.6%. Setting the read repetition
    requirement to values above 5 is not advisable for its wasteful
    consequences. This parameter can be set to a comma delimited list
    of values. In this case, one by one of these values is used as
    filter setting and the resulting contigs are merged into a single
    assembly. During the merge, the algorithm uses first generated
    contigs first and extends or merges them afterwards with contigs
    from subsequent runs.  It is therefore recommended to start the
    list with higher values, since contigs assembled with more
    stringent filtering are expected to be more reliable. If quality
    scores are provided by specifying a quality file, the values given
    to -r are interpreted as minimal quality scores.

    If the -r option is not specified SHARCGS chooses settings
    according to the input data. It estimates the length N of the
    target sequence from the observed read repetitions according to a
    maximum likelihood model. Then, SHARCGS chooses read repetition
    levels for which N/2 to N different reads pass the filter.

=head2 -g <gap-span> (default according to input data)

    This value sets the maximal length of a gap that can be safely
    spanned by the assembly algorithm. Setting this value too low
    leads to small contigs, because gaps that result from missing reads
    cannot be bridged. Moreover, SHARCGS cannot guarantee that
    contigs are joined correctly due to the lack of reads needed to
    detect ambiguous extensions. Setting this value too high causes
    the algorithm to detect spurious ambiguities and will also lead to
    short contigs due to resulting minimal overlaps beeing not
    unique. The best value for this parameter depends on the amount of
    repetitive sequence and the fraction of missing data.  It is
    actually translated into the minimal overlap we use to join a read
    to a contig.  Large gap-span values also imply a deeper prefix
    trie that needs more system memory. If the user sets this parameter
    to a specific value, the algorithm performs exactly one assembly
    using this particular gap span.  By default, the algorithm
    iterates through values between the estimated maximum gap size
    (see section OUTPUT) and the depth of the prefix trie to find the
    smallest gap span which leads to sufficiently few "open ends" (see
    below).

=head2 -q <quality-file> (default is none)

    In a quality file, quality scores for each base call in the main
    input file are specified. The read must be reported in the same
    order in both files and the quality file must conform to the
    Illumina 1G format.  For each base call, four integer values
    separated by space characters are reported. Each base call is
    separated by a tab and each read is reported on one line. Quality
    scores are interpreted as Q = -10 log10 (p/(1-p)) with p .


=head2 -d <max-dangling-ends> (default 2)

    Whenever SHARCGS stops the assembly of a contig, because no
    extension is found with a minimal overlap, we call this an open or
    dangling end. Open ends are expected at both end of a target
    sequence if the target is in linear form. If the target is
    circular, no open ends are expected, i.e. all contigs should be
    terminated due to detected ambiguities. This is the best guarantee
    that no wrong connections between contigs are generated. Thus for
    linear targets, this parameter should be set to 2, while for
    circular ones it should be set to 0.

=head2 -p <prefix-trie-depth> (default according to input data)

    All reads passing the filter are stored in a prefix trie. This
    data structure allows for very fast searching of extensions given
    a prefix. All reads are split into a prefix of minimal length and
    a suffix containing the rest of the read. The suffix is stored in
    a tree structure associated to the prefix using a hash table. In
    this manner, prefixes of any length larger than the minimal prefix
    length can be searched quickly and the corresponding suffix is
    determined easily. The prefix trie requires a large amount of
    memory. Setting this parameter to low values makes SHARCGS less
    greedy but limits it to span gaps shorter than this value.  If the
    parameter is not specified, it is set to a value such that large
    gaps can be spanned, as they are expected when reads from as many
    as 50% of target sequence positions are missing.

=head2 -l <min-contig-length> (default 50)

    The minimal length of contigs to be reported.

=head2 -t <end-trimming> (default 0)

    The number of nucleotides to be trimmed from the end of each read.

=head2 -o <output-file> (default contigs.seq)

    The output file where the contigs are to be stored in FASTA
    format.

=head1 OUTPUT

    The most important output, namely the collection of contigs
    assembled, is written to the output file specified with the '-o'
    parameter. In addition, the following messages generated on the
    screen may be of interest:

=head2 Read input data

    On this line the number of reads in the input file and the length
    of each individual read is given.

=head2 Add reverse complements

    In a first step of the assembly reverse complements of all reads
    are generated and stored into a single hash table. Here the number
    of unique reads from both strands is given.

=head2 Estimate sequence length/error rate

    The length of the target sequence and the error rate per base with
    which the reads are generated is estimated using a maximum
    likelihood method.  The information used is the distribution of
    how often reads are generated.  The probability that a read at a
    given position in the target sequence is read exactly x times, can
    be deduced from a binomial distribution Binom(p, k) using the
    number of reads k in the input file and the probability p that a
    read at a particular position occurs without error as
    parameters. Thereby, p depends on the error rate e, the read
    length l, and the length of the target sequence n:

=over 4

    p = (1 - e)^l / n

=back

    The number of reads which are expected to occur exacly x times in
    the input file is also deduced from a binomial distribution. This
    time the parameters are the length of the target sequence k and
    the probability to read x times without errors computed as
    mentioned above. Maximizing this probability over a grid of
    parameters deemed to be likely for a sequencing machine, we obtain
    the values given here. This method has shown to be accurate on
    simulated data but estimates are imprecise on real data.

=head2 Prefix trie depth

    This line reports, which fraction of missing reads can be handled
    and what depth of the prefix trie is needed to do so. If the -p or
    -g parameters are specified, the fraction of missing reads is
    adapted, otherwise the depth of the trie is adapted such that 50%
    missing reads can be handled.

=head2 Remove rare reads

    The first filter step removes reads which are not confirmed a
    minimal number of times or with a minimal quality score for any of
    its bases. This line states the number of reads that are kept
    after this filtering step.

=head2 Store frequent reads to prefix trie

    This line gives the number of suffix tree nodes allocated during
    the storage of all reads in the prefix trie.

=head2 Remove reads with no partners

    Here SHARCGS reports the number of reads kept after the removal of
    frequent reads which have no partner.

=head2 Guess read fraction/best gap span

    The read fraction given here is the estimated fraction of
    positions from the target sequence, from which reads survived the
    filtering steps. It is deduced from the estimated target sequence
    length and the number of reads kept after both filtering
    steps. The initial gap size to be spanned is chosen such that we
    expect only one gap of that length or longer. The probability that
    a gap of length g starts at a given position is f(1-f)^g, provided
    the read fraction f. This probability is multiplied by the length
    of the target sequence n to obtain the expected number of such
    gaps. Setting this value to 1 gives the first guess for a
    reasonable gap size for our algorithm:

=over 4

    n f (1 - f)^g = 1

    log(n) + log(f) + g log(1 - f) = 0

    g = -(log(n) + log(f)) / log(1 - f)

=back

    This gap size is used to initiate the assembly and increased if
    needed due to open ends.

=head2 Assemble contigs

    For each run of the core assembly algorithm one progress bar is
    shown. After an assembly is finished, the algorithm reports the
    number of contigs generated and the number of open ends
    detected. If this second number is larger than the corresponding
    parameter, the gap span is increased and another run of the core
    assembly algorithm follows.

=head2 Assembly statistics

    Finally, a summary of the generated assembly is given. It provides
    simple characteristics like the number of contigs generated, the
    length of the largest contig found, and the average length of
    contigs. 

=head1 AUTHORS

=over 4

=item Juliane Dohm <dohm@molgen.mpg.de>

=item Claudio Lottaz <Claudio.Lottaz@klinik.uni-regensburg.de>

=back 

=head1 REFERENCE

    Dohm JC, Lottaz C, Borodina T, Himmelbauer H. SHARCGS, a fast and
    highly accurate short-read assembly algorithm for de novo genomic
    sequencing. Genome Research, 2007, 17:1697-1706. PMID 17908823.

=head1 URL

=over 4

=item http://sharcgs.molgen.mpg.de

=back

=head1 EXAMPLES

    # trust defaults
    sharcgs.pl my.reads -o my.seq

    # circular targets 
    sharcgs.pl -d 0 my.reads -o my.seq

    # try to limit memory usage
    sharcgs.pl -p 10 my.reads -o my.seq

    # focus on large contigs
    sharcgs.pl -l 500 my.reads -o my.seq

    # tell the algorithm exactly, what to do
    sharcgs.pl -r 2 -g 10 my.reads -o my.seq

    # tell the algorithm, what option settings to explore
    sharcgs.pl -r 2,3,5 my.reads -o my.seq

=cut

########## End of file
