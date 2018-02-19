#!/usr/bin/perl

my $usage="Usage: $0 [-hv] [-t genetic code] [-p PAUP-Commands] [-s seqNames] [-n sampling] [fastaFile]\n" .
  "  -s: use the sequences listed only in the given file (not implemented)\n" .
  "  -p: PAUP commands to calculate distances\n" .
  "      By default, 'set criterion=distance; dset distance=k2p' is assumed\n".
  "      Internal distance engine: JC, JCm, P\n" .
  "  -n: number of times sampling will be repeated (not implemented)\n" .
  "  -t: supply alternative genetic code table\n" .
  "  -r: index of the reference sequences\n" .
  "  -v: verbose output\n" .
  " STDIN is used as the input if no fastaFile is given\n";

my $sep = "\t";  # if you use tab in the sequence name, change this to
                 # other characters such as ","

my $EXTRACT_PAIR_DIST_EXE = "~/myProgram/extractPairDist.pl";
my $PAIRDIST_EXE = "~/myProgram//cntPairwiseDiffs.pl";

use IO::File;
use POSIX qw(tmpnam);
use Getopt::Std;
use Data::Dumper;
getopts('hn:p:s:t:r:v') || die "$usage\n";

if (defined($opt_h)) {
    die "$usage\n";
}

die "$usage\n" if (@ARGV > 1);

# initialize the hash %aminoAcid
my $refSeqIndex;
my %aminoAcid;
if (defined($opt_r)) {  # -r specify the order of the reference sequence
    $refSeqIndex = $opt_r;
}
if (defined($opt_t)) {  # -t genCodeTbl was given
    open (CODE_TBL, "<$opt_t") || die "ERROR: Can't open $opt_t\n";
    %aminoAcid = InitGenCode(*CODE_TBL);
    close (CODE_TBL);
} else {                # use the default standard table supplied at the end
    %aminoAcid = InitGenCode(*DATA);
}

if ($opt_p eq 'JC' || $opt_p eq 'jc' || $opt_p eq 'JCm' || $opt_p eq 'p' || $opt_p eq 'P') {
    $distEngine = 'internal';
    $opt_p = 'JC' if ($opt_p eq 'jc');
    $opt_p = 'P' if ($opt_p eq 'p');
} else {
    $distEngine = 'PAUP';
}

@ARGV = ('-') unless @ARGV; # take STDIN when no arg.
my $seqFile = shift @ARGV;
my @dat = ReadInFASTA($seqFile);
my $numSeq = @dat;

@dat = AdjustSeqLength(@dat);    # attach '-' for shorter sequences
@dat = RemoveGapOnlySites(@dat);

@allFourFoldSite = FourFoldDegenerateSite(@dat);


my @seqName = GetSeqName(@dat);
for my $i (0 .. $#seqName) {
    print ">$seqName[$i]\n";
    print "$allFourFoldSite[$i]\n";
}

exit (0);

my @bootDat = (defined($opt_a)) ? BootstrapProtein(@dat) : BootstrapPAUP(@dat);

## check for NA
my @naNum = CntNA(@bootDat);
my $problemCnt = 0;
for $i (0..$#naNum) {
    if ($naNum[$i] > 0) {
	if ($problemCnt == 0) {
	    warn "Following pairwise distance estimates (1st column) " .
		"becamse infinity in n bootstraps (2nd column). " .
		    "These infinite distances are ignored in calculation " .
			"of cov and var.\n";
	}
	my $tmp = $i + 1;
	warn "$tmp\t $bootDat[$i]\n";
	$problemCnt++;
    }
}

my @covMat = MkCovMat (@bootDat);
print join "\n", @covMat;
print "\n";

exit(0);

# takes an arg; name of a file from which data are read Then read in
# the data and make an array.  Each element of this array corresponds
# to a sequence, name tab data.
sub ReadInFASTA {
    my $infile = shift;
    my @line;
    my $i = -1;
    my @result = ();
    my @seqName = ();
    my @seqDat = ();

    open (INFILE, "<$infile") || die "Can't open $infile\n";

    while (<INFILE>) {
        chomp;
        if (/^>/) {  # name line in fasta format
            $i++;
            s/^>\s*//; s/^\s+//; s/\s+$//;
            $seqName[$i] = $_;
            $seqDat[$i] = "";
        } else {
            s/^\s+//; s/\s+$//;
	    s/\s+//g;                  # get rid of any spaces
            next if (/^$/);            # skip empty line
            s/[uU]/T/g;                  # change U to T
            $seqDat[$i] = $seqDat[$i] . uc($_);
        }

	# checking no occurence of internal separator $sep.
	die ("ERROR: \"$sep\" is an internal separator.  Line $. of " .
	     "the input FASTA file contains this charcter. Make sure this " . 
	     "separator character is not used in your data file or modify " .
	     "variable \$sep in this script to some other character.\n")
	    if (/$sep/);

    }
    close(INFILE);

    foreach my $i (0..$#seqName) {
	$result[$i] = $seqName[$i] . $sep . $seqDat[$i];
    }
    return (@result);
}

# Initialize the hashTbl, %aminoAcid{codon},
# by reading in the FH given as the argument, and return the hash table
sub InitGenCode {  # take a typeglob of FILEHANDLE as an argument, e.g. *INFILE
    local *FH = shift;
    my $type;
    my (@aa, @b1, @b2, @b3);
    my $i;
    my $codon;
    my %aminoAcid = ();

    while (<FH>) {
	chomp;
	s/\s+$//;

	if (/^\s*(.*)\s*=/) {  # extract the type of the line
	    $type = $1;
	} else {
	    next;
	}

	s/^\s*(.*)=\s*//;      # get rid of characters before "="

	if ($type =~ /AAs/) {
	    @aa = split (//);
	} elsif ($type =~ /Base1/) {
	    @b1 = split (//);
	} elsif ($type =~ /Base2/) {
	    @b2 = split (//);
	} elsif ($type =~ /Base3/) {
	    @b3 = split (//);
	} else {
	    next;
	}
    }

    if (@aa + @b1 + @b2 + @b3 != 64 * 4) {  # checking the length of arrays
	die "ERROR, Please check the genetic code table is well formatted\n";
    }

    # making a hash table, %aminoAcid, Note all upper cases are used
    for ($i = 0; $i < 64; $i++) {
	$codon = uc ($b1[$i] . $b2[$i] . $b3[$i] );
	$aminoAcid{$codon} = uc $aa[$i];
    }

    return (%aminoAcid);
}


# Each element of a returned array is a result of one replication.
# Within each element, the distances are tab delimited. 
# n * (n-1) /2 distances per replicate.
sub FourFoldDegeneratePWDist {
    my @data = @_;
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);

    my ($tmpOutFile, $tmpSeqFileName);
    my @resultDistTbl = ();
    my @dist = ();

    # getting tmpfilenames
    my $fh;
    do { $tmpSeqFileName = tmpnam() }
    until $fh = IO::File->new($tmpSeqFileName, O_RDWR|O_CREAT|O_EXCL);
    close $fh;
    do { $tmpOutFile = tmpnam() }
    until $fh = IO::File->new($tmpOutFile, O_RDWR|O_CREAT|O_EXCL);
    close $fh;

    END { if (-e $tmpSeqFileName) {
	unlink($tmpSeqFileName) or die "Couldn't unlink $tmpSeqFileName : $!"}}
    END { if (-e $tmpOutFile) {
	unlink($tmpOutFile) or die "Couldn't unlink $tmpOutFile : $!"}}

    # prepare PAUP cmd
    if (defined($opt_p)) {
	$setting = $opt_p;
    } else {
	$setting = "set criterion=distance; dset distance=k2p";
    }
    my $paupCmd = "execute $tmpSeqFileName; $setting; " .
	"log start file=$tmpOutFile replace=yes ; showdist; " .
	    "log stop; quit WarnTSave=no;";
    warn "PAUP commands:\n$paupCmd\n" if ($distEngine eq 'PAUP');

    # opt_s isn't working now WORK HERE
    my $extractCmd = (defined ($opt_s)) ? 
	"$EXTRACT_PAIR_DIST_EXE -i -s $opt_s $tmpOutFile" :
	    "$EXTRACT_PAIR_DIST_EXE -i $tmpOutFile" ;

    # start looping through all pairs
    my ($i, $j);
    if(defined($refSeqIndex)){
	for my $j (0..$#data) {
	    next if $j == $refSeqIndex;
	    ## get 4-fold sites
	    my @sampledDat = ($data[$i], $data[$j]);
	    my @degenSites = FourFoldDegenerateSitePositions(@sampledDat);
	    print Dumper(@degenSites);
	    if (@degenSites < 1) { # no degenerate sites found
		push @resultDistTbl, "No4Degen";
		next;
	    }
	    @degenSites = map { $_ - 1 } @degenSites; # convert to 0-offset
	    @sampledDat = SubsetSites (\@sampledDat, \@degenSites);
	    push @resultDistTbl, join($sep, scalar(@degenSites), @dist);
	}
    }
    for my $i (1..$#data) {
	for my $j (0..($i-1)) {
	    ## get 4-fold sites
	    my @sampledDat = ($data[$i], $data[$j]);
	    my @degenSites = FourFoldDegenerateSitePositions(@sampledDat);
	    print Dumper(@degenSites);
	    if (@degenSites < 1) { # no degenerate sites found
		push @resultDistTbl, "No4Degen";
		next;
	    }
	    @degenSites = map { $_ - 1 } @degenSites; # convert to 0-offset
	    @sampledDat = SubsetSites (\@sampledDat, \@degenSites);

	    if ($distEngine eq 'PAUP') {  # check global var
		@dist = CalcDistWithPAUP ($paupCmd, $tmpSeqFileName, 
					  $extractCmd, @sampledDat);
	    } else { # internal
		@dist = CalcDistInternal (@sampledDat);
	    }
	    if (@dist != 1) {
		warn "WARN: For, $i and $j, number of distance returned by " .
		    "PAUP is not 1\n";
	    }

	    push @resultDistTbl, join($sep, scalar(@degenSites), @dist);
	}
    }
#    if (@dist != @seqName * (@seqName - 1) / 2) {
#	warn "## DANGER ## PAUP didn't out put correct number of " .
#	    "pairwise dists\n"
#    }
    return (@resultDistTbl);
}

sub FourFoldDegenerateSite{
    my @data = @_;
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);

    my @resultSeq = ();

    my @degenSites;
    # start looping through all pairs
    if(defined($refSeqIndex)){
	for my $j (0..$#data) {
	    next if $j == $refSeqIndex;
	    ## get 4-fold sites
	    my @sampledDat = ($data[$i], $data[$j]);
	    @degenSites = FourFoldDegenerateSitePositions(@sampledDat);
	    last;
	  }
    }
    foreach my $seq (@seqDat){
	  my @tmpSeq = ();
	  print Dumper(@degenSites);
	  foreach my $pos (@degenSites){
		print substr($seq, $pos-1,1);
		push(@tmpSeq, substr($seq, $pos-1, 1));
	  }
	#--------------------------------------------------
	#   print Dumper(@tmpSeq);
	#-------------------------------------------------- 
	  push(@resultSeq, join("",@tmpSeq));
    }
    return (@resultSeq);
}

sub CalcDistWithPAUP {
    my ($paupCmd, $tmpSeqFileName, $extractCmd, @sampledDat) = @_;
    my @dist = ();

    # write NEXUS to feed to PAUP
    my @seqName = GetSeqName(@sampledDat);
    my @degenDat = GetSeqDat(@sampledDat);
    if (! ContainAllBasesQ(join("", @degenDat))) {
	@dist = ("NotAllBases");
	return @dist;
    }
    WriteNEXUS ($tmpSeqFileName, \@seqName, \@degenDat);
    
    # run PAUP
    if (defined($opt_v)) {
	open (PAUP, "|paup -n");
    } else {
	open (PAUP, "|paup -n > /dev/null 2>&1");
    }
    print PAUP $paupCmd;
    close(PAUP);
    
    open (GETDIST, "$extractCmd |");
    while(<GETDIST>) {
	my @line = split;
	if ($. == 1) {
	    if ($line[$#line] ne "dist") {
		warn "## WARN: using the last column ($line[$#line]) ".
		    "of output from $EXTRACT_PAIR_DIST_EXE as " .
			"the distance\n"; 
	    }
	    next;
	}
	push @dist, $line[$#line];
    }
    close(GETDIST);

    return (@dist);
}

sub CalcDistInternal  {
    my @dat = @_;
    my @seqDat = GetSeqDat (@_);

    my @result = ();
    for my $i (1..$#seqDat) {
	for my $j (0..($i-1)) {
	    my $d = PWDist(($seqDat[$i], $seqDat[$j])) . "\t" . 
		VarPWDist(($seqDat[$i], $seqDat[$j]));
	    push @result, $d;
	}
    }
    return (@result);
}

# Two sequences are given as an array (no sequence names included).
sub PWDist {
    my @seqs = @_;
    my %subst = BaseChangeTbl (@seqs);

    my $pdist = $subst{'diff'}/$subst{'noGap'};
    if ($opt_p eq 'P') { # P - distance
	return $pdist;
    } elsif ($opt_p eq 'JC') { # JC for DNA
	return (($pdist == 0) ? 0 :
		(($pdist < 0.75) ? (- log(1 - 4 * $pdist / 3) * 3/4) : 'NA'));
    } elsif ($opt_p eq 'JCm') { # modified JC by Tajima
	return JCDistMod ($subst{'diff'}, $subst{'noGap'});
    }
}

# Two sequences are given as an array (no sequence names included).
sub VarPWDist {
    my @seqs = @_;
    my %subst = BaseChangeTbl (@seqs);

    my $pdist = $subst{'diff'}/$subst{'noGap'};
    if ($opt_p eq 'P') { # P - distance
	return $pdist * (1 -$pdist) / $subst{'noGap'} ;
    } elsif ($opt_p eq 'JC') { # JC for DNA
	return (($pdist < 0.75) ? 
	   (9 * $pdist * (1-$pdist)/((3 - 4 * $pdist)**2 * $subst{'noGap'})) : 
		'NA');
    } elsif ($opt_p eq 'JCm') {
	return VarJCDistMod ($subst{'diff'}, $subst{'noGap'});
    }
}

# Evolutionary distance with modified JC correction
# 1st arg is number of sites which differ, 2nd arg number of total sites.
# Tajima (1993) Mol Biol Evol 10:677-688
sub JCDistMod {
    my ($k, $n) = @_;
    my $b = 0.75;
    my $i;
    my $result = 0;
    for ($i = $k; $i >= 1; $i--) { # large to small $i to reduce round off err
	$result += ProductRatio ($k, $n, $i) / $i / $b ** ($i-1);
    }
    return $result;
}

sub MyPow {
    my ($a, $b) = @_;
    return ($b == 0) ? 1 : $a ** $b;
}

# returns n_i / d_i, where n_i is n (n-1) ... (n-i+1)
# all arguments should be integers
# Following function avoid integer overflow, but round off error become bigger
sub ProductRatio {
    my ($n, $d, $i) = @_;
    my $result = 1;
    my $j;
    for ($j = $i-1; $j >= 0; $j--) {
	$result *= ($n - $j)/($d - $j);
    }
    return $result;
}

# variance of evolutionary distance with modified JC correction
# 1st arg is number of sites which differ, 2nd arg number of total sites.
# Tajima (1993) Mol Biol Evol 10:677-688
sub VarJCDistMod {
    my ($k, $n) = @_;
    my $b = 0.75;
    my $dist = JCDistMod($k, $n);
    my $pdist = $k / $n;

    return $pdist * (1-$pdist) * exp(2 * $dist / $b) / ($n-1);
}


# Two sequences are given as an array (no sequence names included).
# returns an hash table, whose key is 'AC', 'CA', 'GG' etc.
# The value of 'AC' is the number of sites where seq1 is A and seq2 is C.
# Additionally, 4 special keys:
#   total: total number of sites (from longer sequences)
#   noGap: total number of sites without any gaps
#   gap1:  number of sites with 1 sequences having a gap
#   gap2:  number of sites with both sequences having a gap
#   diff:
#   same:
#   transition:
#   transversion:
# Ambiguous characters are considered as a gap.
sub BaseChangeTbl {
    my @seqs = @_;
    my %chTbl = ();

    # init the hash table
    my @keys =  qw(total noGap gap1 gap2 diff same transition transversion);
    @chTbl{@keys} = (0) x @keys;
    for my $i ((A,T,G,C)) {
	for my $j((A,T,G,C)) {
	    $chTbl{$i . $j} = 0;
	}
    }

    $seqs[0] =~ s/\?/-/g;  # converting unknow characters, ?, to -
    $seqs[1] =~ s/\?/-/g;

    my @seq1 = split //, $seqs[0];
    my @seq2 = split //, $seqs[1];

    my $len = Min(scalar(@seq1), scalar(@seq2));
    my $maxLen = Max(scalar(@seq1), scalar(@seq2));
    $chTbl{'total'} = $maxLen;

    if ($len != $maxLen) { # deal with unequal lengths
	my @extra = (@seq1 == $maxLen) ? splice(@seq1,$len) : 
	    splice(@seq2, $len);
	foreach my $b (@extra) {
	    if ($b eq '-') { $chTbl{'gap2'} ++; } else {$chTbl{'gap1'} ++;};
	}
    }

    for my $i (0..($len-1)) {
	my $b1 = $seq1[$i];
	my $b2 = $seq2[$i];
	if ($b1 ne '-' && $b2 ne '-') {
	    $chTbl{ $b1 . $b2 } ++;
	    $chTbl{'noGap'} ++;
	    if ($b1 ne $b2) {
		$chTbl{'diff'}++;
	    } else {
		$chTbl{'same'}++;
	    }
	} elsif ($b1 eq '-' && $b2 eq '-') {
	    $chTbl{'gap2'} ++;
	} else {
	    $chTbl{'gap1'} ++;
	}
    }

    $chTbl{'transition'} = $chTbl{'CT'} + $chTbl{'TC'} + 
	$chTbl{'AG'} + $chTbl{'GA'};
    $chTbl{'transversion'} = $chTbl{'diff'} - $chTbl{'transition'};

    return %chTbl;
}

# Each element of a returned array is a result of one replication.
# Within each element, the distances are tab delimited. 
# n * (n-1) /2 distances per replicate.
sub BootstrapPAUP {
    my @data = @_;
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);

    my $maxLen = MaxSeqLen(@data);
    my ($tmpOutFile, $tmpSeqFileName);
    my @resultDistTbl = ();

    # getting tmpfilenames
    my $fh;
    do { $tmpSeqFileName = tmpnam() }
    until $fh = IO::File->new($tmpSeqFileName, O_RDWR|O_CREAT|O_EXCL);
    close $fh;
    do { $tmpOutFile = tmpnam() }
    until $fh = IO::File->new($tmpOutFile, O_RDWR|O_CREAT|O_EXCL);
    close $fh;

    END { if (-e $tmpSeqFileName) {
	unlink($tmpSeqFileName) or die "Couldn't unlink $tmpSeqFileName : $!"}}
    END { if (-e $tmpOutFile) {
	unlink($tmpOutFile) or die "Couldn't unlink $tmpOutFile : $!"}}

    # prepare PAUP cmd
    if (defined($opt_p)) {
	$setting = $opt_p;
    } else {
	$setting = "set criterion=distance; dset distance=k2p";
    }
    my $paupCmd = "execute $tmpSeqFileName; $setting; " .
	"log start file=$tmpOutFile replace=yes ; showdist; " .
	    "log stop; quit WarnTSave=no;";
    warn "PAUP commands:\n$paupCmd\n";

    my $max = defined($opt_n) ? $opt_n : 1000;
    my $extractCmd = (defined ($opt_s)) ? 
	"$EXTRACT_PAIR_DIST_EXE -i -s $opt_s $tmpOutFile" :
	    "$EXTRACT_PAIR_DIST_EXE -i $tmpOutFile" ;

    # start sampling loop
    for my $repeat (0..($max-1)) {
	# make a sampled Data file
	my @sampledDat = SampleSites($maxLen, @seqDat);
	WriteNEXUS ($tmpSeqFileName, \@seqName, \@sampledDat);
	
	# run PAUP
	if (defined($opt_v)) {
	    open (PAUP, "|paup -n");
	} else {
	    open (PAUP, "|paup -n > /dev/null 2>&1");
	}
	print PAUP $paupCmd;
	close(PAUP);
	
	open (GETDIST, "$extractCmd |");
	my @dist = ();
	while(<GETDIST>) {
	    my @line = split;
	    if ($. == 1) {
		if ($line[$#line] ne "dist") {
		    warn "## WARN ## using the last column ($line[$#line]) " .
			"of output from $EXTRACT_PAIR_DIST_EXE as " .
			    "the distance\n"; 
		}
		next;
	    }
	    push @dist, $line[$#line];
	}
	close(GETDIST);
	push @resultDistTbl, join($sep, @dist);
    }

#    if (@dist != @seqName * (@seqName - 1) / 2) {
#	warn "## DANGER ## PAUP didn't out put correct number of " .
#	    "pairwise dists\n"
#    }
    return (@resultDistTbl);
}

# Each element of a returned array is a result of one replication.
# Within each element, the distances are tab delimited. 
# n * (n-1) /2 distances per replicate.
sub BootstrapProtein {
    my @data = @_;
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);

    my $maxLen = MaxSeqLen(@data);
    my $tmpSeqFileName;
    my @resultDistTbl = ();

    my $max = defined($opt_n) ? $opt_n : 1000;

    # getting tmpfilenames
    my $fh;
    do { $tmpSeqFileName = tmpnam() }
    until $fh = IO::File->new($tmpSeqFileName, O_RDWR|O_CREAT|O_EXCL);
    close $fh;

    END { if (-e $tmpSeqFileName) {
	unlink($tmpSeqFileName) or die "Couldn't unlink $tmpSeqFileName : $!"}}

    if ($opt_a ne "p" && $opt_a ne "poissonProt" && $opt_a ne "kimuraProt" ) {
	die "ERROR: -a takes poissonProt or kimuraProt\n";
    }

    my $distCmd = (defined ($opt_s)) ? 
	"$PAIRDIST_EXE -n -d $opt_a -s $opt_s $tmpSeqFileName" :
	    "$PAIRDIST_EXE -n -d $opt_a $tmpSeqFileName" ;

    # start sampling loop
    for my $repeat (0..($max-1)) {
	# make a sampled Data file
	my @sampledDat = SampleSites($maxLen, @seqDat);
	
	WriteFASTA ($tmpSeqFileName, \@seqName, \@sampledDat);

	# calculate pairwise distance
	open (PAIRDIST, "$distCmd|");
	my @dist=(); my $line = 1;
	while (<PAIRDIST>) {
	    chomp;
	    s/\s+//g;
	    next if (/^$/);
	    if (/Inf/) {
		$_ = "NA";
		warn "In $repeat-th sampling, infinite distance encountered " .
		    "at $line-th row\n";
	    }
	    push @dist, $_;
	    $line++;
	}
	close(PAIRDIST);
	# check for infitity
	push @resultDistTbl, join($sep, @dist);
    }

#    if (@dist != @seqName * (@seqName - 1) / 2) {
#	warn "## DANGER ## PAUP didn't out put correct number of " .
#	    "pairwise dists\n"
#    }
    return (@resultDistTbl);
}

sub CntNA {
    my @bootDat = TransposeMat(@_); # each row contains x bootstrapped 
                                    # pairwise dist for a particular pair
    my @result =();
    for my $i (@bootDat) {
	my @line = split /$sep/, $i;
	my @na = grep /NA/, @line;
	push @result, scalar @na;
    }
    return @result;
}

sub MkCovMat {
    my @bootDat = TransposeMat(@_); # each row contains x bootstrapped 
                                    # pairwise dist for a particular pair
    my @rowMean = RowMeans(@bootDat);
    
    my $numRow = @bootDat;
    my $numCol = -1;
    my @tmpMat = ();

    # make 2-dimensonal matrix
    for my $row (0..($numRow-1)) {
	my @thisRow = split(/$sep/, $bootDat[$row]);
	if ($numCol < 0) {
	    $numCol = @thisRow;
	} elsif (@thisRow != $numCol) {
	    die "ERROR: the number of columns are not equal in MkCovMat\n";
	}
	push @tmpMat, [ @thisRow ]
    }

    # subtract the row means from each element
    for my $row (0..($numRow-1)) {
	for my $col (0..($numCol-1)) {
	    if ($tmpMat[$row][$col] ne "NA") {
		$tmpMat[$row][$col] -= $rowMean[$row];
	    }
	}
    }

    my @result = ();
    for my $i (0..($numRow-1)) {
	my @covRow = ();
	for my $j (0..($numRow-1)) {
	    my $ss = 0;
	    my $cntr = 0;
	    for my $rep (0..($numCol-1)) {
		if ($tmpMat[$i][$rep] ne "NA" &&  $tmpMat[$j][$rep] ne "NA") {
		    $ss += $tmpMat[$i][$rep] * $tmpMat[$j][$rep];
		    $cntr++;
		}
	    }
	    push @covRow, $ss / ($cntr-1);
	}
	push @result, join("$sep", @covRow);
    }
    return @result;
}

sub TransposeMat {
    my @mat = @_;
    my @tmpMat =();
    my $numRow = @mat;
    my $numCol = -1;

    foreach my $row (0..($numRow-1)) {
	my @thisRow = split(/$sep/, $mat[$row]);
	if ($numCol < 0) {
	    $numCol = @thisRow;
	} elsif (@thisRow != $numCol) {
	    die "WARN: the number of columns are not equal in TransposeMat\n";
	}
	push @tmpMat, [ @thisRow ]
    }

    my @result = ();
    for my $i (0..($numCol-1)) {
	my @thisCol = ();
	for my $j (0..($numRow-1)) {
	    push @thisCol, $tmpMat[$j][$i];
	}
	push @result, join "$sep", @thisCol;
    }
    return @result;
}

sub RowMeans {
    my @array = @_;
    my @meanArray = ();

    foreach my $row (@array) {
	my @thisRow = split(/$sep/, $row);
	my $sum = 0;
	my $num = 0;
	foreach my $i (@thisRow) {
	    next if ($i eq "NA");
	    $sum += $i;
	    $num++;
	}
	push @meanArray, $sum / $num;
    }
    return @meanArray;
}



# note this function take only @seqDat 
sub SampleSites {
    my $maxLen = shift;
    my @seqDat = @_;
    my @result = ();

    my @randSites = RandIntArray($maxLen, $maxLen-1);

    for my $seqNumber (0 .. $#seqDat) {
	my @line = split (//, $seqDat[$seqNumber]);
	@line = SelectSites (\@line, \@randSites);
	my $randomized = join ("", @line);
	push @result, $randomized;
    }
    return (@result);
}

# rand integers between 0 and $max (both ends inclusive)
sub RandIntArray {
    my ($size, $max) = @_;
    my @result = ();

    for my $i (0 .. $size - 1) {
	push @result, int(rand ($max + 1));  # rand returns [0, $max + 1)
    }
    return (@result);
}

sub WriteNEXUS {
    my ($fileName, $nameArrayRef, $datArrayRef) = @_;

    my @nameArray = @$nameArrayRef;
    my @datArray = @$datArrayRef;

    die "Error in WriteNEXUS\n" if (@nameArray != @datArray);
    my $numSeq = @nameArray;
    my $seqLen = CharLen($datArray[0]);

    my $type = "nucleotide";
    if (defined ($opt_a)) {
	$type = "protein";
    }

    open (FP, ">$fileName") || die "Can't open a tmpFile $fileName";

    print FP "#NEXUS\nBegin data;\n" .
	"    Dimensions ntax=$numSeq nchar=$seqLen;\n" .
	    "    Format datatype=$type gap=- missing=? matchchar=.;\n" .
		"    Matrix\n";

    for my $i (0 .. $numSeq - 1) {
	print FP "\'$nameArray[$i]\' $datArray[$i]\n";
    }
    print FP "    ;\nEnd;\n";

    close(FP);
    return (0);
}

sub WriteFASTA {
    my ($fileName, $nameArrayRef, $datArrayRef) = @_;

    my @nameArray = @$nameArrayRef;
    my @datArray = @$datArrayRef;

    die "Error in WriteNEXUS\n" if (@nameArray != @datArray);
    my $numSeq = @nameArray;
    my $seqLen = CharLen($datArray[0]);

    open (FP, ">$fileName") || die "Can't open a tmpFile $fileName";

    for my $i (0 .. $numSeq - 1) {
	print FP ">$nameArray[$i]\n$datArray[$i]\n";
    }
    close(FP);
    return (0);
}

sub GetSeqDat {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
	@line = split (/$sep/, $i);
	push @result, $line[1];
    }

    return (@result)
}

sub GetSeqName {
    my @data = @_;
    my @line;
    my @result = ();

    foreach my $i (@data) {
	@line = split (/$sep/, $i);
	push @result, $line[0];
    }
    return (@result)
}

# take a sequence data array as the argument
# Returns an array of nucleotide positions which are 4-fold degenerate.
#
# A nucleotide position is 4-fold degenerate if the translation of the codon 
# which contains this site is same for all sequences, and changes of
# nucleotide at this sites do not change the AA for all sequences.
# If any sequences contain an ambiguous amino acid (? or -) at a site, the site
# is not four fold degenerate.
# unit offset (i.e. 1 means the first nucleotide)
sub FourFoldDegenerateSitePositions {
    my @seqDat = GetSeqDat(@_);
    my @aaSeq = Translate(@_);
    @aaSeq = GetSeqDat(@aaSeq);

    die "ERROR: Translation failed\n" if (scalar(@seqDat) != scalar(@aaSeq));
    my $numRow = @seqDat;

#    print join "\n", @seqDat, "\n";
#    print join "\n", @aaSeq, "\n";

    # make 2-dimensonal matrix
    my @aaMat = ();
    my $numCol = CharLen($aaSeq[0]);
    for my $row (0..$#aaSeq) {
	my @thisRow = split(//, $aaSeq[$row]);
	if (@thisRow != $numCol) {
	    die "ERROR: the number of columns of AATbl are not equal in " .
		"FourFoldDegenerate Sites\n";
	}
	push @aaMat, [ @thisRow ];
    }

    my @codonMat = ();
    for my $row (0..$#seqDat) {
	my @thisRow = MkTripletArray($seqDat[$row]);

	if (@thisRow != $numCol) {
	    die "ERROR: the number of columns of codonTbl (" . 
		scalar(@thisRow) . ") should be $numCol in " .
		    "FourFoldDegenerate Sites in row $row\n";
	}
	push @codonMat, [ @thisRow ];
    }

    # go through each site
    my @degenerateSites = ();
    for my $col (0..($numCol-1)) {
	my $thisAA = $aaMat[0][$col];
	for my $row (0..($numRow-1)) {
	    my $tmpAA = $aaMat[$row][$col];
#	    print "$col:$row:  $codonMat[0][$col]:$codonMat[$row][$col]:  $thisAA:$tmpAA\n"; #debug

	    if ($tmpAA eq '-' || $tmpAA eq '?' || $tmpAA ne $thisAA) {
		$thisAA = 0 ;  # signal that this site isn't conserved
		last;
	    }
	}
	## THIS DOESN"T HANDLE ? in 4-fold degenerate site.
	next if ($thisAA eq 0);
	# when reached here the all amino acids are conserved at this site
	for my $codonPosi (1..3) {
	    my $fourDegenQ = 1;
	    for my $row (0..($numRow-1)) {
		$tmp = Degeneracy($codonMat[$row][$col], $codonPosi);
#		print "DEGEN: $tmp  "; #debug
		if (Degeneracy($codonMat[$row][$col], $codonPosi) != 4) {
		    $fourDegenQ = 0;
		    last;
		}
	    }
#	    print "$fourDegenQ\n"; #debug
	    push @degenerateSites, $col * 3 + $codonPosi if ($fourDegenQ == 1);
	}
    }
    return @degenerateSites;
}

# 2 args, codon (e.g. AGG), and position (1, 2, or 3)
# returns 1, 2, or 4 (or some other values with non-standard genetic code)
# assumes global %aminoAcid is already initialized
# If the nucleotide at the position is ambiguous '?', 4 is returned 
# when it is 4-fold degenerate, -1 is returned otherwise.
# -1 is also returned when the other positions contains umbiguity.
sub Degeneracy {
    my ($codon, $posi) = @_;
    my @base = split //, $codon;
    die "ERROR: 1st arg. of Degeneracy() should be codon\n" if (@base != 3);
    die "ERROR: 2nd arg. of Degeneracy() should be 1, 2, or 3, not \"$posi\"\n"
	if ($posi !~ /^[123]$/);

    $posi--;               # making it to 0, 1 or 2 (i.e. zero offset)

    my $ambiguous = 0;
    if ($base[$posi] eq '?') {     # deal with ambiguous character
	$base[$posi] = 'A';        # temporarily putting A
	$ambiguous = 1;
    }

    return -1 if (! defined($aminoAcid{$codon}));
	   # we can't determine if the other positions have ambiguity
    my $thisAA = $aminoAcid{$codon};

    my $cnt = 0;
    foreach my $b ('A', 'T', 'G', 'C') {
	$base[$posi] = $b;
	my $tmpCodon = join "", @base;
	$cnt++ if ($thisAA eq $aminoAcid{$tmpCodon});
    }
    
    # when the position is ambiguous, we can be sure only when the site
    # is four-fold degenerate.
    if ($ambiguous == 1 && $cnt != 4) {
	$cnt = -1; 
    }
    return $cnt;
}

sub MkTripletArray {
    my $seq = shift;
    $seq =~ s/\s+//g;
    $seq =~ s/(.{3})/$1 /g;
    $seq =~ s/\s+$//;
    my @result = split(/ /, $seq);
    my $lastCharLen = CharLen($result[$#result]);
    # adjusting the last codon with "?"
    $lastCharLen = 3 - $lastCharLen;
    while ($lastCharLen > 0) {
	$result[$#result] = $result[$#result] . "?";
	$lastCharLen--;
    }
    return @result;
}


# assumes global %aminoAcid is already initialized
# returns an array of amino acid sequences (no name attached).
sub Translate {
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);

    my @result = ();
    for my $seqNum (0..$#seqDat) {
	my $seq = $seqDat[$seqNum];
	$seq = uc($seq);
	$seq =~ s/U/T/g;                 # converting RNA to DNA
	my $tmpAA = $seqName[$seqNum] . $sep;
	my @line = split //, $seq;
	my $seqLen = @line;
	for my $i (0 .. int($seqLen/3 - 1)) {
	    my $thisCodon = join ('', splice @line, 0, 3);
	    if ($thisCodon =~ /^[ATGC]{3}$/) {
		$tmpAA = $tmpAA . $aminoAcid{$thisCodon};
	    } elsif ($thisCodon =~ /^-{3}$/) { # contains all gaps
		$tmpAA = $tmpAA . "-";
	    } else {  # contains some gaps
		$tmpAA = $tmpAA . "?";
	    }
	}
	if (($seqLen % 3) != 0) {
		$tmpAA = $tmpAA . "?";	    
	}
	push @result, $tmpAA;
    }
    return @result;
}

# 1st arg is a REF to std. Seq data (name\tseq)
# 2nd arg is a REF to index of sites to be selected
# Returns the std seq array
sub SubsetSites {
    my ($arrayRef, $indexRef) = @_;
    unless (@_ == 2 && ref($arrayRef) eq 'ARRAY' && ref($indexRef) eq 'ARRAY'){
	die "args to SubsetSites() should be ARRAY REF, ARRAY REF\n";
    }
    my @seqName = GetSeqName(@$arrayRef);
    my @seqDat = GetSeqDat(@$arrayRef);
    my @result = ();

    for my $seqNumber (0 .. $#seqDat) {
	my @line = split (//, $seqDat[$seqNumber]);
	@line = SelectSites (\@line, $indexRef);
	my $selected = join ("", @line);
	push @result, "$seqName[$seqNumber]\t$selected";
    }
    return (@result);
}

sub SelectSites {
    my ($arrayRef, $indexRef) = @_;
    unless (@_ == 2 && ref($arrayRef) eq 'ARRAY' && ref($indexRef) eq 'ARRAY'){
	die "args to SelectSites() should be ARRAY REF, ARRAY REF\n";
    }

    my @result = ();
    foreach my $posi (@$indexRef) {
	push @result, $$arrayRef[$posi];
    }
    return @result;
}

sub MaxSeqLen {
    my @data = GetSeqDat(@_);
    my $maxLen = 0;
    foreach $i (@data) {
	my $len = CharLen($i);
	$maxLen = $len if ($len > $maxLen);
    }
    return ($maxLen);
}

sub ContainAllBasesQ {
    my $string = shift;
    if ($string =~ /A/  && $string =~ /T/ && $string =~ /G/ && $string =~ /C/){
	return 1;
    } else {
	return 0;
    }
}

# take std seq data (name\tseq), and attach "-" for the shorter sequences
sub AdjustSeqLength {
    my @data = @_;
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);
    my $maxLen = MaxSeqLen(@_);
    
    foreach $i (0 .. $#seqDat) {
	my $thisLen = CharLen ($seqDat[$i]);
	if ($thisLen == $maxLen)  {
	    ; # do nothing
	} elsif ($thisLen < $maxLen) {
	    my $diff = $maxLen - $thisLen;
	    warn "WARN: $seqName[$i] shorter.  " .
		"$diff '?' (missing character) were added at the end\n";
	    for ($j=0; $j < $diff; $j++) {
		$data[$i] = $data[$i] . "-";
	    }
	} else {
	    die "ERROR: the length of sequence $seqName[$i] is $thisLen, " .
		"longer than \$maxLen = $maxLen.  Weird!!";
	}
    }
    return (@data);
}

sub RemoveGapOnlySites {
    my @seqDat = GetSeqDat(@_);
    my @seqName = GetSeqName(@_);
    my $maxLen = MaxSeqLen(@_);
    my @gapSites = ();
    my @notGapSites = ();
    my ($posi, $seqNumber);
    my @seqMat = ();

    # make 2 dimensional matrix
    foreach $seqNumber (0..$#seqDat) {
	my @tmpArray = split(//, $seqDat[$seqNumber]);
	# Check the length
	if (@tmpArray != $maxLen)  {
	    die "ERROR: the sequence $seqName[$i] is not same length " .
		"as \$maxLen = $maxLen.  Weird!!";
	}
	push @seqMat, [ @tmpArray ];
    }

    # now identify the all gap sites
    for $posi (0 .. ($maxLen-1)) {
	my $gap = 1;
	for $seqNumber (0 .. $#seqMat){
	    if ($seqMat[$seqNumber][$posi] !~ /^[-\?]$/) {
		$gap = 0;
		last;
	    }
	}
	if ($gap == 1) {  # all sequences have a gap at these sites
	    push (@gapSites, $posi+1); # now unit-offset
	} else {          # there are some non-gap character at these sites
	    push (@notGapSites, $posi);
	}
    }

    # select sites and make results
    my @result = ();
    for $seqNumber (0 .. $#seqMat) {
	my @thisSeq = SelectSites($seqMat[$seqNumber], \@notGapSites);
	my $line = $seqName[$seqNumber] . $sep . (join("", @thisSeq));
	push (@result, $line);
    }

    if (@gapSites > 0) {
	warn ("Following sites consist of all gaps, removed from analysis\n");
	print STDERR join(" ", @gapSites);
	print STDERR "\n";
    }
    return (@result);
}

# count the number of characters in a string
sub CharLen {
    my $string = shift;
    my @charString = split (//, $string);
    return scalar(@charString);
}

# this function take two scalars and return the larger value
sub Larger {
    my ($a, $b) = @_;

    return (($a > $b) ? $a : $b);
}

sub Min {
    my $min = shift;
    for my $i (@_) {
	$min = $i if ($i < $min);
    }
    return $min;
}

sub Max {
    my $max = shift;
    for my $i (@_) {
	$max = $i if ($i > $max);
    }
    return $max;
}


sub sortByColumn {
# numerical sort by a column, return an array
#    sortbyColumn ($col_num, $order, @record)
# @record is an array with each element representing a space delimited record
# example
#    ("473 p1 S     0:06 -bash", "541 p2 SW    0:00 ps-a", ....)
# $col_num -- the column by which the record is sorted by (left-most col is 0)
# $order can be "a" (ascending) or "d" (descending),
# sort column can be hyphnated numbers (e.g. 10-4-150)

    local $col_num = shift(@_);
    local $order = shift(@_);
    local @record = @_ ;
    local ($sortCol);
    
    ## test if the sort column is hyphnated or plain number
    local $sortMethod = "number";
    foreach $sortCol (@record) {
	if ( (split(/\s+/,$sortCol))[$col_num] =~ /\d+-\d+/ ) {
	    $sortMethod = "hyphnated";
	    last ;
	}
    }

    return sort $sortMethod @record;

## two sub-sub-routines
    sub number {
	# $col_num, $order are the given arguments
	# the left-most column is 0 
	local $first = (split(/\s+/, $a))[$col_num];
	local $second = (split(/\s+/, $b))[$col_num];
# argument errors not well trapped here
	($first,$second) = ($second, $first) if ($order eq "d");
	
	return ($first <=> $second);
    }

#probably I don't need the "sub number"
    sub hyphnated {
	# $col_num, $order are the given arguments
	local ($each_number, $cmp_result, @temp_swap);

	## separte the hyphnated numbers and put them in the following arrays
        local @first = split(/-/, ((split(/\s+/, $a))[$col_num]));
	local @second = split(/-/, ((split(/\s+/, $b))[$col_num]));

	## ascending (default) or descending order
	if ($order eq "d") {
	    @temp_swap = @first;
	    @first = @second;
	    @second = @temp_swap;
	}
	
	## comparison of array elements
	for ($each_number = 0; $each_number <=
	     (($#first < $#second) ? $#first : $#second) ; $each_number++) {
	    $cmp_result = ($first[$each_number] <=> $second[$each_number]);
	    last if ($cmp_result);
	}

	## if the size of two arrays differ
	if ( ($cmp_result == 0) && ($#first != $#second) ) {
	    return (($#first < $#second) ? -1 : 1);
	} else {
	    return $cmp_result;
	}
    }
}

# The Standard Code (transl_table=1) from NCBI
# http://www3.ncbi.nlm.nih.gov/htbin-post/Taxonomy/wprintgc?mode=c#SG1
__DATA__
  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
