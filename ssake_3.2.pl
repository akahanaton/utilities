#!/usr/bin/perl 

#AUTHOR
#   Rene Warren (c) 2006-2008
#   Short Sequence Assembly by K-mer search and 3' Extension (SSAKE)
#   rwarren at bcgsc.ca

#NAME
#   SSAKE v3.2  Rene Warren, March 2008
#   Adjust contig ends to find new extension possibilities.  Important bug fix - addressed issues in failure to explore all read space (README for details)
#   Allows users to specify a seed sequence file.  This feature can be used to allow extension of existing DNA sequences using short reads
#   Supports mate pairs 
#   Handles errors in reads [implemented from a published idea - VCAKE v1.0  William Jeck, May 2007]

#SYNOPSIS
#   Progressive assembly of millions of short DNA sequences by k-mer search through a prefix tree and 3' extension

#DOCUMENTATION
#   SSAKE.readme distributed with this software @ www.bcgsc.ca
#   Warren RL, Sutton GG, Jones SJM, Holt RA.  2007.  Assembling millions of short DNA sequences using SSAKE.  Bioinformatics. 23(4):500-501
#   http://www.bcgsc.ca/platform/bioinfo/software/ssake
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca
#   If you use either the SSAKE code or ideas, please cite our work appropriately and accurately

#LICENSE
#   SSAKE Copyright (c) 2006-2008 Canada's Michael Smith Genome Science Centre.  All rights reserved.
#   Using a complete re-write of error-handling by consensus derivation (VCAKE) with its Copyright (c) 2007 University of North Carolina at Chapel Hill. All rights Reserved.

#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   note: insert size and distance between pairing reads are used interchangeably


use strict;
use Data::Dumper;
require "getopts.pl";
use vars qw($opt_f $opt_m $opt_o $opt_v $opt_r $opt_p $opt_d $opt_k $opt_a $opt_z $opt_e $opt_g $opt_s $opt_t $opt_b);
&Getopts('f:m:l:o:v:r:p:d:k:a:z:e:g:s:t:b:');
my ($base_overlap,$min_overlap,$verbose,$MIN_READ_LENGTH,$SEQ_SLIDE,$min_base_ratio,$paired,$insert_size,$min_links,$max_link_ratio,$contig_size_cutoff,$insert_stdev,$unpaired_file,$seed_file,$max_trim,$base_name)=(3,16,0,16,1,0.7,0,200,2,0.70,50,0.75,"no-g","no-s",0,"");
my $per;
my $MAX = 0; 
my $MAX_TOP = 100; # this is the very maximum anchoring edge sequence that will be searched

#-------------------------------------------------

if(! $opt_f){
   print "Usage: $0\n";
   print "-f  Fasta file containing all the [paired (-p 1) / unpaired (-p 0)] reads (required)\n";
   print "\t! paired reads must now be separated by \":\"\n";
   print "-s  Fasta file containing sequences to use as seeds exclusively (specify only if different from read set, optional)\n";
   print "-m  Minimum number of overlapping bases with the seed/contig during overhang consensus build up (default -m 16)\n";
   print "-o  Minimum number of reads needed to call a base during an extension (default -o 3)\n";
   print "-r  Minimum base ratio used to accept a overhang consensus base (default -r 0.7)\n";
   print "-t  Trim up to -t base(s) on the contig end when all possibilities have been exhausted for an extension (default -t 0)>\n";
   print "-p  Paired-end reads used? (-p 1=yes, -p 0=no, default -p 0)\n";
   print "-v  Runs in verbose mode (-v 1=yes, -v 0=no, default -v 0, optional)\n";
   print "-b  Base name for your output files (optional)\n";
   print "============ Options below only considered with -p 1 ============\n";
   print "-d  Mean distance expected/observed between paired-end reads (default -d 200, optional)\n";
   print "-e  Error (%) allowed on mean distance   e.g. -e 0.75  == distance +/- 75% (default -e 0.75, optional)\n";
   print "-k  Minimum number of links (read pairs) to compute scaffold (default -k 2, optional)\n";
   print "-a  Maximum link ratio between two best contig pairs *higher values lead to least accurate scaffolding* (default -a 0.70, optional)\n";
   print "-z  Minimum contig size to track paired-end reads (default -z 50, optional)\n";
   die "-g  Fasta file containing unpaired sequence reads (optional)\n";
}

my $file = $opt_f;
$min_overlap = $opt_m if ($opt_m);
$base_overlap = $opt_o if ($opt_o);
$min_base_ratio = $opt_r if ($opt_r);
$max_trim = $opt_t if ($opt_t);
$verbose = $opt_v if ($opt_v);
$paired = $opt_p if ($opt_p);
$insert_size = $opt_d if ($opt_d);
$min_links = $opt_k if ($opt_k);
$max_link_ratio = $opt_a if ($opt_a);
$contig_size_cutoff = $opt_z if ($opt_z);
$insert_stdev = $opt_e if ($opt_e);
$unpaired_file = $opt_g if ($opt_g);
$seed_file = $opt_s if($opt_s);
$base_name = $opt_b if($opt_b);

my $display_unpaired_file = $1 if ($unpaired_file=~/([^\/]*)$/);
my $display_seed_file = $1 if ($seed_file=~/([^\/]*)$/);
my $min_allowed = -1 * ($insert_stdev * $insert_size);
my ($low_iz, $up_iz) = ($insert_size + $min_allowed, $insert_size - $min_allowed);

#-------------------------------------------------

if(! -e $file){
   die "Invalid file: $file -- fatal\n";
}

if($opt_s && ! -e $opt_s){
   die "The file $opt_s you specified does not exist -- fatal\n";
}

if($paired && $insert_size <= 0){
   die "You specified that your input consist of paired-end reads (-p $paired), but have not specified a valid distance (-d $insert_size) -- fatal\n";
}

### Naming output files
if ($base_name eq ""){

   $base_name = $file . ".ssake_m" . $min_overlap . "_o" . $base_overlap . "_r" . $min_base_ratio . "_t" . $max_trim;

   if($paired){
      $base_name .= "_d" . $insert_size . "_e" . $insert_stdev . "_k" . $min_links . "_a" . $max_link_ratio . "_z" . $contig_size_cutoff . "_g-" . $display_unpaired_file;
   }
   if($opt_s){
      $base_name .= "_s-" . $display_seed_file;
   }

   my $pid_num = getpgrp(0);
   $base_name .= "_pid" . $pid_num;
}

my $contig = $base_name .  ".contigs";
my $singlet = $base_name . ".singlets";
my $short = $base_name . ".short";
my $log = $base_name . ".log";
my $scaffold = $base_name . ".scaffolds" if ($paired);
my $issues = $base_name . ".pairing_issues" if ($paired);
my $distribution = $base_name . ".pairing_distribution.csv" if ($paired);

open (LOG, ">$log") || die "Can't write to $log -- fatal\n";

if($min_overlap < 11 || $min_overlap > 50){
   my $outofbound_message = "-m must be a number between 11-50 ...Exiting.\n";
   print $outofbound_message;
   print LOG $outofbound_message;
   close LOG;
   exit;
}

if($base_overlap < 1){
   my $outofbound_message = "-o must be set to 1 or higher ...Exiting.\n";
   print $outofbound_message;
   print LOG $outofbound_message;
   close LOG;
   exit;
}

if($min_base_ratio <= 0.5 || $min_base_ratio > 1){
   my $outofbound_message = "-r must be a number between 0.51 and 1.00 ...Exiting.\n";
   print $outofbound_message;
   print LOG $outofbound_message;
   close LOG;
   exit;
}

#-------------------------------------------------

my $init_message = "\nRunning: $0\n-f $file\n-m $min_overlap\n-o $base_overlap\n-r $min_base_ratio\n-t $max_trim\n";
if($paired){$init_message .= "-p $paired\n-d $insert_size\n-e $insert_stdev\n-k $min_links\n-a $max_link_ratio\n-z $contig_size_cutoff\nUnpaired reads (optional) -g $unpaired_file\n\nScaffolds file: $scaffold\nPairing issues file: $issues\nPairing distance distribution file: $distribution\n";}
$init_message .= "Contigs file: $contig\nSinglets file: $singlet\nExcluded short reads file: $short\nLog file: $log\n";

print $init_message;
print LOG $init_message;

#-------------------------------------------------

my $date = `date`;
chomp($date);

my $reading_reads_message = "\nReading sequences initiated $date\n";
print $reading_reads_message;
print LOG $reading_reads_message;

my($set,$bin,$seed);
($set,$bin) = &readFasta($set,$bin,$file,$short,$paired);
($set,$bin) = &readFasta($set,$bin,$unpaired_file,$short,0) if (-e $opt_g);

#-------------------------------------------------
### Allow user to specify a fasta file containing sequences to use as seeds, exclusively
if(-e $opt_s){
   my $use_seed_message = "Using seed sequence file $opt_s for this assembly.\nNote: ONLY sequences in $opt_s will be used as seeds (i.e. -f $opt_f and -g $opt_g will NOT be used as seeds, only used for extension)\n";
   print LOG $use_seed_message;
   print $use_seed_message if ($verbose);

   $seed = &loadSeed($opt_s); 
}else{
   $seed = $set;
}
my $seed_number_message = "Number of unique seed sequences: " . keys( %$seed ) . "\n";
printf $seed_number_message;
print LOG $seed_number_message;
#-------------------------------------------------

$date = `date`;
chomp($date);

my $ssake_start_message = "\nSequence assembly initiated $date\n";
print $ssake_start_message;
print LOG $ssake_start_message;

#-------------------------------------------------
my ($sgl_count,$tig_count,$previous_index) = (1,1,0);

open (TIG, ">$contig") || die "Can't write to $contig -- fatal\n";
open (SIN, ">$singlet") || die "Can't write to $singlet -- fatal\n";
if ($paired){open (SC, ">$scaffold") || die "Can't write to $scaffold -- fatal\n";}

my ($tig_length,$track_all);

eval{

my $status_bar = "+";
for(my $i=1;$i<=99;$i++){
   $per->{$i}++;
   my $ct = $i /10;
   if($ct == int($ct)){$status_bar .= $ct;}else{$status_bar .= "-";}
}
$status_bar .= "+ x 10 [% complete]";
print "$status_bar\n.";

my $keys_start = keys ( %$seed );
#--------------------------------------------
ASSEMBLY:
foreach my $seq (sort {$seed->{$b}{'count'}<=>$seed->{$a}{'count'}} keys %$seed){#cycle through the input [normalized] reads

   my $track;

   if(defined $seed->{$seq}){#sequence read hasn't been used, is longer than 11 nt and the user-defined overlap minimum -m

      my $seed_name = "";
      if(defined $seed->{$seq}{'name'}){$seed_name = "|seed:" . $seed->{$seq}{'name'};}

      my $orig_mer = length($seq);

      if($paired){$track->{$seq}{'start'} = 1;$track->{$seq}{'end'} = $orig_mer;$track->{$seq}{'tig'} = $tig_count;}

      #### Delete keys ref
      my @o=split(//,$seq);                               

      my $start_sequence = $seq;
      my $reads_needed = $set->{$seq}{'count'};                      #tracks coverage
      my $total_bases = $orig_mer * $reads_needed;

      ($bin,$set,$seed)=deleteData($bin,$set,$seed,$seq);   #remove k-mer from hash table and prefix tree
     
      print "\n\n>>> START SEED SEQUENCE :: $seq <<<\n\n" if ($verbose);

      ($seq, $set, $bin, $reads_needed, $total_bases, $track) = doExtension("3", $orig_mer, $seq, $set, $bin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $track, $paired, $tig_count, $max_trim);
      ####end of 3' extension, beginning of 5' extension  (via 3' RC)

      my $seqrc = reverseComplement($seq);
      ($seqrc, $set, $bin, $reads_needed, $total_bases, $track) = doExtension("5", $orig_mer, $seqrc, $set, $bin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $track, $paired, $tig_count, $max_trim);
      ####end of 5' extension

      my $leng = length($seqrc);
      my $reversetig = reverseComplement($seqrc);                   ### return to sequence, as inputted
      my $trimmed_length = length($start_sequence) - 2*($max_trim);

      if($leng > $trimmed_length){ ### commented out: && $start_sequence ne $seqrc && $start_sequence ne $reversetig
         $tig_length->{$tig_count} = $leng;
         my $cov =  $total_bases / $leng;
         printf TIG ">contig%i|size%i|read%i|cov%.2f$seed_name\n%s\n", ($tig_count,$leng,$reads_needed,$cov,$reversetig);    #print contigs to file

         if($paired && $leng >= $contig_size_cutoff){
            $track_all = &trackReads($track,$track_all);
         }
         $tig_count++;
      }else{
         my $cov = $reads_needed;
         printf SIN ">singlet%i|size%i|read%i|cov%.2f$seed_name\n%s\n", ($sgl_count,$leng,$reads_needed,$cov,$start_sequence);    #print singlets to file
         $sgl_count++;
      }
   }

   my $keys_left = keys( %$seed );
   my $index = (int((($keys_start-$keys_left)/$keys_start)*100));
   if(defined $per->{$index}){
      print "." x ($index - $previous_index);
      $|=1; ###clear buffer
      delete $per->{$index};
   }
   $previous_index = $index;

   last ASSEMBLY if (! $keys_left);
}
print ".";
};###end eval block

$date = `date`;
chomp($date);

if($@){
   my $message = $@;
   my $failure = "\nSomething went wrong running $0 $date\n$message\n";
   print $failure;
   print LOG $failure; 
}else{
   my $success = "\nContig assembly executed normally $date\n";
   print $success;
   print LOG $success;
}

close TIG;
close SIN;
close SHO;

#------------------------------------
$date = `date`;
chomp($date);

if($paired){

   my $sc_start_message = "\nScaffolding initiated $date\n";
   print $sc_start_message;
   print LOG $sc_start_message;

   my $clone = &readClone($file);
   my $pair = &pairContigs($clone, $track_all, $tig_length, $issues, $distribution, $verbose);
   &buildScaffolds($pair, $tig_length, $contig_size_cutoff, $verbose);
 
   close SC;
   $date = `date`;
   chomp($date);

   my $sc_end_message = "\nScaffolding ended $date\n";
   print $sc_end_message;
   print LOG $sc_end_message;
}

close LOG;
exit;


#------------------------------------
#Order and orient contigs into scaffolds
sub buildScaffolds{

   my ($pair, $tig_length, $contig_size_cutoff, $verbose) = @_;

   my $seen_it;
   my $sc_ct = 0;
 
   #print SC "Scaffold Number,Scaffold Size (only contig lengths considered),Scaffold Chain: e.g. _f127z7068k12a0.58m42_r3090z62k7r0.14m76_  means: contig127(+ strand=f), size 7068 (z) has 12 links (k), link ratio of 0.58 (a) and with a mean gap/overlap of 42nt (m)  with reverse (r) of contig3090 (size 62) on the right.\n";

   SEED:
   foreach my $tig (sort {$tig_length->{$b}<=>$tig_length->{$a}} keys %$tig_length){

      my $ftig = "f" . $tig;
      my $rtig = "r" . $tig;

      if(! defined $seen_it->{$tig}){##should prevent re-using a contig as seed if it's already been incorporated into a scaffold

         $sc_ct++;

         my $chainleft = "";
          
         my $ori_chainright = $ftig . "Z" . $tig_length->{$tig};
         my $chainright = $ori_chainright;
         my $total = $tig_length->{$tig};

         ($total, $chainright, $seen_it) = &computeLayout("R", $chainright, $ftig, $pair, $tig_length, $total, $seen_it, $tig);
         ($total, $chainleft, $seen_it) = &computeLayout("L", $chainleft, $rtig, $pair, $tig_length, $total, $seen_it, $tig);

         $seen_it->{$tig}++;

         delete $pair->{$ftig};
         delete $pair->{$rtig};
         delete $tig_length->{$tig};

         my $scaffold = $chainleft . $chainright;
         print SC "scaffold" . $sc_ct . ",$total,$scaffold\n" if($total >= $contig_size_cutoff);
      }
   }
}

#------------------------------------
# links contigs together into a chain - must satisfy user-defined criterions (-k -a)
sub computeLayout{

   my ($ext, $chain, $tig, $pair, $tig_length, $total, $seen_it, $orig_tig_number) = @_;

   my $orig_tig = $tig;
   my $extension = 1;

   EXTENSION:
   while($extension){

      my $tnum = $1 if($tig=~/[fr](\d+)/);
      my $tnumf = "f" . $tnum;
      my $tnumr = "r" . $tnum;

      if(! defined $seen_it->{$tnum}){

         $seen_it->{$tnum}++ if($tnumf ne $orig_tig);

         print "Attempt to extend $tig\n" if ($verbose);      
         my $list = $pair->{$tig};
         my ($match1,$link1,$gaps1,$match2,$link2,$gaps2,$cntloop)=("",0,0,"",0,0,0);

         LINK:
         foreach my $match (sort {$list->{$b}{'links'}<=>$list->{$a}{'links'}} keys %$list){

            if($cntloop){
               ($match2,$link2,$gaps2) = ($match,$list->{$match}{'links'},$list->{$match}{'gaps'});
               print "$tig links second best $match2 (links:$link2 total sz:$gaps2)\n" if ($verbose);
               last LINK;
            }else{
               ($match1,$link1,$gaps1) = ($match,$list->{$match}{'links'},$list->{$match}{'gaps'});
               print "$tig links best $match1 (links:$link1 total sz:$gaps1)\n" if ($verbose);
            }
            $cntloop++;
         }

         ###ratio
         my $ratio = 0.00;
         $ratio = $link2 / $link1 if ($link1);        ## relative ratio of the two most abundant contig pairs
         if ($ratio =~ /(\d+\.\d{2})/){$ratio = $1;}
         ###mean
         my $mean = 0;
         $mean = int($gaps1 / $link1) if ($link1);

         my $tempnum = $1 if($match1 =~ /[fr](\d+)/);

         #### Assessment
         if(defined $seen_it->{$tempnum} || $link1 < $min_links || $ratio > $max_link_ratio || $tempnum == $orig_tig_number){
            $extension = 0;
            print "defined seen_it->{ $tempnum } || $link1 < $min_links || $ratio > $max_link_ratio\n L1:$link1 L2:$link2  M1:$match1 M2:$match2 G1:$gaps1 G2:$gaps2 "  if ($verbose);

            last EXTENSION;
         }{### pass filter.. does this contig 
            print "$ext extension.  mean: $mean links:$link1 linkratio:$ratio\n" if ($verbose);

            if($ext eq "R"){
               $chain .= "k" . $link1 . "a" . $ratio . "m" . $mean . "_" . $match1 . "z" . $tig_length->{$tempnum};
            }else{
               my $temp_match = "";
               if($match1 =~ /^r(\d+)/){$temp_match = "f" . $1;}else{$temp_match = "r". $1;}            
               $chain = $temp_match . "z" . $tig_length->{$tempnum} . "k" . $link1 . "a" . $ratio . "m" . $mean . "_" . $chain;
            }   
            $total += $tig_length->{$tempnum};

            print "NEXT TIG TO LOOK AT= $match1\n" if ($verbose);
            $tig = $match1;
            $extension = 1; 
          
            print "Will flag $tnum as seen  (only if $tnumf != $orig_tig)." if ($verbose);
   
            if($tnumf ne $orig_tig){
               delete $pair->{$tnumf};
               delete $pair->{$tnumr};
               delete $tig_length->{$tnum};
            }else{
               delete $pair->{$tnumf};
            }
         }
      }else{
         print "NO MORE MATCH FOR $tig in hash: pair>>\n" if ($verbose);
         $extension = 0;
         last EXTENSION;
      }
   }### pair is defined
   return $total, $chain, $seen_it;
}

#------------------------------------
sub trackReads{

   my ($track, $track_all) = @_;

   foreach my $rd (keys %$track){
      if(! defined $track_all->{$rd}{'tig'}){
         $track_all->{$rd}{'tig'} = $track->{$rd}{'tig'};
         $track_all->{$rd}{'start'} = $track->{$rd}{'start'};
         $track_all->{$rd}{'end'} = $track->{$rd}{'end'};
         delete $track->{$rd};
      }
   }
   return $track_all;
}

#------------------------------------
sub getDistance{

   my ($insert_size, $length_i, $start_i, $start_j) = @_;

   # L  ------  --------- R
   # i    ->        <-    j
   #      ....  ......    insert_span
   #      ============    insert_size

   my $insert_span = ($length_i - $start_i) + $start_j;
   my $gap_or_overlap = $insert_size - $insert_span;

   return $gap_or_overlap;
}

#-----------------
#build contig pairs based on template information  -  must satisfy user-defined criterions (-d -e)
sub pairContigs{

   my ($clone,$track,$tig_length,$issues,$distribution,$verbose) = @_;
   my ($ct_illogical, $ct_ok_contig, $ct_ok_pairs, $ct_problem_pairs, $ct_iz_issues, $ct_single, $ct_both)= (0,0,0,0,0,0,0);

   my ($pair,$err,$track_insert);

   print "Pairing contigs...\n" if ($verbose);

   open(PET, ">$issues") || die "Can't open $issues for writing -- fatal\n";

   foreach my $template (keys %$clone){ 

      my $read_a = $clone->{$template}{'a'};
      my $read_b = $clone->{$template}{'b'};

      print "Pair#$template read1=$read_a read2=$read_b\n" if ($verbose);

      if(defined $track->{$read_a}{'tig'} && $track->{$read_b}{'tig'}){### both pairs assembled

         $ct_both++;

         my $tig_a = $track->{$read_a}{'tig'};
         my $tig_b = $track->{$read_b}{'tig'};

         my $ftig_a = "f" . $tig_a;
         my $ftig_b = "f" . $tig_b;

         my $rtig_a = "r" . $tig_a;
         my $rtig_b = "r" . $tig_b;

         my $A_length = $tig_length->{$tig_a};
         my $A_start = $track->{$read_a}{'start'};
         my $A_end = $track->{$read_a}{'end'};

         my $B_length = $tig_length->{$tig_b};
         my $B_start = $track->{$read_b}{'start'} ;
         my $B_end = $track->{$read_b}{'end'};

         if ($tig_a != $tig_b){####paired reads located on <> contigs

            ####Determine most likely possibility
            if ($track->{$read_a}{'start'} < $track->{$read_a}{'end'}){

               if ($track->{$read_b}{'end'} < $track->{$read_b}{'start'}){####-> <- :::  A-> <-B  /  rB -> <- rA
                   my $d = &getDistance($insert_size, $A_length, $A_start, $B_start);
                   print "A-> <-B  WITH $tig_a -> <- $tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen, Astart,Bstart\n" if($verbose);
                   if($d >= $min_allowed){
                      $pair->{$ftig_a}{$ftig_b}{'links'}++;
                      $pair->{$ftig_a}{$ftig_b}{'gaps'} += $d;                  
                      $pair->{$rtig_b}{$rtig_a}{'links'}++;
                      $pair->{$rtig_b}{$rtig_a}{'gaps'} += $d;
                      $ct_ok_pairs++;
                   }else{
                      my $err_pair = $ftig_a . "-". $ftig_b;
                      $err->{$err_pair}{'links'}++;
                      $err->{$err_pair}{'gaps'} += $d;
                      $ct_problem_pairs++;
                      print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#$tig_a -> $d <- tig#$tig_b, A=$A_length nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                   }
                }else{#### -> -> ::: A-> <-rB  / B-> <-rA 
                   my $rB_start = $B_length - $B_start;
                   my $d = &getDistance($insert_size, $A_length, $A_start, $rB_start);
                   print "A-> <-rB  WITH $tig_a -> <- r.$tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen,Astart,rBstart\n" if($verbose);
                   if($d >= $min_allowed){
                      $pair->{$ftig_a}{$rtig_b}{'links'}++;
                      $pair->{$ftig_a}{$rtig_b}{'gaps'} += $d;
                      $pair->{$ftig_b}{$rtig_a}{'links'}++;
                      $pair->{$ftig_b}{$rtig_a}{'gaps'} += $d;
                      $ct_ok_pairs++;
                   }else{
                      my $err_pair = $ftig_a . "-". $rtig_b;
                      $err->{$err_pair}{'links'}++;
                      $err->{$err_pair}{'gaps'} += $d;
                      $ct_problem_pairs++;
                      print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-rB  WITH tig#$tig_a -> $d <- tig#r.$tig_b, A=$A_length  nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                   }
                }
            }else{

               if ($track->{$read_b}{'end'} > $track->{$read_b}{'start'}){####<-  -> ::: B-> <-A / rA -> <- rB
                  my $d = &getDistance($insert_size, $B_length, $B_start, $A_start);
                  print "B-> <-A  WITH $tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,Bstart,Astart\n" if($verbose);
                  if($d >= $min_allowed){
                     $pair->{$ftig_b}{$ftig_a}{'links'}++;
                     $pair->{$ftig_b}{$ftig_a}{'gaps'} += $d;
                     $pair->{$rtig_a}{$rtig_b}{'links'}++;
                     $pair->{$rtig_a}{$rtig_b}{'gaps'} += $d;
                     $ct_ok_pairs++;
                  }else{
                     my $err_pair = $ftig_b . "-". $ftig_a;
                     $err->{$err_pair}{'links'}++;
                     $err->{$err_pair}{'gaps'} += $d;
                     $ct_problem_pairs++;
                     print PET "Pairs unsatisfied in distance within a contig pair.  B-> <-A  WITH tig#$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                  }
               }else{                          ####<- <-  :::  rB-> <-A / rA-> <-B
                  my $rB_start = $B_length - $B_start;
                  my $d = &getDistance($insert_size, $B_length, $rB_start, $A_start);
                  print "rB-> <-A WITH r.$tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,rBstart,Astart\n" if($verbose);
                  if($d >= $min_allowed){
                     $pair->{$rtig_b}{$ftig_a}{'links'}++;
                     $pair->{$rtig_b}{$ftig_a}{'gaps'} += $d;
                     $pair->{$rtig_a}{$ftig_b}{'links'}++;
                     $pair->{$rtig_a}{$ftig_b}{'gaps'} += $d;
                     $ct_ok_pairs++;
                  }else{
                     my $err_pair = $rtig_b . "-". $ftig_a;
                     $err->{$err_pair}{'links'}++;
                     $err->{$err_pair}{'gaps'} += $d;
                     $ct_problem_pairs++;
                     print PET "Pairs unsatisfied in distance within a contig pair.  rB-> <-A WITH tig#r.$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
                  }
               }
            }
            #print Dumper($pair);
         }else{###Clone, paired reads located on the same contig -- could be used to investigate misassemblies
           
            print "Pair ($read_a and $read_b) located on same contig $tig_a ($A_length nt)\n" if ($verbose);
            my $pet_size = 0;

            if ($A_start > $B_start && ($B_start < $B_end) && ($A_start > $A_end)){    # B --> <-- A
               $pet_size = $A_start - $B_start;
               $track_insert->{$pet_size}++;
               if($pet_size >= $low_iz && $pet_size <= $up_iz){
                  $ct_ok_contig++;
               }else{
                  print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
                  $ct_iz_issues++;
               }
            }elsif($B_start > $A_start && ($B_start > $B_end) && ($A_start < $A_end)){ # A --> <-- B
               $pet_size = $B_start - $A_start;
               $track_insert->{$pet_size}++;
               if($pet_size >= $low_iz && $pet_size <= $up_iz){
                  $ct_ok_contig++;
               }else{
                  print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
                  $ct_iz_issues++;
               }
            }else{
               $ct_illogical++;
               print PET "Pairs unsatisfied in pairing logic within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end\n";
            }
         }
      }else{###both pairs assembled
         $ct_single++;
      }
   }###each template

   ### summary of contig pair issues
   print PET "------------- Putative issues with contig pairing - Summary  ----------------\n";
   foreach my $err_pair (sort {$err->{$b}{'links'}<=>$err->{$a}{'links'}} keys %$err){
      my $mean_iz = 0;
      $mean_iz = $err->{$err_pair}{'gaps'} / $err->{$err_pair}{'links'} if ($err->{$err_pair}{'links'});
      print PET "Pair $err_pair has $err->{$err_pair}{'links'} links and mean distance = $mean_iz\n";
   }
   close PET;
 
   my $satisfied = $ct_ok_pairs + $ct_ok_contig;
   my $unsatisfied = $ct_problem_pairs + $ct_iz_issues + $ct_illogical;
   my $ct_both_reads = $ct_both * 2;

   print LOG "\nPAIRED-END READS STATS\n";
   print LOG "-" x 20, "\n";
   print LOG "At least one sequence/pair missing from contigs >= $contig_size_cutoff bp (user-defined -z): $ct_single\n";
   print LOG "Assembled pairs: $ct_both ($ct_both_reads sequences)\n";
   print LOG "\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: $insert_size +/$min_allowed): $ct_ok_contig\n";
   print LOG "\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds): $ct_iz_issues\n";
   print LOG "\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): $ct_illogical\n";
   print LOG "\t---\n";
   print LOG "\tSatisfied in distance/logic within a given contig pair (pre-scaffold): $ct_ok_pairs\n";
   print LOG "\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): $ct_problem_pairs\n";
   print LOG "\t---\n";
   print LOG "Total satisfied: $satisfied\tunsatisfied: $unsatisfied\n\n";

   open (CSV, ">$distribution") || die "Can't open $distribution for writing -- fatal";

   foreach my $is (sort {$a<=>$b} keys %$track_insert){
      print CSV "$is,$track_insert->{$is}\n";
   }

   close CSV;
   return $pair;
}

#-----------------
# SSAKE contig extension
sub doExtension{

   my ($direction, $orig_mer, $seq, $set, $bin, $reads_needed, $total_bases, $min_overlap, $base_overlap, $min_base_ratio, $verbose, $track, $paired, $tig_count, $max_trim) = @_;

   my $previous = $seq;
   my $extended = 1;
   my $trim_ct = 0;     #trim counter - keeps track of 3'-end trim

   TRIM:
   while($trim_ct <= $max_trim){

      while($extended){

         ### added 31Jan08 R.Warren
         if($orig_mer > $MAX){$orig_mer=$MAX;}  ### Deals with special cases where the seed sequences are different from the read set (and possibly very large) - goal here is not to increase sequence coverage of seed, but rather to extend it.

         my ($ct,$pos,$current_reads,$current_bases,$span) = (0,0,0,0,$orig_mer);

         my $overhang = {};
         my @overlapping_reads = ();
         for (my $x=1;$x <= ($orig_mer *2);$x++){
            ($overhang->{$x}{'A'},$overhang->{$x}{'C'},$overhang->{$x}{'G'},$overhang->{$x}{'T'}) = (0,0,0,0);
         }

         EXTENSION:
         while ($span > $min_overlap){  #will slide the subseq, until the user-defined min overlap size

            ### Added 19March08
            if(length($seq) >= $MAX){   # $seq is length of contig being extended -- if larger than largest read, make sure the largest read could align and all subsequent rds. 
               $span = $MAX - $ct;
            }else{
               $span = length($seq) - $ct;
            }

            $pos = length($seq) - $span;

            my $subseq = substr($seq, $pos, $span);                                                            #make a sub-sequence of length l-(1..i) for searching

            my @s=split(//,$subseq);
            my $subset = $bin->{$s[0]}{$s[1]}{$s[2]}{$s[3]}{$s[4]}{$s[5]}{$s[6]}{$s[7]}{$s[8]}{$s[9]}{$s[10]}; #Will grab everything even the reverse complement ones

            print "####$direction' SEARCH Counter:$ct Position:$pos Span:$span - Subseq:$subseq Previous:$previous\n" if ($verbose);

            SEARCH:   #this cycles through limited k-mer space
            foreach my $pass (sort {$subset->{$b} <=> $subset->{$a}} keys %$subset){
               if($pass =~ /^$subseq([ACGT]*)/){ 
                  #can we align perfectly that subseq to another rd start?
                  my $dangle = $1;
                  print "\n", "=" x 80, "\n$direction'- FOUND sequence: $pass -> subset: $subseq -> overhang: $dangle\n", "=" x 80, "\n\n" if ($verbose);

                  # Collect all overhangs
                  push @overlapping_reads, $pass;                  ### all overlapping reads
                  my @over = split(//,$dangle);
                  my $ct_oh = 0;
 
                  foreach my $bz(@over){
                     $ct_oh++;                                     ### tracks overhang position passed the seed  
                     if(defined $set->{$pass}){
                        $overhang->{$ct_oh}{$bz} += $set->{$pass}{'count'};      ### reflects read coverage (often real duplicates)
                     }else{
                        my $pass_rc = reverseComplement($pass);
                        $overhang->{$ct_oh}{$bz} += $set->{$pass_rc}{'count'};
                     } 
                     print "$ct_oh - $bz = $overhang->{$ct_oh}{$bz}\n" if($verbose);
                  }
               }elsif($subseq =~ /$pass/){ ###cases where the read is fully embedded in the search sequence - want to include for coverage calculations
                  my $complement_pass = reverseComplement($pass);

                  print "$pass found in $subseq ($set->{$pass}{'count'}) - deleting read: $pass and complement ($set->{$complement_pass}): $complement_pass\n\n" if ($verbose);

                  if(defined $set->{$pass}){
                     $current_reads = $set->{$pass}{'count'};
                     $current_bases = length($pass) * $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$pass);
                     if($paired && ! defined $track->{$pass}){
                        $track->{$pass}{'tig'} = $tig_count;
                        $track->{$pass}{'start'} = $pos + 1;
                        $track->{$pass}{'end'} = $pos + length($pass);
                     }
                  }
                  if(defined $set->{$complement_pass}){
                     $current_reads = $set->{$complement_pass}{'count'};
                     $current_bases = length($complement_pass) * $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$complement_pass);
                     if($paired && ! defined $track->{$complement_pass}){ 
                        $track->{$complement_pass}{'tig'} = $tig_count;
                        $track->{$complement_pass}{'end'} = $pos + 1;
                        $track->{$complement_pass}{'start'} = $pos + length($complement_pass);
                     }
                  }
               }  
            }
            $ct += $SEQ_SLIDE;
         }#while overlap >= user-defined -m minimum
 
         my $consensus = "";
         print "Finished Collecting Overlapping Reads - BUILDING CONSENSUS...\n" if ($verbose);
         print Dumper(@overlapping_reads) if ($verbose);

         ### Build consensus
         CONSENSUS:
         foreach my $ohpos (sort {$a<=>$b} keys %$overhang){
            if($ohpos){

               my $coverage = $overhang->{$ohpos}{'A'}+$overhang->{$ohpos}{'C'}+$overhang->{$ohpos}{'G'}+$overhang->{$ohpos}{'T'};
               print "pos:$ohpos cov:$coverage A:$overhang->{$ohpos}{'A'} C:$overhang->{$ohpos}{'C'} G:$overhang->{$ohpos}{'G'} T:$overhang->{$ohpos}{'T'}\n" if($verbose);

               if ($coverage < $base_overlap){
                  print "COVERAGE BELOW THRESHOLD: $coverage < -o $base_overlap @ $ohpos :: will extend by: $consensus\n" if ($verbose);
                  last CONSENSUS;
               }
               my $baselist = $overhang->{$ohpos};

               my $ct_dna=0;
               my $previous_bz = "";

               BASE:
               foreach my $bz (sort {$baselist->{$b}<=>$baselist->{$a}} keys %$baselist){
                  #print "\t$ct_dna -> $bz..$baselist->{$previous_bz} > $baselist->{$bz}\n";
                  if($ct_dna){## the two most abundant bases at that position
                     #print "\t\t$ct_dna\n";
                     if($previous_bz ne "" && ($baselist->{$previous_bz} / $coverage) >= $min_base_ratio && $baselist->{$previous_bz} > $baselist->{$bz}){### a simple consensus btw top 2 
                        $consensus .= $previous_bz;                                         ### build consensus
                        print "Added base $previous_bz (cov = $baselist->{$previous_bz}) to $consensus **\n" if ($verbose);
                        last BASE;
                     }else{
                        print "ISSUES EXTENDING: best base = $previous_bz (cov=$baselist->{$previous_bz}) at $ohpos.  Second-Best: $bz (cov=$baselist->{$bz}) (ratio best=$baselist->{$previous_bz} / total=$coverage) >= $min_base_ratio (-r) -- will terminate with $consensus\n" if($verbose);
                        last CONSENSUS;
                     }
                  }
                  $previous_bz = $bz;                 
                  $ct_dna++;
               }
            }
         }

         ### deal with sequence reads making up the consensus/newly formed contig
         if($consensus ne ""){
            print "Will extend $seq\nwith: $consensus\n\n" if($verbose);
            my $temp_sequence = $seq . $consensus;  ## this is the contig extension
            my $integral = 0;
            foreach my $ro (@overlapping_reads){

               while($temp_sequence =~ /$ro/g){                                   ### read found integral in the newly built sequence
                  my $complement_ro = reverseComplement($ro);
                  $integral=1;

                  print "$ro found in $seq ($set->{$ro}{'count'}) - deleting read: $ro and complement ($set->{$complement_ro}{'count'}): $complement_ro\n\n" if ($verbose); 
                  if(defined $set->{$ro}){          
                     $current_reads = $set->{$ro}{'count'};  
                     #print "fwd SET:$current_reads BIN $subset->{$ro}\n";
                     $current_bases = length($ro) * $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$ro);
                     if($paired && ! defined $track->{$ro}){
                        $track->{$ro}{'tig'} = $tig_count;
                        $track->{$ro}{'start'} = (pos $temp_sequence) - length($ro) + 1;
                        $track->{$ro}{'end'} =  (pos $temp_sequence);
                     }
                  }

                  if(defined $set->{$complement_ro}){          
                     $current_reads = $set->{$complement_ro}{'count'}; 
                     #print "rc SET:$current_reads BIN $subset_rc->{$complement_ro}\n";
                     $current_bases = length($complement_ro) * $current_reads;
                     $reads_needed += $current_reads;
                     $total_bases += $current_bases;
                     ($bin,$set,$seed) = deleteData($bin,$set,$seed,$complement_ro);
                     if($paired && ! defined $track->{$complement_ro}){
                        $track->{$complement_ro}{'tig'} = $tig_count;
                        $track->{$complement_ro}{'end'} = (pos $temp_sequence) - length($complement_ro) + 1;
                        $track->{$complement_ro}{'start'} = (pos $temp_sequence);
                     }
                  }
               }
            }
            if(! $integral){### no reads are found overlapping with the consensus might be indicative of low complexity regions -- Stop the extension
               print "No overlapping reads agree with the consensus sequence.   Stopping extension" if ($verbose);
               $extended = 0;
            }else{
               $seq = $temp_sequence;
               print "New Contig is: $seq\n" if ($verbose);
               $extended = 1;
            }
            $previous = $seq;
         }else{### no consensus built, will stop the extension
            $extended = 0;
         }

      }###while get the OK for extension

      $trim_ct++;
      if ($trim_ct <= $max_trim){
         last TRIM if (length($seq) <= $MIN_READ_LENGTH); #terminate assembly if trimming becomes too agressive
         $seq = substr($seq, 0, -1);
         $extended = 1;
         print "\n$direction prime EXTENSION ROUND $trim_ct COMPLETE UNTIL $max_trim nt TRIMMED OFF => TRIMMED SEQUENCE:$seq\n\n" if ($verbose);
      }
      
   }### while trimming within bounds
   #### Adjust the position if tracking paired reads in assembly
   if($paired){
      foreach my $rd (keys %$track){
         $track->{$rd}{'start'} = length($seq) - $track->{$rd}{'start'};
         $track->{$rd}{'end'} = length($seq) - $track->{$rd}{'end'};
      }
   }
   
   print "\n*** NOTHING ELSE TO BE DONE IN $direction prime- PERHAPS YOU COULD DECREASE THE MINIMUM OVERLAP -m (currently set to -m $min_overlap) ***\n\n" if ($verbose);

   return $seq, $set, $bin, $reads_needed, $total_bases, $track;
}


#-----------------------
sub deleteData {
   my ($bin,$set,$seed,$sequence) = @_;
   
   my @o=split(//,$sequence);
   my $comp_seq = reverseComplement($sequence);
   my @c=split(//,$comp_seq);

   #remove k-mer from hash table and prefix tree
   delete $bin->{$o[0]}{$o[1]}{$o[2]}{$o[3]}{$o[4]}{$o[5]}{$o[6]}{$o[7]}{$o[8]}{$o[9]}{$o[10]}{$sequence};
   delete $bin->{$c[0]}{$c[1]}{$c[2]}{$c[3]}{$c[4]}{$c[5]}{$c[6]}{$c[7]}{$c[8]}{$c[9]}{$c[10]}{$comp_seq};
   delete $set->{$sequence};
   delete $seed->{$sequence};

   return $bin,$set,$seed;
}
#-----------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGC/TACG/;
   return (reverse());
}

#-----------------------
#Seems redundant, but try to save space in memory, esp since in this first release, PET isn't used to assist the overlap stage
sub readClone{

   my $file = shift;

   my $clone_number = 0;
   my $clone;

   open(IN,$file) || die "Can't open $file -- fatal\n";
   while(<IN>){
      chomp;

      if(/^([^\>]*)$/i){
         my $sdna = $1;

         if($sdna =~/([ACGT]*)\:([ACGT]*)/i){
            $clone_number++;
            $clone->{$clone_number}{'a'} = uc($1);
            $clone->{$clone_number}{'b'} = uc($2);
         }
      }
   }
   close IN;
   return $clone;
}

#-----------------
sub readFasta{
   my ($set,$bin,$file,$short,$paired) = @_;

   my $ctrd = 0;

   open(IN,$file) || die "Can't open $file -- fatal\n";
   open (SHO, ">>$short") || die "Can't write to $short -- fatal\n";

   while(<IN>){
      chomp;
      if(/^([^\>]*)$/i){
         my $sdna = $1;

         if($paired){
            if($sdna =~/^([ACGT]*)\:([ACGT]*)$/i){
               ($set,$bin,$ctrd) = &loadSequence($set,$bin,$ctrd,$1);
               ($set,$bin,$ctrd) = &loadSequence($set,$bin,$ctrd,$2);
            }else{
               my $pairing_failure_message = "The sequence \"$sdna\" is not in the right format for paired-end reads  -- Fatal\nMake sure your input is in the form (input sequences can be of variable lengths):\n\n>test\nGCTACGACTATGACATACAGT:GTAGATTGATCGCATGCACGCT\n\nWhere : separates paired reads.  Spaces in your input file might have caused this error.\n";
               print $pairing_failure_message;
               print LOG $pairing_failure_message;
               close LOG;
               exit;
            }
         }else{
            ($set,$bin,$ctrd) = &loadSequence($set,$bin,$ctrd,$1) if ($sdna =~ /^([ACGT]*)$/i);
         }
      }
   }
   close IN;
   close SHO;

   my $read_number_message = "$ctrd total sequences (" . keys( %$set ) . " unique)\n";
   printf $read_number_message;
   print LOG $read_number_message;

   return $set,$bin;
}

#-----------------
### added 31Jan08 R.Warren
sub loadSeed{

   my $file = shift;
   my $seed;  

   open(IN,$file) || die "Can't open $file -- fatal\n";
   
   my ($subseq,$prev)=('','');

   while(<IN>){
      chomp;

      if (/\>(\S+)/){
         my $head=$1;
         my $subseq_length = length($subseq);
         if($head ne $prev && $subseq ne '' && $subseq_length >= $MIN_READ_LENGTH && $subseq_length >= $min_overlap){
            $seed->{$subseq}{'count'}++;
            $seed->{$subseq}{'name'} = $prev;
            if($subseq=~/([NX])/i){print "WARNING: the fasta sequence >$prev in your seed file contains characters other than ACGT (i.e. $1) and may prevent proper contig extension.\n";}
         }
         $subseq='';
         $prev=$head;
      }elsif(/^([ACGTNX]*)$/i){
         $subseq .= uc($_);
      }
   }
   my $subseq_length = length($subseq);
   if($subseq ne '' && $subseq_length >= $MIN_READ_LENGTH && $subseq_length >= $min_overlap){
      $seed->{$subseq}{'count'}++;
      $seed->{$subseq}{'name'} = $prev;
      if($subseq=~/([NX])/i){print "WARNING: the fasta sequence >$prev in your seed file contains characters other than ACGT (i.e. $1) and may prevent proper contig extension.\n";}
   }
   
   close IN;

   return $seed;
}

#-----------------
sub loadSequence{

   my ($set,$bin,$ctrd,$seq) = @_;

   my $orig=uc($seq);
   my $orig_mer = length($orig);

   if ($orig ne '' && $orig_mer >= $MIN_READ_LENGTH && $orig_mer >= $min_overlap){
      ####show progress
      my $s10k = $ctrd / 10000;
      print "." if ($s10k == int($s10k) && $s10k && $ctrd);
      $|=1; ###clear buffer

      my $s100k = $ctrd / 100000;
      printf "%i sequences inputted\n", ($ctrd) if ($s100k == int($s100k) && $s100k && $ctrd);

      my @f=split(//,$orig);
      $set->{$orig}{'count'}++;

      $_ = $orig;
      tr/ACTG/TGAC/;
      my $rc=reverse();

      my @r=split(//,$rc);

      ### added 31Jan08 R.Warren
      $MAX=$orig_mer if ($orig_mer > $MAX);

      $bin->{$f[0]}{$f[1]}{$f[2]}{$f[3]}{$f[4]}{$f[5]}{$f[6]}{$f[7]}{$f[8]}{$f[9]}{$f[10]}{$orig}++;
      $bin->{$r[0]}{$r[1]}{$r[2]}{$r[3]}{$r[4]}{$r[5]}{$r[6]}{$r[7]}{$r[8]}{$r[9]}{$r[10]}{$rc}++;
      $ctrd++;
      
   }elsif($orig ne ''){
      if($orig_mer < $MIN_READ_LENGTH){
         print SHO "$seq\tInput sequence shorter than minimum read length allowed ($orig_mer < $MIN_READ_LENGTH nt)\n";
      }elsif($orig_mer < $min_overlap){
         print SHO "$seq\tInput sequence shorter than minimum overlap specified($orig_mer < -m $min_overlap)\n";
      }
   }
  
   $MAX = $MAX_TOP if ($MAX > $MAX_TOP);
   return $set,$bin,$ctrd;
}

### We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca
