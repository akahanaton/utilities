#!/usr/bin/perl
use strict;
use Getopt::Long;
my $gb_file;
my $cds_nofile;
my $mismatch;
my $seed;
GetOptions("i:s"=> \$gb_file,"u:s"=> \$cds_nofile,"l:s"=> \$seed,"m:s"=> \$mismatch);
if(!defined $gb_file or !defined $cds_nofile or !defined $mismatch or !defined $seed)
{
	die "
	***************************************************
	Script: get_gene_info.pl
	Useage: get_gene_info.pl [options]
	Available options: -i : Input cds info sequence File (fasta format)
	Available options: -u : Input reads no File (number of cds txt file)
	Available options: -l : reads lengths (example: 35bps)
	Available options: -M : erro rate (one mismatch per  bases,example: put 0.0001 if there is one mismatch per 10000bps)	
	***************************************************
	\n"; 
}
print "-i:$gb_file ";
print "-u:$cds_nofile ";
print "-m:$mismatch ";
$mismatch=1/$mismatch;
open(IN2,"$gb_file")||die"can not open fasta format file!";
#$gb_file=~/(\w+)\./;
$gb_file=~/(NC\_\d+)/;
my $dir=$1;
#open(IN,"cds_reads_NO.txt");
open(IN,"$cds_nofile")||die"can not open txt format file!";
my $output = "reads.fa";
my $output2 = "reads_info.txt";
open(OUT2,">./$dir/$output")||die"can not open output file!";
open(OUT,">./$dir/$output2")||die"can not open output file!";
#$mismatch=10000;
my @readnu;
@readnu=<IN>;
my $count=0;
my $readid=0;
#$seed=35;
my $i=0;
while(<IN2>)
{
	my $readseq="";
	my $line=$_;
	chomp($line);
	my @b=split(/\|/,$line);
	my $accnum=$b[1];
	my $loc=$b[6];
	chomp($loc);
  my	$seq=<IN2>;
	chomp($seq);
	my $len=length($seq);
  my 	$ggg=$readnu[$i];
	for(my $n=0;$n<$ggg;$n++)
	{
		$count+=$seed;
		srand;
		my $kk=$len-$seed;
	my $fl=abs(int(rand($len-$seed-2)));	
	my $readseq=substr($seq,$fl,$seed);
   $fl=$fl+1;
	 $readid++;
  if($count>=$mismatch)
  {
	my $loc2="";
	$count=0;
	my $tras=int(rand($seed));
	my $readseq1=substr($readseq,0,$tras-1);
  my $readseq2=substr($readseq,$tras,$seed);
	my $char1=substr($readseq,$tras-1,1);
	my $char2=get_mismatch($char1);
	$readseq=$readseq1."$char2".$readseq2;
	print OUT ">read$readid|cds",$i+1,"|";
	print OUT2 ">read$readid\n";
	$readseq=substr($readseq,0,$seed);
	if($b[2] eq -1)
	{
		$loc2=getloc_rev($loc,$fl,$seed);
		
	}
	else{
		$loc2=getloc($loc,$fl,$seed);
		}
	print OUT "$loc2|$accnum|$tras:$char1->$char2\n";
	print OUT2 "$readseq\n";
}
else{
	my $loc2="";
	print OUT ">read$readid|cds",$i+1,"|";
	print OUT2 ">read$readid\n";
	$readseq=substr($readseq,0,$seed);
	if($b[2] eq -1)
	{
		$loc2=getloc_rev($loc,$fl,$seed);
	}
	else{$loc2=getloc($loc,$fl,$seed);}
	print OUT "$loc2\n";
	print OUT2 "$readseq\n";
}
}
$i++;
#print "$i---->\n";
}
sub getloc()
{
	my($loc,$f,$seed)=@_;
my	@a=split(/\)/,$loc);
	 $loc="";
	my @b;
	foreach my $a(@a)
	{
		$a=~/(\d+)..(\d+)/;
				push(@b,$1,$2);
			}
			my $m=@b;
			for(my $n=0;$n<$m;$n=$n+2)
			{
				my $ff=$f+$b[$n]-1;
				if($ff+$seed-1<=$b[$n+1])
				{
					my $ll=$ff+$seed-1;
					$loc="($ff..$ll)";
					return $loc;
					#print OUT "$loc";
					last;
				}
				elsif($ff<=$b[$n+1] and ($ff+$seed-1)>$b[$n+1])
				{
					my $c=$b[$n+2]+$seed-2+$ff-$b[$n+1];
					$loc="($ff..$b[$n+1])($b[$n+2]..$c)";
					return $loc;
				#	print OUT "$loc";
					last;
				}
				$f=$f-$b[$n+1]+$b[$n]-1;
			}
}
sub getloc_rev()
{
	my($loc,$f,$seed)=@_;
	my	@a=split(/\)/,$loc);
 $loc="";
	my @b;
	foreach my $a(@a)
	{
		$a=~/(\d+)..(\d+)/;
				push(@b,$1,$2);
			}
			my $m=@b;
			for(my $n=0;$n<$m;$n=$n+2)
			{
				my $ff=$b[$n]-$f+1;
				if($ff-$seed+1>=$b[$n+1])
				{
					my $ll=$ff-$seed+1;
					$loc="($ll..$ff)";
					return $loc;					
					last;
				}
				elsif($ff>=$b[$n+1] and $ff-$seed+1<$b[$n+1])
				{
					#$c=$b[$n+1]-33+$ff-$b[$n+2];
					my $c=$b[$n+2]-$seed+2-$b[$n+1]+$ff;
					$loc="($c..$b[$n+2])($b[$n+1]..$ff)";
					return $loc;
					last;
				}
				$f=$f-$b[$n]+$b[$n+1]-1;
			}
}
sub get_mismatch()
{
	my($char1)=@_;
	my $char2;
	if ($char1 eq "A" or  $char1 eq "a")
	{
	my	$i=int(rand(2));
	
		if($i==0)
		{
			$char2="T";
		}
		elsif($i==1)
		{
			$char2="G";
		}
		else
		{
			$char2="C";
		}
	}
	elsif ($char1 eq "T" or  $char1 eq "t")
	{
	my	$i=int(rand(2));
		if($i==0)
		{
			$char2="A";
		}
		elsif($i==1)
		{
			$char2="G";
		}
		else
		{
			$char2="C";
		}
	}
	elsif ($char1 eq "G" or  $char1 eq "g")
	{
	my	$i=int(rand(2));
		if($i==0)
		{
			$char2="A";
		}
		elsif($i==1)
		{
			$char2="T";
		}
		else
		{
			$char2="C";
		}
	}
	else
	{
	my	$i=int(rand(2));
		if($i==0)
		{
			$char2="A";
		}
		elsif($i==1)
		{
			$char2="G";
		}
		else
		{
			$char2="T";
		}
	}
	return $char2;
}
=script
sub get_nomber()
{
	open(IN,"cds_normornd0_1_8233.txt");
open(IN2,"NC_003070.gbk.cds2.txt");
open(OUT,">cds_reads_NO.txt");
@a=<IN>;
$reads=10000000;
$sum=0;
while(<IN2>)
{
	$line=$_;
	chomp($line);
	@b=split(/\|/,$line);
	$len=@b;
	push(@c,$b[$len-1]);
	$sum+=$b[$len-1];
	$s=<IN2>;
	#print "$b[$len-1] ";
}
$i=0;
foreach $a(@a)
{
	push(@d,$a*$c[$i]/$sum);
	$sum2+=$a*$c[$i]/$sum;
	$i++;
}
print $sum2," ";
$sum3=0;
foreach $d(@d)
{
	$ff=$d*$reads/$sum2;
	$ff=$ff+0.5;
$ff=~/(\d+).\d+/;
$e=$1;
	print OUT $e,"\n";
	$sum3+=$e;
}
print "there are total :$sum3 reads created";
}
=cut