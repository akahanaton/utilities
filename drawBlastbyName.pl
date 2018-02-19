#!/usr/bin/perl	-w
use strict;
use	FindBin	qw($Bin);
use	lib	$Bin;
use	SVG;
use	Font;

unless 	(@ARGV){
	print	"Usage:	$0	blast_result  ref_length \n";
	exit;
}
open IN, "$ARGV[0]" or die"Can't	open	$ARGV[0]\n";
open LINE, ">$ARGV[0].SVG" or die "can not write into file";

 my $ref_length=$ARGV[1];
 #my $max=$ref_length/1000;
 my	$svg=	SVG->new('width',$ref_length,'height',1000);

 $svg->line('x1',10,'y1',60,'x2',$ref_length,'y2',60,'stroke','black','stroke-width',1);
my $i=0;
my $pre_name="";
my $color='red';
my $nextcolor='blue';
while (my $line =<IN>)
{ 
	my @info=split(/\s+/,$line);
	my $cur_name=$info[0];
	my $contigStart=$info[6]/100;
	my $contigEnd=$info[7]/100;
  
  if($pre_name ne $cur_name)
  { my $temp;
  	$temp=$color;
  	$color=$nextcolor;
  	$nextcolor=$temp;
  }
  
  $svg->line('x1',$contigStart,'y1',80+$i,'x2',$contigEnd,'y2',80+$i,'stroke',$color,'stroke-width',1);
  #$svg->text('x',$contigStart,'y',80+$i,'stroke','black','fill','black','-cdata',$cur_name,"transform"=>"rotate(90,$contigStart,80)");
	$i+=5;
	if($i>=50)
	{$i=$i%5;}
  $pre_name=$cur_name;
}
print LINE	$svg->xmlify();
close LINE;
close IN;
