#!/usr/bin/perl -w
if(@ARGV<1)
  {
    print "$0 fqfilelist\n";
    ##two paired fq should be placed one after another
    exit;
  }
 open(LIST,"$ARGV[0]") or die;
 while($line=<LIST>)
   {
     $line=~s/\s+$//;
     $line1=<LIST>;
     $line1=~s/\s+$//;
     $i=0;
     open(IN,"$line") or die;
     open(IN1,"$line1") or die;
     while($seq=<IN>)
       {
         $seq1=<IN1>;
         $i++;
         if($i%4==2)
           {
            # $tag=substr($seq,0,8).substr($seq,16,8).substr($seq1,0,8).substr($seq1,16,8);  
            # $tag1=substr($seq1,0,8).substr($seq1,16,8).substr($seq,0,8).substr($seq,16,8);
            $tag=$seq.$seq1;
            $tag1=$seq1.$seq;
             if(defined($fre{$tag}))
               {
                 $fre{$tag}++;
               }elsif(defined($fre{$tag1}))  
               {
                 $fre{$tag1}++;
               }else
               {
                 $fre{$tag}++;
               }
           }
       } 
     close (IN);
     close (IN1);
   }
  seek (LIST,0,0);
  while($line=<LIST>)
   {
      $line=~s/\s+$//;
      $line1=<LIST>;
      $line1=~s/\s+$//; 
      $i=0;
      open(IN,"$line") or die;
      open(IN1,"$line1") or die;
      @array1=split /\//,$line;
      @array2=split /\//, $line1;
      $outfile1=$array1[$#array1].".fasta";
      $outfile2=$array2[$#array2].".fasta";
      open(OUT1,">$outfile1") or die;
      open(OUT2,">$outfile2") or die; 
      while($seq=<IN>)
       {
         $seq1=<IN1>;
         $i++;
         if($i%4==2)
           {
            # $tag=substr($seq,0,8).substr($seq,16,8).substr($seq1,0,8).substr($seq1,16,8);  
            $tag=$seq.$seq1;
             if($fre{$tag}==1)
               {
               	 print OUT1 "$beforeline$seq";
               	 print OUT2 "$beforeline1$seq1";
               }elsif($fre{$tag}>1)
               {
               	  print OUT1 "$beforeline$seq";
               	  print OUT2 "$beforeline1$seq1";
          #     	  print "$fre{$tag}\n";
               	  $fre{$tag}=0;
               	  
               }
           }elsif($i%4==1)
           {
           	 $beforeline=$seq;
           	 $beforeline1=$seq1;
           	 $beforeline=~s/@/>/;
           	 $beforeline1=~s/@/>/;
           }	 
       }
      close IN ;
      close IN1 ;
      close OUT1 ;
      close OUT2 ;
   }           	    	 	 
