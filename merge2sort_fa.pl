#! /usr/bin/perl -w

open (READ1, $ARGV[0]) && open (READ2, $ARGV[1]) || die ($!);
open (OUT,"| msort | uniq -c | awk '{if(\$2~/AAAAAAAAAAAAAAAAAAAA/)c+=\$1;if(\$1>1)b+=\$1;a+=\$1;}END{print b,\"/\",a,\" \",c}'") || die($!);
while (my $read1 =<READ1>){
	chomp $read1;
	if($read1=~/^\@/){$read1=~s/^\@/>/;}
	elsif($read1=~/^\+/){$read1=<READ1>;next;}
	if ($.%2==1){
		#print "$read1\n";
	}
	elsif ($.%2==0){
		print OUT "$read1";
		for(my $i=1; $i!=3;$i++){
			my $read2 = <READ2>;
			chomp $read2;
			if($read2=~/^\@/){$read2=~s/^\@/>/;}
        		elsif($read2=~/^\+/){$read2=<READ2>;$i--;next;}
		   	if ($i%2==1){
				#print "$read2\n"; 
			}
			elsif ($i%2==0){
				print OUT "$read2\n";
			}
		}
	}
}
close READ1;
close READ2;
close OUT;
