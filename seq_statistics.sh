#!/bin/sh
if [ $# -lt 1 ]; then
	echo "Usage : seq_statistics.sh  seqfile"
	exit;
fi

seq=$1
len_info=``
for n in 90 80 70 60 50
do
	len_info=`seq_n50.pl $seq 100 $n`
	Nn0=`echo $len_info | awk '{print $10}'`
	seq_num=`seq_length_filter -f $seq -l $Nn0 | wc -l`
	printf "%s\t%s\n" $Nn0 $seq_num
done
	echo $len_info | awk '{print $2}'
	seq_length_filter -f $seq -l 100 | wc -l
	seq_length_filter -f $seq -l 2000 | wc -l
