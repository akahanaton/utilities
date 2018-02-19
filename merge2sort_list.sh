#!/bin/sh
if [ $# -lt 1 ]; then
	echo "Usage : $0 list"
	echo "Notice: list only contain the pair_read 1"
	exit;
fi

list=$1

while read read1
do
	read2=${read1/1.fa/2.fq}
	dup_info=`merge2sort.pl $read1 $read2`
	dup_reads_num=`echo $dup_info | awk '{print $1 }'`
	all_reads_num=`echo $dup_info | awk '{print $6 }'`
	echo $dup_info
	echo "duplication ratio:"
	let "dup_ratio = dup_reads_num * 100 / all_reads_num"
	echo $dup_ratio
done < $list
