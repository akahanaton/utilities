#!/bin/sh
if [ $# -lt 2 ]; then
	echo "Usage : pair_reads.sh in_dir out_dir"
	exit;
fi

in_dir=$1
out_dir=$2

if [ ! -e $out_dir ]; then mkdir $out_dir; fi

for file in `ls $in_dir`
do
	if [ -e tmp_id ]; then rm tmp_id id1 id2 id3; fi

	grep '^>' $in_dir/$file | cut -d '/' -f1 | uniq -d > tmp_id

	sed -e 's/$/\/1/' tmp_id > id1
	sed -e 's/$/\/2/' tmp_id > id2
	cat id1 id2 > id3

	f1_file=$out_dir/${file/_1.fq.unmap/_1.fa}
	f2_file=$out_dir/${file/_1.fq.unmap/_2.fa}
	s_file=$out_dir/${file/_1.fq.unmap/.single}

	fnafile -i id1 -o $f1_file $in_dir/$file
	fnafile -i id2 -o $f2_file $in_dir/$file
	fnafile -e id3 -o $s_file $in_dir/$file
done

if [ -e tmp_id ]; then rm tmp_id id1 id2 id3; fi
