#!/bin/sh
if [ $# -lt 1 ]; then
	echo "Usage : make_grape_lib.sh  lib_prefix   SolexaReads_dir   SolexaReadsFilterN_dir"
	exit;
fi
echo 'max_rd_len=75'
prefix=$1
SolexaReadsDir=$2
SolexaReadsFilterNDir=$3
if [ ! -e $SolexaReadsFilterNDir ]; then
	mkdir $SolexaReadsFilterNDir
fi

Lib=`find $SolexaReadsDir/ -name "$prefix*"`
for LibDir in $Lib; do
	echo '[LIB]'
	SubDirName=`basename $LibDir`
	LibName="${SubDirName/$prefix/}"
	echo "name="$LibName
	if [ -e $SolexaReadsDir/InsertSize.txt ]; then
		InsertSize=`grep "$SubDirName" $SolexaReadsDir/InsertSize.txt | awk '{print $2}'`
		echo "min_ins="$InsertSize
		echo "avg_ins="$InsertSize
		echo "max_ins="$InsertSize
		if [ $InsertSize -ge 2000 ]; then
			echo "reverse_seq=1"
		else
			echo "reverse_seq=0"
		fi
	else
		echo "I can't find InsertSize.txt in "$SolexaReadsDir" exting now." 
		exit;
	fi
	SeqFileList=`find $LibDir/ -name "*fq"`
	for SeqFile in $SeqFileList; do
		if [[ $SeqFile = *_1.fq ]]; then
			q1=`basename $SeqFile`
			q2="${q1%*_1.fq}_2.fq"
			SeqDir="${SeqFile%/*}"
			OutDir="${SeqDir/$SolexaReadsDir/$SolexaReadsFilterNDir}"
			if [ ! -e $OutDir ]; then
				mkdir $OutDir
			fi
			if [ -e $SeqDir/$q2 ]; then
				#--------------------------------------------------
				# python N_fileter.py -a $SeqDir/$q1 -c $OutDir/filterN_$q1 -b $SeqDir/$q2 -d $OutDir/filterN_$q2 -q
				#-------------------------------------------------- 
				echo "q1="$OutDir"/filterN_"$q1
				echo "q2="$OutDir"/filterN_"$q2
			else
				#--------------------------------------------------
				# python N_fileter.py -a $SeqDir/$q1 -c $OutDir/filterN_$q1 -q
				#-------------------------------------------------- 
				echo "q="$OutDir"/filterN_"$q1
			fi
		elif [[ $SeqFile = *_2.fq ]]; then
			q2=`basename $SeqFile`
			q1="${q2%*_2.fq}_1.fq"
			SeqDir="${SeqFile%/*}"
			OutDir="${SeqDir/$SolexaReadsDir/$SolexaReadsFilterNDir}"
			if [ ! -e $OutDir ]; then
				mkdir $OutDir
			fi
			if [ -e $SeqDir/$q1 ]; then
				continue	
			else
				#--------------------------------------------------
				# python N_fileter.py -a $SeqDir/$q2 -c $OutDir/filterN_$q2 -q
				#-------------------------------------------------- 
				echo "q="$OutDir"/filterN_"$q2
			fi
		else
			q=`basename $SeqFile`
			SeqDir="${SeqFile%/*}"
			OutDir="${SeqDir/$SolexaReadsDir/$SolexaReadsFilterNDir}"
			if [ ! -e $OutDir ]; then
				mkdir $OutDir
			fi
			#--------------------------------------------------
			# python N_fileter.py -a $SeqDir/$q -c $OutDir/filterN_$q -q
			#-------------------------------------------------- 
			echo "q="$OutDir"/filterN_"$q
		fi
	done
done
