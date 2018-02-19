#!/usr/bin/python
import time
import os
import re
from Bio import SeqIO

def print_seq(id,seq):
	print '>%s' % id
	seq_len = len(seq)
	if seq_len < 100:
		print seq
	else:
		i = 0
		for i in range(seq_len/100):
			print seq[i*100:(i+1)*100]
		print seq[(i+1)*100:]
	

if __name__ == '__main__':
	import sys
	def usage():
		print 'usage:  python  fna_index.py  contig.fasta prefix' 
		sys.exit(1) 

	argc = len(sys.argv)
	if argc < 2:
		usage()

	InFile = sys.argv[1]
	Prefix = sys.argv[2]

	index = 1
	
	InFile_h = open(InFile)
	for seq_record in SeqIO.parse(InFile_h, "fasta"):
		new_id = Prefix + str(index)
		print_seq(new_id, seq_record.seq)
		index += 1 
