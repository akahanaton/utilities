#!/usr/bin/python
import time
import os
from Bio import SeqIO
from pylab import *
import math




def PlotQual(out_dir, format, db_id, db_seq, all_kmre_result,save_tag):
	db_len = len(all_kmre_result[db_id])	
	db_depth = all_kmre_result[db_id]
	window_size = 10
	#--------------------------------------------------
	# window_num = db_len / window_size + 1
	# window_cov=[]
	# for i in range(window_num):
	# 	single_window_cov = sum(db_depth[i*window_size:(i+1)*window_size -1]) /window_size
	# 	window_cov.append(single_window_cov)
	#-------------------------------------------------- 
	
	#--------------------------------------------------
	# db_qual=[]
	# for i in all_kmre_result[db_id]:
	# 	if i > 100:
	# 		first = i / 100
	# 		end = i % 100
	# 		db_qual.append(first)
	# 		db_qual.append(end)
	# 	else:
	# 		db_qual.append(i)
	#-------------------------------------------------- 

	for i in range(0, db_len):
		if db_depth[i] > 0:
			db_depth[i] = math.log(db_depth[i],10)
	seq_length = range(0, db_len)
	plot(seq_length, db_depth)
	hold(True)
	plot([0, db_len],[1.47,1.47], 'k', linewidth =1 )
	plot([0, db_len],[1.78,1.78], 'k', linewidth =1 )

	#--------------------------------------------------
	# seq_length = range(db_len)
	# plot(seq_length, all_kmre_result[db_id])
	#-------------------------------------------------- 
	#--------------------------------------------------
	# db_len = len(db_depth)
	# seq_length = range(db_len)
	# plot(seq_length,db_depth)
	# hold(True)
	# plot([0, db_len ],[30,30], 'k', linewidth =1 )
	# plot([0, db_len ],[60,60], 'k', linewidth =1 )
	#-------------------------------------------------- 

	locs,labes=yticks()
	locs = range(int(max(locs)+0.5) + 1)
	locs.insert(1,1.48)
	#--------------------------------------------------
	# locs.insert(1,1.78)
	#-------------------------------------------------- 
	labes=[]
	for j in range(0, len(locs)):
		labes.append( str(int(10**locs[j])))
	labes[0] = '0'
	yticks(locs, labes)
	text(1.02, 0.5, 'Heterozygosis',horizontalalignment='left',	verticalalignment='center',	rotation='vertical',transform=gca().transAxes,clip_on=False)
	#--------------------------------------------------
	# ylabel('Heterozygosis')
	#-------------------------------------------------- 
	xlim_buffer = db_len * 0.2
	xlim(-xlim_buffer, db_len + xlim_buffer)
	if save_tag:
		if not os.path.exists(out_dir):
			os.mkdir(out_dir)	
		FigName = out_dir + '/' + db_id + '.' + format
		savefig(FigName)
		hold(False)

if __name__ == '__main__':
	import sys
	def usage():
		print 'plot depth info according to soap result'
		print 'usage:  python  plot_blast.py  soap.file database.fasta  out_dir  picture.format'
		sys.exit(1)

	argc = len(sys.argv)
	if argc < 5:
		usage()

	KmerFile = sys.argv[1]
	DbSeqFile = sys.argv[2]
	OutDir = sys.argv[3]
	PicFormat = sys.argv[4]

	DbSeqFile_h = open(DbSeqFile)
	DbSeq={}
	for seq_record in SeqIO.parse(DbSeqFile_h, "fasta"):
		DbSeq[seq_record.id] = seq_record.seq

	AllKmerResult={}
	DbIdUniq=[]
	KmerFile_h = open(KmerFile, 'r')
	while(1):
		id = KmerFile_h.readline().rstrip()[1:]
		Eachdepth = []
		if id == '':
			break
		DbIdUniq.append(id)
		kmerline=KmerFile_h.readline().rstrip().split()
		for depth in kmerline:
			Eachdepth.append(int(depth))
		AllKmerResult[id]=Eachdepth		
	KmerFile_h.close()
	for DbId in DbIdUniq: 
		PlotQual(OutDir,PicFormat,DbId, DbSeq, AllKmerResult, True)
