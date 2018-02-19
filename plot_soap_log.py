#!/usr/bin/python
import time
import os
from Bio import SeqIO
from pylab import *
import math

def PlotDepth(out_dir, format, db_id, db_len, all_db_depth,save_tag):
	db_len = db_len[db_id]
	db_depth = all_db_depth[db_id]
	for i in range(0, db_len):
		if db_depth[i] > 0:
			db_depth[i] = math.log(db_depth[i], 10)
	seq_length = range(0, db_len)
	plot(seq_length, db_depth)
	hold(True)
	plot([0, db_len],[1.47,1.47], 'k', linewidth =1 )
	plot([0, db_len],[1.78,1.78], 'k', linewidth =1 )

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


	xlim_buffer = db_len * 0.2
	xlim(-xlim_buffer, db_len + xlim_buffer)
	text(1.02, 0.5, 'Depth',horizontalalignment='left',	verticalalignment='center',	rotation='vertical',transform=gca().transAxes,clip_on=False)
	#--------------------------------------------------
	# ylabel('Depth')
	#-------------------------------------------------- 
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

	SoapResultFile = sys.argv[1]
	DbSeqFile = sys.argv[2]
	OutDir = sys.argv[3]
	PicFormat = sys.argv[4]

	DbSeqFile_h = open(DbSeqFile)
	DbDepth={}
	DbLen={}
	for seq_record in SeqIO.parse(DbSeqFile_h, "fasta"):
		DbLen[seq_record.id] = len(seq_record.seq)
		DbDepth[seq_record.id] = [0]*DbLen[seq_record.id]

	SoapResult ={}.fromkeys(('read_len','pos','db_name'))
	DbIdUniq = []
	tmpSoapResult=SoapResult.copy()
	SoapResult_h = open(SoapResultFile, 'r')
	for eachline in SoapResult_h:
		SoapList = eachline.rstrip().split()
		tmpSoapResult['read_len'] = int(SoapList[5])
		tmpSoapResult['db_name'] = SoapList[7]
		tmpSoapResult['pos'] = int(SoapList[8])
		for i in range(0, tmpSoapResult['read_len']) :
			if tmpSoapResult['pos'] < DbLen[tmpSoapResult['db_name']]- 75:
				DbDepth[tmpSoapResult['db_name']][tmpSoapResult['pos'] + i - 1] += 1
		if tmpSoapResult['db_name'] not in DbIdUniq :
			DbIdUniq.append(tmpSoapResult['db_name'])
	SoapResult_h.close()

	for DbId in DbIdUniq: 
		PlotDepth(OutDir,PicFormat,DbId,DbLen,DbDepth,True)
