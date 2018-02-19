#!/usr/bin/python

if __name__ == '__main__':

	import sys
	def usage():
		print 'plot alignment info ref/bac VS contigs/scaffold'
		print 'usage:  python  plot_blast.py  BlastResult.tab SoapResult.file database.fasta query.fasta out_dir  picture.format'
		print 'Notice: blast result must be sort according to and sbjct and sbjct-start pos'
		sys.exit(1)

	argc = len(sys.argv)
	if argc < 6:
		usage()

	from plot_blast import *
	#--------------------------------------------------
	# from plot_soap_line import *
	#-------------------------------------------------- 
	from plot_soap_log import *
	import time
	import os
	from pylab import *
	from Bio import SeqIO

	BlastResultFile = sys.argv[1]
	SoapResultFile = sys.argv[2]
	DbSeqFile = sys.argv[3]
	QuerySeqFile = sys.argv[4]
	OutDir = sys.argv[5]
	PicFormat = sys.argv[6]
	
	DbSeqFile_h = open(DbSeqFile)
	DbSeq={}
	DbDepth={}
	DbLen={}
	for seq_record in SeqIO.parse(DbSeqFile_h, "fasta"):
		DbSeq[seq_record.id] = seq_record.seq
		DbLen[seq_record.id] = len(seq_record.seq)
		DbDepth[seq_record.id] = [0]*DbLen[seq_record.id]

	QuerySeqFile_h = open(QuerySeqFile)
	QuerySeq={}
	for seq_record in SeqIO.parse(DbSeqFile_h, "fasta"):
		QuerySeq[seq_record.id] = seq_record.seq

	BlastResult_h = open(BlastResultFile, 'r')
	AllBlastResult=[]
	SbjctIdUniq = []
	for eachline in BlastResult_h:
		tmpBlastResult = BlastResult()
		tmpBlastResult.Set(eachline)
		AllBlastResult.append(tmpBlastResult)
		if tmpBlastResult.Sbjct not in SbjctIdUniq:
			SbjctIdUniq.append(tmpBlastResult.Sbjct)
	BlastResult_h.close()

	#read soap result
	SoapResult ={}.fromkeys(('read_len','pos','db_name'))
	DbIdUniq = []
	tmpSoapResult=SoapResult.copy()
	SoapResult_h = open(SoapResultFile, 'r')
	for eachline in SoapResult_h:
		SoapList = eachline.rstrip().split()
		tmpSoapResult['read_len'] = int(SoapList[5])
		tmpSoapResult['db_name'] = SoapList[7]
		tmpSoapResult['pos'] = int(SoapList[8])
		if tmpSoapResult['pos'] < DbLen[tmpSoapResult['db_name']]- 75:
			for i in range(0, tmpSoapResult['read_len']) :
				DbDepth[tmpSoapResult['db_name']][tmpSoapResult['pos'] + i - 1] += 1
		else
			print eachline
		if tmpSoapResult['db_name'] not in DbIdUniq :
			DbIdUniq.append(tmpSoapResult['db_name'])
	SoapResult_h.close()

	print "SbjctID\tSbjctLen\tSbjctCov %\tCheckedQueryNum\tDirectionError\tMismatchError"
	if not os.path.exists(OutDir):
		os.mkdir(OutDir)	
	for each_sbjct_id in DbLen:

		subplot(211)
		if each_sbjct_id in SbjctIdUniq:
			PlotSbjct(OutDir,PicFormat, each_sbjct_id, QuerySeq, DbSeq, AllBlastResult,50,2,False)
		hold(False)

		subplot(212)
		if each_sbjct_id in DbIdUniq:
			PlotDepth(OutDir,PicFormat,each_sbjct_id,DbLen,DbDepth,False)
		hold(False)

		FigName = OutDir + '/' + each_sbjct_id + '.' + PicFormat 
		savefig(FigName)
