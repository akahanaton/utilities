#!/usr/bin/python

if __name__ == '__main__':

	import sys
	def usage():
		print 'plot alignment info ref/bac VS contigs/scaffold'
		print 'usage:  python  plot_blast.py  BlastResult.tab   SoapResult.file   Qual.file  database.fasta   query.fasta   out_dir   picture.format'
		print 'Notice: blast result must be sort according to and sbjct and sbjct-start pos'
		sys.exit(1)

	argc = len(sys.argv)
	if argc < 7:
		usage()

	import plot_blast
	#--------------------------------------------------
	# from plot_soap_line import *
	#-------------------------------------------------- 
	import plot_soap_log 
	import plot_qual
	import time
	import os
	from pylab import *
	from Bio import SeqIO

	BlastResultFile = sys.argv[1]
	SoapResultFile = sys.argv[2]
	QualFile = sys.argv[3]
	DbSeqFile = sys.argv[4]
	QuerySeqFile = sys.argv[5]
	OutDir = sys.argv[6]
	PicFormat = sys.argv[7]
	
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
	for seq_record in SeqIO.parse(QuerySeqFile_h, "fasta"):
		QuerySeq[seq_record.id] = seq_record.seq

	BlastResult_h = open(BlastResultFile, 'r')
	AllBlastResult=[]
	SbjctIdUniq = []
	for eachline in BlastResult_h:
		tmpBlastResult = plot_blast.BlastResult()
		tmpBlastResult.Set(eachline)
		AllBlastResult.append(tmpBlastResult)
		if tmpBlastResult.Sbjct not in SbjctIdUniq:
			SbjctIdUniq.append(tmpBlastResult.Sbjct)
	BlastResult_h.close()

	#read soap result
	SoapResult ={}.fromkeys(('read_len','pos','db_name'))
	SoapIdUniq = []
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
		#--------------------------------------------------
		# else:
		# 	print eachline
		#-------------------------------------------------- 
		if tmpSoapResult['db_name'] not in SoapIdUniq :
			SoapIdUniq.append(tmpSoapResult['db_name'])
	SoapResult_h.close()

	AllQualResult={}
	QualIdUniq=[]
	QualFile_h = open(QualFile, 'r')
	while(1):
		id = QualFile_h.readline().rstrip()[1:]
		Eachdepth = []
		if id == '':
			break
		QualIdUniq.append(id)
		kmerline=QualFile_h.readline().rstrip().split()
		for depth in kmerline:
			Eachdepth.append(int(depth))
		AllQualResult[id]=Eachdepth		
	QualFile_h.close()

	if not os.path.exists(OutDir):
		os.mkdir(OutDir)	
	for each_sbjct_id in SbjctIdUniq:
		print each_sbjct_id
		print "SbjctID\tSbjctLen\tSbjctCov %\tCheckedQueryNum\tDirectionError\tMismatchError"
		subplot(311)
		if each_sbjct_id in SbjctIdUniq:
			plot_blast.PlotSbjct(OutDir,PicFormat, each_sbjct_id, QuerySeq, DbSeq, AllBlastResult,50,2,False)
		hold(False)
		subplot(312)
		if each_sbjct_id in SoapIdUniq:
			plot_soap_log.PlotDepth(OutDir,PicFormat,each_sbjct_id,DbLen,DbDepth,False)
		hold(False)
		subplot(313)
		if each_sbjct_id in QualIdUniq: 
			plot_qual.PlotQual(OutDir,PicFormat,each_sbjct_id, DbSeq, AllQualResult, False)
		hold(False)
		FigName = OutDir + '/' + each_sbjct_id + '.' + PicFormat 
		savefig(FigName)
