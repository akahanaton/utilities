#!/usr/bin/python
import time
import os
from pylab import *
from Bio import SeqIO


class BlastResult:
	""" A single Blast Hit Result """
	def Clear(self): 
		self.Query 		= ''
		self.Sbjct 		= ''
		self.Identity 	= 0
		self.AlignLen 	= 0
		self.Mismatch 	= 0
		self.Gaps 		= 0
		self.QLen 		= 0
		self.QStart 	= 0
		self.QEnd 		= 0
		self.SLen 		= 0
		self.SStart 	= 0
		self.SEnd 		= 0
		self.Score 		= 0
		self.RC 		= 0
		self._empty    	= True 

	def __init__(self):
		self.Clear()
	def Empty(self):
		return self._empty

	def Set(self, hit_line): 
		Alist = hit_line.rstrip().split()
		if (len(Alist) == 14):
			self.Query 		= Alist[0]
			self.Sbjct 		= Alist[1]
			self.Identity 	= int(Alist[2])
			self.AlignLen 	= int(Alist[3])
			self.Mismatch 	= int(Alist[4])
			self.Gaps 		= int(Alist[5])
			self.QLen 		= int(Alist[6])
			self.QStart 	= int(Alist[7])
			self.QEnd 		= int(Alist[8])
			self.SLen 		= int(Alist[9])
			self.SStart 	= int(Alist[10])
			self.SEnd 		= int(Alist[11])
			self.Score 		= float(Alist[12])
			self.RC 		= int(Alist[13])
			self._empty    	= False 
		else:
			print "Set BlastResult Error, exiting"
			exit(0)
	def SwapQuerySbjct(self):
		def swap(t1, t2):
			return t2, t1
		self.Query, self.Sbjct = swap( self.Query, self.Sbjct )
		self.QLen, self.SLen = swap( self.QLen, self.SLen)
		self.QStart, self.SStart = swap( self.QStart, self.SStart)
		self.QEnd, self.SEnd = swap( self.QEnd, self.SEnd )


	def Print(self):
		if not self.Empty():
			out_str = '%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d' %\
					   (self.Query, self.Sbjct, self.Identity, self.AlignLen, self.Mismatch, self.Gaps,\
						self.QLen, self.QStart, self.QEnd, self.SLen, self.SStart,self.SEnd, self.Score, self.RC)
			print out_str

def N_count(seq):
	N_num = 0
	for c in seq:
		if(c == 'n' or c == 'N'):
			N_num += 1
	return N_num


def PlotSbjct(out_dir, format, sbjct_id, all_query_seq, all_sbjct_seq, all_blast_result,query_cov_cutoff, sbjct_cov_cutoff, save_tag):
	direction_err_num = 0
	position_err_num = 0
	SeqLineWidth = 8
	GapLineWidth = ErrorLineWidth = 3

	sbjct_len = 0
	all_blast_result_num = len(all_blast_result)
	for i in range(0,all_blast_result_num):
		if (all_blast_result[i].Sbjct == sbjct_id):
			sbjct_len = all_blast_result[i].SLen
			break

	QueryForSameSbjct = {}	
	sbjct_seq=[0] * sbjct_len
	for j in range(i,all_blast_result_num):
		if (all_blast_result[j].Sbjct == sbjct_id):
			if(all_blast_result[j].Query not in QueryForSameSbjct):
				ResultList=[]
				ResultList.append(all_blast_result[j])
				QueryForSameSbjct[all_blast_result[j].Query] = ResultList 
				sbjct_seq[ all_blast_result[j].SStart : all_blast_result[j].SEnd] = [1] * (all_blast_result[j].SEnd - all_blast_result[j].SStart)
			else:
				QueryForSameSbjct[all_blast_result[j].Query].append(all_blast_result[j])
				sbjct_seq[ all_blast_result[j].SStart : all_blast_result[j].SEnd] = [1] * (all_blast_result[j].SEnd - all_blast_result[j].SStart)
		else:
			break
	total_sbjct_cov = 100 * sum(sbjct_seq) / sbjct_len 

	units_len = sbjct_len / 1000
	QueryNum = len(QueryForSameSbjct)
	print "QueryNum " + str(QueryNum) 
	#plot reference seq
	SbjctY = 0
	#--------------------------------------------------
	# plot([0,sbjct_len],[SbjctY, SbjctY],'k', linewidth=SeqLineWidth)
	#-------------------------------------------------- 
	ref_seq = all_sbjct_seq[sbjct_id]
	if 'N' not in ref_seq and 'n' not in ref_seq:
		 plot([0,sbjct_len],[SbjctY, SbjctY],'k', linewidth=SeqLineWidth)
		 hold(True)
 	else:
		 N_Set = S_Set = False
		 N_begin = N_end = 0;
		#--------------------------------------------------
		#  plot([0,sbjct_len],[SbjctY, SbjctY],'k', linewidth=SeqLineWidth)
		#  hold(True)
		#-------------------------------------------------- 
		 for i in range( 0, sbjct_len ):
		    if ( ref_seq[i] == 'n' or ref_seq[i] == 'N' ):
				if not N_Set:
					N_Set = True
					N_begin = i
					plot([N_end, N_begin -1],[SbjctY, SbjctY],'k', linewidth=SeqLineWidth)
					hold(True)
					S_Set = False
		    else:
				if N_Set:
					N_Set = False
					S_Set = True
					if sbjct_len > 10000:
						N_end = i -1 + 400
					else:
						N_end = i -1
					plot([N_begin, N_end],[SbjctY, SbjctY ],'r', linewidth=GapLineWidth)
					hold(True)
		 if S_Set:
			 N_Set = True
			 N_begin = i
			 plot([N_end, N_begin -1],[SbjctY , SbjctY],'k', linewidth=SeqLineWidth)
	t_begin = time.localtime() 
	HighCovQuery = {}
	LowCovQuery = {}
	HighCov2Sbjct={}
	high_query_cov_index = 0
	GapBetweenQuery = 5
	LowCov2Sbjct_Y= SbjctY - GapBetweenQuery 
	HighCov2Sbjct_Y = SbjctY - GapBetweenQuery * 2 
	YticksPos=[SbjctY, LowCov2Sbjct_Y]
	YticksStr=['BAC', 'Repeat']
	cur_query_y = 0 
	sbjct_seq=[0] * sbjct_len
	for each_query in QueryForSameSbjct:
		query_seq=[0] * QueryForSameSbjct[each_query][0].QLen
		for each_blast_hit in QueryForSameSbjct[each_query]:
			query_seq[ each_blast_hit.QStart : each_blast_hit.QEnd] = [1] * (each_blast_hit.QEnd - each_blast_hit.QStart)
			sbjct_seq[ each_blast_hit.SStart : each_blast_hit.SEnd] = [1] * (each_blast_hit.SEnd - each_blast_hit.SStart)
		query_cov = 100 * sum(query_seq) / len(query_seq)
		if query_cov > query_cov_cutoff:
			HighCovQuery[ each_query ] = query_cov
		else:
			LowCovQuery[ each_query ] = query_cov
		sbjct_cov = 100 * sum(sbjct_seq) / sbjct_len
		if sbjct_cov > sbjct_cov_cutoff:
			print each_query
			HighCov2Sbjct[each_query] = sbjct_cov
			cur_query_y = HighCov2Sbjct_Y - high_query_cov_index * GapBetweenQuery
			high_query_cov_index += 1
			YticksStr.append(each_query)
			YticksPos.append(cur_query_y)
			blast_hit_num = len(QueryForSameSbjct[each_query])
			if(blast_hit_num > 1):
				cur_blast_hit = QueryForSameSbjct[each_query][0]
				#--------------------------------------------------
				# cur_blast_hit.Print()
				#-------------------------------------------------- 
				plot([cur_blast_hit.SStart, cur_blast_hit.SEnd], [cur_query_y, cur_query_y], 'b', linewidth=SeqLineWidth)
				if (cur_blast_hit.RC):
					if (cur_blast_hit.QEnd < cur_blast_hit.QLen):
						plot([cur_blast_hit.SStart-cur_blast_hit.QLen+cur_blast_hit.QEnd, cur_blast_hit.SStart -1 ], [cur_query_y, cur_query_y], 'y', linewidth=SeqLineWidth)
				else:
					if (cur_blast_hit.QStart > 1):
						plot([cur_blast_hit.SStart-cur_blast_hit.QStart, cur_blast_hit.SStart -1 ], [cur_query_y, cur_query_y], 'y', linewidth=SeqLineWidth)
				i = 1 
				while (i < blast_hit_num ):
					pre_blast_hit = QueryForSameSbjct[each_query][i-1]
					# remove the nest result according to sbjct pos
					cur_blast_hit = QueryForSameSbjct[each_query][i]
					while(i < blast_hit_num-1 and cur_blast_hit.SStart >= pre_blast_hit.SStart and cur_blast_hit.SEnd <= pre_blast_hit.SEnd):
						i += 1
						cur_blast_hit = QueryForSameSbjct[each_query][i]
					#--------------------------------------------------
					# cur_blast_hit.Print()
					#-------------------------------------------------- 
					plot([cur_blast_hit.SStart, cur_blast_hit.SEnd], [cur_query_y, cur_query_y], 'b', linewidth=SeqLineWidth)
					if (cur_blast_hit.RC == pre_blast_hit.RC == 0):
						#check werther it is a real gap
						gap_len = cur_blast_hit.QStart - pre_blast_hit.QEnd
						if(gap_len > 100):
							n_ratio = 100 * N_count(all_query_seq[each_query][pre_blast_hit.QEnd:cur_blast_hit.QStart]) / (cur_blast_hit.QStart - pre_blast_hit.QEnd)
							if(n_ratio > 50):
								if(pre_blast_hit.SEnd < cur_blast_hit.SStart):
									plot([pre_blast_hit.SEnd, cur_blast_hit.SStart], [cur_query_y, cur_query_y], 'cyan', linewidth=GapLineWidth)
							else:
								plot([pre_blast_hit.SEnd, cur_blast_hit.SStart], [cur_query_y, cur_query_y], 'magenta', linewidth=ErrorLineWidth)
								position_err_num += 1

					elif (cur_blast_hit.RC == pre_blast_hit.RC == 1):
						#check werther it is a real gap
						gap_len = pre_blast_hit.QStart - cur_blast_hit.QEnd
						if(gap_len > 100):
							n_ratio = 100 * N_count(all_query_seq[each_query][cur_blast_hit.QEnd:pre_blast_hit.QStart]) / (pre_blast_hit.QStart - cur_blast_hit.QEnd)
							if(n_ratio > 50):
								if(pre_blast_hit.SEnd < cur_blast_hit.SStart):
									plot([pre_blast_hit.SEnd, cur_blast_hit.SStart], [cur_query_y, cur_query_y], 'cyan', linewidth=GapLineWidth)
							else:
								plot([pre_blast_hit.SEnd, cur_blast_hit.SStart], [cur_query_y, cur_query_y], 'magenta', linewidth=ErrorLineWidth)
								position_err_num += 1
					else:
						direction_err_num += 1
					i += 1
				cur_blast_hit = QueryForSameSbjct[each_query][-1]
				plot([cur_blast_hit.SStart, cur_blast_hit.SEnd], [cur_query_y, cur_query_y], 'b', linewidth=SeqLineWidth)
				if (cur_blast_hit.RC):
					if (cur_blast_hit.QStart > 1):
						plot([cur_blast_hit.SEnd+1, cur_blast_hit.SEnd+cur_blast_hit.QStart-1 ], [cur_query_y, cur_query_y], 'y', linewidth=SeqLineWidth)
				else:
					if (cur_blast_hit.QEnd < cur_blast_hit.QLen):
						plot([cur_blast_hit.SEnd+1, cur_blast_hit.SEnd + cur_blast_hit.QLen - cur_blast_hit.QStart], [cur_query_y, cur_query_y], 'y', linewidth=SeqLineWidth)
			else:
				cur_blast_hit = QueryForSameSbjct[each_query][0]
				plot([cur_blast_hit.SStart, cur_blast_hit.SEnd], [cur_query_y, cur_query_y], 'b', linewidth=SeqLineWidth)
				if (cur_blast_hit.RC):
					if (cur_blast_hit.QEnd < cur_blast_hit.QLen):
						plot([cur_blast_hit.SStart-cur_blast_hit.QLen+cur_blast_hit.QEnd, cur_blast_hit.SStart -1 ], [cur_query_y, cur_query_y], 'y', linewidth=SeqLineWidth)
					if (cur_blast_hit.QStart > 1):
						plot([cur_blast_hit.SEnd+1, cur_blast_hit.SEnd+cur_blast_hit.QStart-1 ], [cur_query_y, cur_query_y], 'y', linewidth=SeqLineWidth)
				else:
					if (cur_blast_hit.QStart > 1):
						plot([cur_blast_hit.SStart-cur_blast_hit.QStart, cur_blast_hit.SStart -1 ], [cur_query_y, cur_query_y], 'y', linewidth=SeqLineWidth)
					if (cur_blast_hit.QEnd < cur_blast_hit.QLen):
						plot([cur_blast_hit.SEnd+1, cur_blast_hit.SEnd + cur_blast_hit.QLen - cur_blast_hit.QStart], [cur_query_y, cur_query_y], 'y', linewidth=SeqLineWidth)
		else:	
			for each_blast_hit in QueryForSameSbjct[each_query]:
				plot([each_blast_hit.SStart, each_blast_hit.SEnd], [LowCov2Sbjct_Y, LowCov2Sbjct_Y], 'r', linewidth=5)
				#each_blast_hit.Print()
		sbjct_seq=[0] * sbjct_len
	t_end = time.localtime()
	#--------------------------------------------------
	# print time.mktime(t_end) - time.mktime(t_begin)
	#-------------------------------------------------- 

	print '%s\t%d\t%d\t%d\t%d\t%d' % (sbjct_id, sbjct_len, total_sbjct_cov, high_query_cov_index,direction_err_num, position_err_num)

	y_lim_up = SbjctY + GapBetweenQuery
	if cur_query_y < LowCov2Sbjct_Y:
		y_lim_down = cur_query_y - GapBetweenQuery
	else:
		y_lim_down = LowCov2Sbjct_Y - GapBetweenQuery
	ylim(y_lim_down, y_lim_up)
	yticks( YticksPos, YticksStr) 
	xlim_buffer = sbjct_len * 0.2
	xlim( -xlim_buffer, sbjct_len + xlim_buffer)
	title(sbjct_id)
	if save_tag:
		if not os.path.exists(out_dir):
			os.mkdir(out_dir)	
		FigName = out_dir + '/' + sbjct_id + '.' + format
		savefig(FigName)
		hold(False)

if __name__ == '__main__':
	import sys
	def usage():
		print 'plot alignment info ref/bac VS contigs/scaffold'
		print 'usage:  python  plot_blast.py  BlastResult.tab  query.fasta sbjct.fasta  out_dir  picture.format'
		print 'Notice: blast result must be sort according to and sbjct and sbjct-start pos'
		sys.exit(1)

	argc = len(sys.argv)
	if argc < 6:
		usage()


	BlastResultFile = sys.argv[1]
	QuerySeqFile = sys.argv[2]
	SbjctSeqFile = sys.argv[3]
	OutDir = sys.argv[4]
	PicFormat = sys.argv[5]

	QuerySeqFile_h = open(QuerySeqFile)
	QuerySeq={}
	for seq_record in SeqIO.parse(QuerySeqFile_h, "fasta"):
		QuerySeq[seq_record.id] = seq_record.seq
	SbjctSeqFile_h = open(SbjctSeqFile)
	SbjctSeq={}
	for seq_record in SeqIO.parse(SbjctSeqFile_h, "fasta"):
		SbjctSeq[seq_record.id] = seq_record.seq

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

	print "SbjctID\tSbjctLen\tSbjctCov %\tCheckedQueryNum\tDirectionError\tMismatchError"
	for each_sbjct_id in SbjctIdUniq:
		PlotSbjct(OutDir,PicFormat,each_sbjct_id, QuerySeq,SbjctSeq, AllBlastResult,50,2, True)
