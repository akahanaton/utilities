#!/usr/bin/python
import time
import os
#--------------------------------------------------
# from pylab import *
#-------------------------------------------------- 
#--------------------------------------------------
# from Bio import SeqIO
#-------------------------------------------------- 


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

	def Set(self, hit_array): 
		if (len(hit_array) == 16):
			self.Query 		= hit_array[0]
			self.Sbjct 		= hit_array[4]
			#--------------------------------------------------
			# self.Identity 	= int(hit_array[8])
			#-------------------------------------------------- 
			self.AlignLen 	= int(hit_array[11])
			#--------------------------------------------------
			# self.Mismatch 	= int(hit_array[4])
			#-------------------------------------------------- 
			#--------------------------------------------------
			# self.Gaps 		= int(hit_array[5])
			#-------------------------------------------------- 
			self.QLen 		= int(hit_array[1])
			self.QStart 	= int(hit_array[2])
			self.QEnd 		= int(hit_array[3])
			self.SLen 		= int(hit_array[5])
			self.SStart 	= int(hit_array[6])
			self.SEnd 		= int(hit_array[7])
			#--------------------------------------------------
			# self.Score 		= float(hit_array[12])
			#-------------------------------------------------- 
			#--------------------------------------------------
			# self.RC 		= int(hit_array[13])
			#-------------------------------------------------- 
			self._empty    	= False 
		else:
			print "Set BlastResult Error, exiting"
			print hit_array
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


def ConstructPhrapList(all_blast_result, sbjct_id_uniq):
	#--------------------------------------------------
	# bac_end_info = {}.fromkeys(('RC','SStart','SEnd','QName','QStart','QEnd'))
	#-------------------------------------------------- 
	for each_sbjct in sbjct_id_uniq:
		hit_sbjct = 0
		hit_query = []
		for each_result in all_blast_result:
			if each_result.Sbjct == each_sbjct:
				hit_sbjct += 1
				hit_query.append(each_result)
				#--------------------------------------------------
				# tmp_bac_end_info = bac_end_info.copy()
				#-------------------------------------------------- 
		if hit_query != []:
			scaffold_info = '%s.  %d reads; %d bp' % (hit_query[0].Sbjct, len(hit_query), hit_query[0].SLen)
			print scaffold_info
			def my_cmp(E1, E2):
				return cmp(E1.SStart, E2.SStart)
			hit_query.sort(my_cmp)
			for result in hit_query:
				if result.SStart > result.SEnd:
					bac_end_info = 'C\t%d\t%d\t%s\t%d\t%d' % (result.SEnd, result.SStart, result.Query, result.QStart, result.QEnd)
					print bac_end_info
				else:
					bac_end_info = '\t%d\t%d\t%s\t%d\t%d' % (result.SStart, result.SEnd, result.Query, result.QStart, result.QEnd)
					print bac_end_info

if __name__ == '__main__':
	import sys
	def usage():
		print 'construct a phrap list for SCRA.sh'
		print 'usage:  python  construct_phrap_list.py   BlastResult.tab'
		print 'Notice: blast result must be sort according to and sbjct and sbjct-start pos'
		sys.exit(1)

	argc = len(sys.argv)
	if argc < 2:
		usage()


	BlastResultFile = sys.argv[1]

	BlastResult_h = open(BlastResultFile, 'r')
	AllBlastResult=[]
	SbjctIdUniq = []
	for eachline in BlastResult_h:
		Alist = eachline.rstrip().split()
		tmpBlastResult = BlastResult()
		tmpBlastResult.Set(Alist)
		AllBlastResult.append(tmpBlastResult)
		if tmpBlastResult.Sbjct not in SbjctIdUniq:
			SbjctIdUniq.append(tmpBlastResult.Sbjct)
	BlastResult_h.close()
	
	ConstructPhrapList(AllBlastResult,SbjctIdUniq)
