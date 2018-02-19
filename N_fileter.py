#!/usr/bin/python

import sys,getopt
import os

def usage():
	print "Filter the reads with N in it. "
	print "\tInput:  Suport both for fq or fa file format."
	print "\tOutput: fq to fa or fq, fa to fa"
	print "\tOptions:"
	print "\t\t-a, a file. "			
	print "\t\t-b, b file for paired end file format, optional."			
	print "\t\t-c, output file for a file. "			
	print "\t\t-d, output file for paired end b file. "			
	print "\t\t-q, fq file format for output, default is fa format."			
	
def FilterNSingle(in_seq_file, out_seq_file, fq_format):
	in_seq_file_h = open(in_seq_file,'r')
	out_seq_file_h = open(out_seq_file,'w')
	name_line = in_seq_file_h.readline()
	# determine werther is fq or fq format
	if(name_line[0] == '>'):
		while(name_line != ""):
			seq_line = in_seq_file_h.readline()
			if 'N' not in seq_line and 'n' not in seq_line:
				out_seq_file_h.write(name_line)
				out_seq_file_h.write(seq_line)
			name_line = in_seq_file_h.readline()
	elif(name_line[0] == '@'):
		while(name_line != ""):
			seq_line = in_seq_file_h.readline()
			name_line_q = in_seq_file_h.readline()
			qual_line = in_seq_file_h.readline()
			if 'N' not in seq_line and 'n' not in seq_line:
				if fq_format:
					out_seq_file_h.write(name_line)
					out_seq_file_h.write(seq_line)
					out_seq_file_h.write(name_line_q)
					out_seq_file_h.write(qual_line)
				else:
					out_seq_file_h.write(name_line.replace('@','>',1))
					out_seq_file_h.write(seq_line)
			name_line = in_seq_file_h.readline()
	else:
		print "Unrecognize file format:  " + in_seq_file + " exiting now.... "
		sys.exit()
	in_seq_file_h.close()
	out_seq_file_h.close()


def FilterNPair(in_seq_file_1, out_seq_file_1, in_seq_file_2, out_seq_file_2, fq_format):
	in_seq_file_1_h = open(in_seq_file_1,'r')
	out_seq_file_1_h = open(out_seq_file_1,'w')
	in_seq_file_2_h = open(in_seq_file_2,'r')
	out_seq_file_2_h = open(out_seq_file_2,'w')
	name_line_1 = in_seq_file_1_h.readline()
	# determine werther is fq or fq format
	if(name_line_1[0] == '>'):
		while(name_line_1 != ""):
			seq_line_1 = in_seq_file_1_h.readline()
			name_line_2 = in_seq_file_2_h.readline()
			seq_line_2 = in_seq_file_2_h.readline()
			if 'N' not in seq_line_1 and 'n' not in seq_line_1 and 'N' not in seq_line_2 and 'n' not in seq_line_2:
				out_seq_file_1_h.write(name_line_1)
				out_seq_file_1_h.write(seq_line_1)
				out_seq_file_2_h.write(name_line_2)
				out_seq_file_2_h.write(seq_line_2)
			name_line_1 = in_seq_file_h.readline()
	elif(name_line_1[0] == '@'):
		while(name_line_1 != ""):
			seq_line_1 = in_seq_file_1_h.readline()
			name_line_1_q = in_seq_file_1_h.readline()
			qual_line_1 = in_seq_file_1_h.readline()
			name_line_2 = in_seq_file_2_h.readline()
			seq_line_2 = in_seq_file_2_h.readline()
			name_line_2_q = in_seq_file_2_h.readline()
			qual_line_2 = in_seq_file_2_h.readline()
			if 'N' not in seq_line_1 and 'n' not in seq_line_1 and 'N' not in seq_line_2 and 'n' not in seq_line_2:
				if fq_format:
					out_seq_file_1_h.write(name_line_1)
					out_seq_file_1_h.write(seq_line_1)
					out_seq_file_1_h.write(name_line_1_q)
					out_seq_file_1_h.write(qual_line_1)
					out_seq_file_2_h.write(name_line_2)
					out_seq_file_2_h.write(seq_line_2)
					out_seq_file_2_h.write(name_line_2_q)
					out_seq_file_2_h.write(qual_line_2)
				else:
					out_seq_file_1_h.write(name_line_1.replace('@','>',1))
					out_seq_file_1_h.write(seq_line_1)
					out_seq_file_2_h.write(name_line_2.replace('@','>',1))
					out_seq_file_2_h.write(seq_line_2)
			name_line_1 = in_seq_file_1_h.readline()
	else:
		print "Unrecognize file format:  " + in_seq_file + " exiting now.... "
		sys.exit()
	in_seq_file_1_h.close()
	out_seq_file_1_h.close()
	in_seq_file_2_h.close()
	out_seq_file_2_h.close()
	
def main():
	try:
		opts,args = getopt.getopt(sys.argv[1:],"ha:b:c:d:q",["help"])
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)
	Input_1 = None
	Input_2 = None
	Output_1 = None
	Output_2 = None
	Out_fq = False
	for o,a in opts:
		if o in("-h","--help"):
			usage()
			sys.exit()
		elif o == "-q":
			Out_fq = True
		elif o == "-a":
			Input_1 = a
		elif o == "-b":
			Input_2 = a
		elif o == "-c":
			Output_1 = a
		elif o == "-d":
			Output_2 = a
		else:
			assert False, "unhandled option"
	if len(opts) < 1:
		usage()
	else:
		if (Input_2 == None):
			FilterNSingle(Input_1, Output_1, Out_fq)
		else:
			FilterNPair(Input_1, Output_1,Input_2, Output_2, Out_fq)
if __name__ == "__main__":
	main()
