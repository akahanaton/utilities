#!/usr/bin/python
import os

def uniq_by_field(infile, field_num):
	infile_h = open(infile,'r')
	alllines= infile_h.readlines()
	infile_h.close()

	pre_key = ''
	for eachline in alllines:
		arr = eachline.rstrip().split()
		if arr[field_num -1] != pre_key:
			print eachline.rstrip()
		pre_key = arr[field_num-1]

if __name__ == '__main__':
	import sys
	def usage():
		print 'uniq eachline by field, field begin with 1'
		print 'usage:  python  uniq_field.py  input_file  field_num'
		sys.exit(1)

	argc = len(sys.argv)
	if argc < 3:
		usage()
	
	Infile = sys.argv[1]
	Field_num = int(sys.argv[2])
	uniq_by_field(Infile,Field_num)


