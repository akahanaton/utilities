#!/usr/bin/python
from Bio import SeqIO
import os


if __name__ == '__main__':
    import sys
    def usage():
        print('usage:  python3',sys.argv[0],'seq.file')
        sys.exit(1)

    #--------------------------------------------------
    # sys.argv = ["../tree.all", "all.longest.fna.aln.fna", "input.tree", "input.seq2", "eriEur"]
    #-------------------------------------------------- 
    argc = len(sys.argv)
    if argc < 2:
        usage()

    seqFile = sys.argv[1]

    bases = ['A','T','C','G','N','-']
    seqHandle = open(seqFile, "rU")

    for record in SeqIO.parse(seqHandle, "fasta"):
        curSeq=str(record.seq).upper()
        counter = dict.fromkeys(bases,0)
        try:
            for i in range(0, len(curSeq)):
                counter[curSeq[i]] = counter[curSeq[i]] + 1
        except ValueError:
            print("unkown seq\t", + record.id + "\t" +  curSeq[i])
        print(str(len(curSeq)) + "\t" + '\t'.join(["%s: %s" % (key, ('%('+key+')s') % counter) for key in bases]) + "\t" + str(counter['N'] + counter['-']) + "\t"  + record.id )
    seqHandle.close()
