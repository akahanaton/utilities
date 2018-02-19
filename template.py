#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
from BCBio import GFF
from collections import defaultdict
import argparse
import sys
import pprint
import GTF                      # https://gist.github.com/slowkow/8101481
from docopt import docopt
import pandas as pd
import numpy as np
import gzip
import time
from contextlib import contextmanager
import re

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

def uniq_by_field(infile, field_num):
# read all line
    infile_h = open(infile,'r')
    alllines= infile_h.readlines()
    infile_h.close()

    pre_key = ''
    for eachline in alllines:
            arr = eachline.rstrip().split()
            if arr[field_num -1] != pre_key:
                    print(eachline.rstrip())
            pre_key = arr[field_num-1]

# read line by line
    with open(infile) as f:
        for line in f:
            print(line)
    f.close()

@contextmanager
def log(message):
    """Log a timestamp, a message, and the elapsed time to stderr."""
    start = time.time()
    sys.stderr.write("{} # {}\n".format(time.asctime(), message))
    yield
    elapsed = int(time.time() - start + 0.5)
    sys.stderr.write("{} # done in {} s\n".format(time.asctime(), elapsed))
    sys.stderr.flush()


def count_bp(df):
    """Given a DataFrame with the exon coordinates from Gencode for a single
    gene, return the total number of coding bases in that gene.
    Example:
        >>> import numpy as np
        >>> n = 3
        >>> r = lambda x: np.random.sample(x) * 10
        >>> d = pd.DataFrame([np.sort([a,b]) for a,b in zip(r(n), r(n))], columns=['start','end']).astype(int)
        >>> d
           start  end
        0      6    9
        1      3    4
        2      4    9
        >>> count_bp(d)
        7

    Here is a visual representation of the 3 exons and the way they are added:

          123456789  Length
        0      ----       4
        1   --            2
        2    ------       6
            =======       7
    """
    start = df.start.min()
    end = df.end.max()
    bp = [False] * (end - start + 1)
    for i in range(df.shape[0]):
        s = df.iloc[i]['start'] - start
        e = df.iloc[i]['end'] - start + 1
        bp[s:e] = [True] * (e - s)
    return sum(bp)

def main(args):

    with log("Reading the Fasta  file: {}".format(args.fastaFile)):
        records = list(SeqIO.parse(args.fastaFile, "fasta"))

    with log("Reading the Gencode annotation file: {}".format(args.wigFile)):
        wig = pd.DataFrame.from_csv(args.wigFile,header=0, sep= " ", index_col=None)
        wig['CpG'] = ["CpG"] * (wig.size/2)
        #--------------------------------------------------
        # for raw in range(0,wig.size):
        #--------------------------------------------------

    wig.to_csv(args.outWigFile, header=True, index=None, sep=' ', mode='a')

    seqStr = dict()
    with log("Reading the fasta file: {}".format(args.fastaFile)):
        seqHandle = open(args.fastaFile, "rU")
        for record in SeqIO.parse(seqHandle, "fasta"):
            seqStr[record.id] = record.seq
    print(seqStr)

    with log("Reading the Gencode annotation file: {}".format(args.gffFile)):
        gc = GTF.dataframe(args.gffFile)
        #--------------------------------------------------
        # print(gc[1:10])
        #--------------------------------------------------

    # Select just exons of protein coding genes, and columns that we want to use.
    idx = (gc.feature == 'exon')
    exon = gc.ix[idx, ['seqname','start','end','ID','Parent']]
    exon['ID'] = exon['ID'].map(lambda x: re.sub(r'-mRNA.*','',x))
    # Convert columns to proper types.
    exon.start = exon.start.astype(int)
    exon.end = exon.end.astype(int)

    # Sort in place.
    exon.sort_values(['seqname','start','end'], inplace=True)

    # Group the rows by the Ensembl gene identifier (with version numbers.)
    groups = exon.groupby('ID')

    with log("Calculating coding region (exonic) length for each gene..."):
        lengths = groups.apply(count_bp)

    print(type(lengths))
    with log("Writing output file: {}".format(args.outFile)):
        lengths.to_csv(args.outFile, sep="\t", encoding="utf-8",index=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Count the number of coding base pairs in each gene.')
    parser.add_argument('wigFile', type=str, help='CpG posistion file')
    parser.add_argument('fastaFile', type=str, help='fasta file')
    parser.add_argument('outWigFile', type=str, help='Output file')
    parser.add_argument('-f', dest='fastaFile2', type=str, help='fasta file matching the gff file')
    #--------------------------------------------------
    # parser.add_argument('-c', dest='cleanCodon', type=str, default='Y', choices=['Y','N'], help='clean stop codon: Y/N')
    # parser.add_argument('-r', dest='cleanByRef', type=str, default=None, help='aligment format')
    #--------------------------------------------------
    args = parser.parse_args()
    main(args)
