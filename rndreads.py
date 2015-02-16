#!/usr/bin/env python3

import sys
import random
import argparse
import SeqUtils
import io
import gzip

parser = argparse.ArgumentParser(description='Script to select reads based on various criteria')
parser.add_argument("readfiles",type=str,nargs="+")
parser.add_argument("-v","--verbose",action="store_true",default=False,help="Print errors to StdErr")
parser.add_argument("-p","--prefix",default="out")
parser.add_argument("-T","--trim",type=int)
parser.add_argument("-t","--tmp",action="store_true",help="Append .tmp to outfiles")
parser.add_argument("-Q","--offset",type=int)
parser.add_argument("-q","--qtrim",type=int,help="Quality trim cutoff")
parser.add_argument("-w","--window",type=int,default=1,help="window to use for quality trim")
parser.add_argument("-s","--single",action="store_true",default=False)
parser.add_argument("-g","--gzip",action="store_true")
parser.add_argument("-m","--minlength",type=int,default=1)
args=parser.parse_args()


if args.single==False:
    postfix=".fastq"
    if args.gzip:
        postfix += ".gz"
    if args.tmp:
        postfix += ".tmp"
    if args.gzip:
        of=[gzip.open(args.prefix+"_1"+postfix,"wb"),gzip.open(args.prefix+"_2"+postfix,"wb")]
    else:
        of=[open(args.prefix+"_1"+postfix,"wb"),open(args.prefix+"_2"+postfix,"wb")]  
else:
    of=[sys.stdout.buffer]

fh=[Utils.FqStream(f,inQualBase=args.offset) for f in args.readfiles]
i=0
while(True):
    try:
        fq=[next(f) for f in fh]
    except StopIteration:
        break
    if args.trim:
        [e.ltrim(args.trim) for e in fq]
    if args.qtrim:
        [e.qtrim(qual=args.qtrim,window=args.window) for e in fq]
    if  not (False in [(len(e) >= args.minlength) for e in fq ]):
        for f in fq:
            try:
                f.write(of[i%len(of)])
                #print(f,file=of[i%len(of)])
            except Utils.QualError:
                pass
            i+=1

[f.close() for f in of]


