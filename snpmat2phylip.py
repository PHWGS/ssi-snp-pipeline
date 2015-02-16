#!/usr/bin/env python3

import argparse
import sys

parser=argparse.ArgumentParser(description='Make a phylip infile from a snp matrix')
parser.add_argument('matrix',type=argparse.FileType(),help="fastafiles")
parser.add_argument('-n','--norgs',type=int,default=None)
parser.add_argument('-p','--posfile',type=str,default=None)
parser.add_argument('-s','--sedscript',type=str,default=None)
args =parser.parse_args()

if args.norgs:
    header=next(args.matrix).strip().split("\t",maxsplit=args.norgs+1)[:-1] 
else:
    header=next(args.matrix).strip().split("\t")
if args.sedscript:
    ss=open(args.sedscript,'w')
    for i in range(1,len(header)):
        replacement = "seq{:06d}".format(i)
        print("s/{}/{}/".format(replacement,header[i]),file=ss)
        header[i]=replacement
    ss.close()
        
if args.posfile:
    pf=open(args.posfile,'w')
aln=list()
for line in args.matrix:
    line=line.strip()
    if args.norgs:
        fields=line.split("\t",maxsplit=args.norgs+1)[:-1]
    else: 
        fields=line.split("\t")
    if len(fields) != len(header):
        print("Wrong number of entries in line (expected {} fields)\n{}".format(len(header),fields),file=sys.stderr)
    else:
        aln.append(fields)
    if args.posfile:
        print(line.split("\t",maxsplit=1)[0],file=pf)
if args.posfile:
    pf.close()
print("{} {}".format(len(header)-1,len(aln)))

for org in range(len(header[1:])):
    print("{: <9s} ".format(header[org+1][:10]),end='')
    for pos in range(len(aln)):
        print(aln[pos][org+1].upper(),end='')
    print("")
