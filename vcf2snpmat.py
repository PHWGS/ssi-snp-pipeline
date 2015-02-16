#!/usr/bin/env python3

import re
import sys
import argparse

parser = argparse.ArgumentParser(description='Converts vcf files to fq files')
parser.add_argument("vcf",type=str,nargs='+')
parser.add_argument("-d","--min_depth",type=int,default=4,help="Minimum depth to accept a snp call")
parser.add_argument("-D","--max_depth",type=int,default=1000,help="Maximum depth to trust read alignment")
parser.add_argument("-1","--haploid",action="store_true")
parser.add_argument("-c","--clean",action="store_true",help="Discard ambiguous base calls")
parser.add_argument("-s","--strip",type=str,default="",help="Strip postfix from headers")
parser.add_argument("-p","--prefix_strip",type=str,default="",help="Strip prefix from headers")
args=parser.parse_args()

prestrip=re.compile('^{}'.format(args.prefix_strip))
poststrip=re.compile('{}$'.format(args.strip))
log=sys.stderr
comment=re.compile('^#')
indel=re.compile('INDEL')
snpref = re.compile('^([A-Za-z.])(,[A-Za-z])*$')
depth=re.compile('DP=(\d+)')
het = {'A' : 'A', 'T' : 'T', 'G' : 'G', 'C' : 'C', 
       'AC' : 'M', 'AG' : 'R', 'AT' : 'W', 'CA' : 'M', 'CG' : 'S', 'CT' : 'Y', 
       'GA' : 'R', 'GC' : 'S', 'GT' : 'K', 'TA' : 'W', 'TC' : 'Y', 'TG' : 'K',
       'CGT' : 'B', 'CTG' : 'B', 'GCT' : 'B', 'GTC' : 'B', 'TCG' : 'B', 'TGC' : 'B',
       'AGT' : 'D', 'ATG' : 'D', 'GAT' : 'D', 'GTA' : 'D', 'TAG' : 'D', 'TGA' : 'D',
       'CAT' : 'H', 'CTA' : 'H', 'ACT' : 'H', 'ATC' : 'H', 'TCA' : 'H', 'TAC' : 'H',
       'CGA' : 'V', 'CAG' : 'V', 'GCA' : 'V', 'GAC' : 'V', 'ACG' : 'V', 'AGC' : 'V'}

def GATKvcf2snpmat(fh,d,D,name):
    for line in fh:
        if comment.match(line):
            continue
        t=line.split('\t')
        chrom=t[0]
        pos=int(t[1])
        ref=t[3]
        alt=t[4]
        info=t[7]
        gap="-"
        bad="N"
        if args.haploid:
            alleles=ref + alt
        else:
            alleles=""
            chars=ref+"".join(alt.split(','))
            for i in range (len(chars)):
                for j in range(i,len(chars)):
                    alleles+= chars[i] if i==j else het[chars[i]+chars[j]] 
        try:
            plidx=t[8].split(":").index("PL")
            PL=t[9].split(":")[plidx].split(",")
            best=PL.index(min(PL))
            result=alleles[best]
        except ValueError:
            result=gap
            print("{}\t{}\t{}\tgap".format(name,chrom,pos),file=log)
        if args.clean and not result in ['A','T','G','C']:
            result=bad
        hit=depth.search(info)
        try:
            dp=int(hit.group(1))
            if(dp<d or dp>D):
                result=bad
                print("{}\t{}\t{}\tbad".format(name,chrom,pos),file=log)
        except AttributeError:
            result=bad
        try:
            if result in [bad,gap]:
                drop.add((chrom,pos))
            mat[(chrom,pos)].append(result)
        except KeyError:
            mat[(chrom,pos)]=[result]

drop=set()
mat=dict()
for f in args.vcf:
    fh=open(f)
    GATKvcf2snpmat(fh,args.min_depth,args.max_depth,f)

print("\t".join(["Pos"] + [prestrip.sub('',poststrip.sub('',i)) for i in args.vcf]))

for i in sorted(mat.keys()):
    if i not in drop:
        if len(set(mat[i])) >1 :
            print("\t".join(["___".join([i[0],str(i[1])])]+mat[i]))

