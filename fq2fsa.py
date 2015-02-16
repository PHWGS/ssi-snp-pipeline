#!/usr/bin/env python3

import argparse
import SeqUtils

parser=argparse.ArgumentParser(description="Script to extract fragments of FastA files.")
parser.add_argument('fqfile',type=str)
args=parser.parse_args()

for fq in Utils.FqStream(args.fqfile):
    print(fq.fsa())

