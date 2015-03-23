#!/usr/bin/env python3

import argparse
import sys
import random

# import numpy as np
 
# class Decoder(object):
#     def __init__(self, initialProb, transProb, obsProb):
#         self.N = initialProb.shape[0]
#         self.initialProb = initialProb
#         self.transProb = transProb
#         self.obsProb = obsProb
#         assert self.initialProb.shape == (self.N, 1)
#         assert self.transProb.shape == (self.N, self.N)
#         assert self.obsProb.shape[0] == self.N
         
#     def Obs(self, obs):
#         return self.obsProb[:, obs, None]
 
#     def Decode(self, obs):
#         trellis = np.zeros((self.N, len(obs)))
#         backpt = np.ones((self.N, len(obs)), 'int32') * -1
                 
#         # initialization
#         trellis[:, 0] = np.squeeze(self.initialProb * self.Obs(obs[0]))
                 
#         for t in xrange(1, len(obs)):
#             trellis[:, t] = (trellis[:, t-1, None].dot(self.Obs(obs[t]).T) * self.transProb).max(0)
#             backpt[:, t] = (np.tile(trellis[:, t-1, None], [1, self.N]) * self.transProb).argmax(0)
#         # termination
#         tokens = [trellis[:, -1].argmax()]
#         for i in xrange(len(obs)-1, 0, -1):
#             tokens.append(backpt[tokens[-1], i])
#         return tokens[::-1]

parser=argparse.ArgumentParser(description='Remove recombination from snp matrix')
parser.add_argument('alnfile',type=argparse.FileType(),help="fastafiles")
parser.add_argument('-p','--posfile',type=argparse.FileType(),default=None)
args =parser.parse_args()


class segment:
    ntest=0
    alnlength=0
    def __init__(self,h,start,end):
        self.hash = h
        self.start = start
        self.end = end
        segment.ntest +=1
        segment.alnlength = l
    def test(self):
        statistic= 1.0*segment.ntest*(hashprob[self.hash]/l*1.0)**(len(self)-1)
        if statistic < 0.05:
            print("Length: {}, stat: {}, ntest: {}, l: {}, hash: {}".format(len(self),statistic,segment.ntest,segment.alnlength,self.hash),file=sys.stderr)
        return statistic >= 0.05
    def __len__(self):
        return self.end-self.start
    def discard(self):
        hashprob[self.hash]-=len(self)
        segment.alnlength-=len(self)
        segment.ntest -=1

def hash(pattern):
    """ Determines if patterns are the same. Note that some nonidentical patterns may still give the same hash. """
    prev=""
    n=1
    hashstring=list()
    for i in pattern:
        if i == prev:
            n+=1
        else:
            hashstring.append(n)
            n=1
        prev=i
    hashstring.append(n)
    return tuple(hashstring[1:])

if args.posfile:
    posList=list()
    for line in args.posfile:
        posList.append(line.strip())
    args.posfile.close()

line=next(args.alnfile)
(norgs,alnlength)=line.split()
(norgs,alnlength)=(int(norgs),int(alnlength))
aln=list()
for i in range(alnlength):
    aln.append(list("N" * norgs))

forspalte=list()
org=0
for line in args.alnfile:
    line=line.strip()
    forspalte.append(line[0:10])
    pos=0
    for char in line[10:]:
        try:
            aln[pos][org]=char
        except IndexError:
            print("norgs=%i\norg=%i\nalnlength=%i\ni=%i" %(norgs,org,alnlength,i))
        pos+=1
    org+=1

""" Calculate frequencies of hash values """
hashprob=dict()

for profile in aln:
    h=hash("".join(profile))
    try:
        hashprob[h]+=1
    except KeyError:
        hashprob[h]=1
l=len(aln)
last=hash("")
lastdeleted=hash("")
sincelast=0

segments=list()
begin = None
for i in range(l):
    sh=hash("".join(aln[i]))
    if sh == last:
        """ A new segment is started """
        pass
    else:
        """ segment ends and is recorded """
        if begin:
            segments.append(segment(last,begin,i))
        last = sh
        begin = i
        
discardedPositions = set()
lastDiscarded=None
segments.sort(key=len,reverse=True)
for i in range(len(segments)):
    s = segments[i]
    if s.test()== False:
        s.discard()
        # if lastDiscarded and (i-lastDiscarded<4) and (s.hash == segments[lastDiscarded].hash):
        #     discardedPositions.update(range(segments[lastDiscarded].start,s.end))
        #     print("Discarding {} SNPs".format(s.end-segments[lastDiscarded].start),file=sys.stderr)
        #else:
        discardedPositions.update(range(s.start,s.end))
        print("Discarding {} SNPs".format(s.end-s.start),file=sys.stderr)
        lastDiscarded = i
print("%i %i" % (norgs,l-len(discardedPositions)))

for i in range(norgs):
    print(forspalte[i],end='')
    for j in range(l):
        if not j in discardedPositions:
            print(aln[j][i],end='')
    print('')

