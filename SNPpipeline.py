#!/usr/bin/env python3

import re
import os
import urllib.request
import subprocess as sp
import sys
import tempfile
import argparse
import threading
import SeqUtils
import time
import glob
#import ThreadPool
from multiprocessing.pool import ThreadPool

javapath     = "java"
GATKpath     = "/home/kiil/Build/GenomeAnalysisTK-2.5-2-gf57256b/GenomeAnalysisTK.jar"
bwapath      = "/home/kiil/Build/bwa-0.7.4/bwa"
picardpath   = "/home/kiil/Build/picard-tools-1.92/"
samtoolspath = "samtools"
bcftoolspath = "bcftools"
vcfutilspath = "/usr/share/samtools/vcfutils.pl"
soappath     = "SOAPdenovo-63mer"
edenapath    = "edena"
velvetpath   = "velvet"
elpreppath   = "elprep"

devnull=open(os.devnull,"w") # Use to discard output

class RunError (Exception):
    def __init__(self,message=""):
        print(message,file=sys.stderr)
        self.message=message

def untmp(path):
    pMv=sp.Popen(["mv",".".join([path,"tmp"]),path])
    retVal=pMv.wait()
    if retVal!=0:
        raise IOError
    if os.path.exists(".".join([path,"tmp","idx"])):
        pMv=sp.Popen(["mv",".".join([path,"tmp","idx"]),".".join([path,"idx"])])
        retVal=pMv.wait()
        if retVal!=0:
            raise IOError


class Read:
    """ Parent class for read classes that defines methods for per read analysis steps """
    def __init__(self):
        self.A_OK=True
    
    def callSNPs1(self):
        self.A_OK=False
        print("Calling SNPs round 1...")
        if args.haploid:
            #Use the haploid model:
            cmd=[javapath,"-jar",GATKpath,"-T","UnifiedGenotyper","-nt",str(self.parent.factor),"-I",".".join([self.aln,"sorted.dedup.bam"]),"-R",self.parent.path,"-o",".".join([self.aln,"vcf","tmp"]),"-ploidy","1","--genotype_likelihoods_model","SNP","--annotateNDA","-maxAltAlleles","3"]
            pGenotyper=sp.Popen(cmd,stderr=devnull,stdout=errorlog)
        else:
            #Use the diploid model, accounting for repeats and sequencing errors
            cmd=[javapath,"-jar",GATKpath,"-T","UnifiedGenotyper","-nt",str(self.parent.factor),"-I",".".join([self.aln,"sorted.dedup.bam"]),"-R",self.parent.path,"-o",".".join([self.aln,"vcf","tmp"]),"-ploidy","2","--genotype_likelihoods_model","SNP","--annotateNDA","-maxAltAlleles","6"]
            pGenotyper=sp.Popen(cmd,stderr=errorlog,stdout=devnull)

        pGenotyper.communicate()
        retVal=pGenotyper.wait()
        if retVal==0:
            untmp(".".join([self.aln,"vcf"]))
        else:
            raise RunError(" ".join(cmd))

    def reducebam(self):
        #Make reduced.bam
        self.A_OK=False
        print("Reducing reads...")
        cmd=[javapath,"-jar",GATKpath,"-T","ReduceReads","-I",".".join([self.aln,"sorted.dedup.bam"]),"-R",self.parent.path,"-o",".".join([self.aln,"reduced.bam.tmp"])]
        pReduce=sp.Popen(cmd,stderr=devnull)
        pReduce.communicate()
        retVal=pReduce.wait()
        if retVal==0:
            untmp(".".join([self.aln,"reduced.bam"]))
            pMv=sp.Popen(["mv",".".join([self.aln,"reduced.bam.tmp.bai"]),".".join([self.aln,"reduced.bam.bai"])])
            retVal=pMv.wait()
            if retVal!=0:
                raise IOError
        else:
            raise RunError(" ".join(cmd))

    def deduplicate(self):
        self.A_OK=False
        print("Deduplicating...")
        cmd=[javapath,"-jar",os.path.join(picardpath,"MarkDuplicates.jar"),"I=",".".join([self.aln,"sorted.bam"]),"O=",".".join([self.aln,"sorted.dedup.bam.tmp"]),"M=marked_duplicates","VALIDATION_STRINGENCY=LENIENT"]
        pDedup=sp.Popen(cmd,stderr=errorlog)
        pDedup.communicate()
        retVal=pDedup.wait()
        if retVal==0:
            untmp(".".join([self.aln,"sorted.dedup.bam"]))
        else:
            raise RunError(" ".join(cmd))
        cmd=[javapath,"-jar",os.path.join(picardpath,"BuildBamIndex.jar"), "I=",".".join([self.aln,"sorted.dedup.bam"]),"VALIDATION_STRINGENCY=LENIENT"]
        pBuildIndex=sp.Popen(cmd,stderr=devnull)
        pBuildIndex.communicate()
        retVal=pBuildIndex.wait()
        if retVal!=0:
            raise RunError(" ".join(cmd))

    def runbwamem(self):
        self.A_OK=False
        print("Running bwa...")
        fastq="fastq"
        if args.gzip:
            fastq+=".gz"
        if self.pairs:
            cmd_bwa=[bwapath,"mem","-M","-R","@RG\tID:1\tPL:illumina\tPU:barcode\tLB:myLibrary\tSM:20043495",self.parent.path,".".join([self.basename+"_1",fastq]),".".join([self.basename+"_2",fastq])]
        else:
            cmd_bwa=[bwapath,"mem","-p","-M","-R","@RG\tID:1\tPL:illumina\tPU:barcode\tLB:myLibrary\tSM:20043495",self.parent.path,".".join([self.basename+"_1",fastq])]
        pBwaMEM=sp.Popen(cmd_bwa,stdout=sp.PIPE,stderr=devnull)
        cmd_unmatched=["perl","-nae",'if($F[1] & 4){print STDERR $_;}else{print $_;}']
        unmatched=open(".".join([self.aln,"unmatched","sam"]),'w')
        pUnmatched=sp.Popen(cmd_unmatched,stdin=pBwaMEM.stdout,stderr=unmatched,stdout=sp.PIPE)
        cmd_sort=[javapath,"-jar",os.path.join(picardpath,"SortSam.jar"),"I=","/dev/stdin","O=",".".join([self.aln,"sorted.bam.tmp"]),"SO=coordinate","VALIDATION_STRINGENCY=LENIENT"]
        pSort=sp.Popen(cmd_sort,stdin=pUnmatched.stdout,stderr=errorlog)
        retVal=pSort.wait()
        if retVal==0:
            untmp(".".join([self.aln,"sorted.bam"]))
        else:
            raise RunError(" | ".join([" ".join(cmd_bwa)," ".join(cmd_unmatched)," ".join(cmd_sort)]))
        cmd=[javapath,"-jar",os.path.join(picardpath,"BuildBamIndex.jar"), "I=",".".join([self.aln,"sorted.bam"])]
        pBuildIndex=sp.Popen(cmd,stderr=errorlog)
        retVal=pBuildIndex.wait()
        if retVal!=0:
            raise RunError(" ".join(cmd))

    def runbwamemelprep(self):
        self.A_OK=False
        print("Running bwa...")
        fastq="fastq"
        if args.gzip:
            fastq+=".gz"
        if self.pairs:
            cmd_bwa=[bwapath,"mem","-M","-R","@RG\tID:1\tPL:illumina\tPU:barcode\tLB:myLibrary\tSM:20043495",self.parent.path,".".join([self.basename+"_1",fastq]),".".join([self.basename+"_2",fastq])]
        else:
            cmd_bwa=[bwapath,"mem","-p","-M","-R","@RG\tID:1\tPL:illumina\tPU:barcode\tLB:myLibrary\tSM:20043495",self.parent.path,".".join([self.basename+"_1",fastq])]
        pBwaMEM=sp.Popen(cmd_bwa,stdout=sp.PIPE,stderr=devnull)
        cmd_elprep=["elprep","--filter-unmapped-reads","--mark-duplicates","remove","--clean-sam","--sorting-order","coordinate","--nr-of-threads",str(3),"/dev/stdin",".".join([self.aln,"sorted.dedup.tmp.bam"])]
        #cmd_unmatched=["perl","-nae",'if($F[1] & 4){print STDERR $_;}else{print $_;}']
        #unmatched=open(".".join([self.aln,"unmatched","sam"]),'w')
        #pUnmatched=sp.Popen(cmd_unmatched,stdin=pBwaMEM.stdout,stderr=unmatched,stdout=sp.PIPE)
        pElprep=sp.Popen(cmd_elprep,stdin=pBwaMEM.stdout,stderr=errorlog)
        retVal=pElprep.wait()
        sp.Popen(["mv",".".join([self.aln,"sorted.dedup.tmp.bam"]),".".join([self.aln,"sorted.dedup.bam.tmp"])]).wait()
        if retVal==0:
            untmp(".".join([self.aln,"sorted.dedup.bam"]))
        else:
            raise RunError(" | ".join([" ".join(cmd_bwa)," ".join(cmd_elprep)]))
        cmd=[javapath,"-jar",os.path.join(picardpath,"BuildBamIndex.jar"), "I=",".".join([self.aln,"sorted.dedup.bam"])]
        pBuildIndex=sp.Popen(cmd,stderr=errorlog)
        retVal=pBuildIndex.wait()
        if retVal!=0:
            raise RunError(" ".join(cmd))

    def bwamem_picard_GATK(self):
        """ Use bwa MEM and GATK/picard for reference assembly """
        if nonexist(".".join([self.aln,"vcf"])) or not self.A_OK:
                if nonexist(".".join([self.aln,"sorted.dedup.bam"])) or not self.A_OK:
                    if nonexist(".".join([self.aln,"sorted.bam"])) or not self.A_OK:
                        self.runbwamem()
                    self.deduplicate()
                    if args.clean:
                        os.unlink(".".join([self.aln,"sorted.bam"]))
                        os.unlink(".".join([self.aln,"sorted.bai"]))
                self.callSNPs1()

    def bwamem_elprep_GATK(self):
        """ Use bwa MEM, elprep and GATK for reference assembly """
        if nonexist(".".join([self.aln,"vcf"])) or not self.A_OK:
                if nonexist(".".join([self.aln,"sorted.dedup.bam"])) or not self.A_OK:
                        self.runbwamemelprep()
                self.callSNPs1()

    def runbwa(self):
## Align reads to reference genome
        if(args.long):
            pBwaSW=sp.Popen([bwapath,"bwasw","-t",str(self.parent.factor),self.parent.path,".".join([self.basename+"_1","fastq"]),self.parent.path,".".join([self.basename+"_2","fastq"])],stdout=sp.PIPE,stderr=devnull)
            pView=sp.Popen([samtoolspath,"view","-Sb",
                            "-"],
                           stdin=pBwaSW.stdout,
                           stdout=sp.PIPE,
                           stderr=open("/dev/null",'w'))
            
## Sort reads in BAM-file - [aln].sorted is a prefix for the output file.
            pSort=sp.Popen([samtoolspath,"sort","-",".".join([self.aln,"sorted"])],stdin=pView.stdout)
            pSort.wait()
        else:
            runStdout([bwapath,"aln","-t",str(self.parent.factor),self.parent.path,".".join([self.basename+"_1","fastq"])],".".join([self.basename+"_1","sai"]))
            runStdout([bwapath,"aln","-t",str(self.parent.factor),self.parent.path,".".join([self.basename+"_2","fastq"])],".".join([self.basename+"_2","sai"]))
            if nonexist(".".join([self.aln,"sorted","bam"])):
                pSampe=sp.Popen([bwapath,"sampe","-A",self.parent.path,
                                 ".".join([self.basename+"_1","sai"]),
                                 ".".join([self.basename+"_2","sai"]),
                                 ".".join([self.basename+"_1","fastq"]),
                                 ".".join([self.basename+"_2","fastq"])],
                                stdout=sp.PIPE,stderr=devnull)

                pView=sp.Popen([samtoolspath,"view","-Sb",
                                "-"],
                               stdin=pSampe.stdout,
                               stdout=sp.PIPE,
                               stderr=devnull)
                
## Sort reads in BAM-file - [aln].sorted is a prefix for the output file.
                pSort=sp.Popen([samtoolspath,"sort","-",
                                ".".join([self.aln,"sorted"])],
                               stdin=pView.stdout)
                pSort.wait()
## Index the sorted BAM-file
        run([samtoolspath,"index",
             ".".join([self.aln,"sorted","bam"])],
            ".".join([self.aln,"sorted","bam","bai"]))
## Make pileup to check coverage and sequence quality
        if nonexist(os.path.join(self.dir,".".join([self.aln,"pileup","gz"]))):
            pPileup=sp.Popen([samtoolspath,"mpileup","-f",self.parent.path,".".join([self.aln,"sorted","bam"])],stdout=sp.PIPE)
            pGz=sp.Popen(["gzip"],stdin=pPileup.stdout,stdout=open(os.path.join(self.dir,".".join([self.aln,"pileup","gz"])),'w'))
## Call snps from alignment
        if nonexist(".".join([self.aln,"vcf"])) and not self.A_OK:
            self.A_OK=False
            pPileup=sp.Popen([samtoolspath,"mpileup",
                              "-gf",self.parent.path,
                              ".".join([self.aln,"sorted","bam"])],stdout=sp.PIPE)
            pView=sp.Popen([bcftoolspath,"view",
                            "-cg","-"],
                           stdin=pPileup.stdout,
                           stdout=open(".".join([self.aln,"vcf"]),'w'))
            pView.wait()
## Generate FastA file of calls
        runStdout([vcfutilspath,
                   "vcf2fq",
                   ".".join([self.aln,"vcf"])],
                  os.path.join(self.dir,".".join([self.aln,"fq"])))
        runStdout(["fq2fsa.py",
                   os.path.join(self.dir,".".join([self.aln,"fq"]))],
                  os.path.join(self.dir,".".join([self.aln,"fasta"])))
    def edena(self):
        run([edenapath,"-nThreads",str(args.cpu),
             "-DRpairs",".".join([self.basename+"_1","fastq"]),
             ".".join([self.basename+"_2","fastq"]),
             "-t",str(args.minlength),
             "-p",os.path.join(self.dir,".".join([self.basename,edenapath]))],
            os.path.join(self.dir,".".join([self.basename,edenapath,"ovl"])))
        run([edenapath,
             "-e",os.path.join(self.dir,".".join([self.basename,edenapath,"ovl"])),
             "-t",str(20),
             "-p",os.path.join(self.dir,".".join([self.basename,edenapath]))],
            os.path.join(self.dir,".".join([self.basename,"edena_contigs","fasta"])))
        runStdout(["rename.pl",os.path.join(self.dir,".".join([self.basename,"edena_contigs","fasta"]))],"tmp.fasta")
        run(["mv","tmp.fasta",os.path.join(self.dir,".".join([self.basename,"edena_contigs","fasta"]))],"xxx")
        #run(["formatdb","-p","F","-i", os.path.join(self.dir,".".join([self.basename,"edena_contigs","fasta"]))],"xxx")
        #runStdout(["blastall","-p","blastn","-F","F","-m","9",
        #           "-d",os.path.join(self.dir,".".join([self.basename,"edena_contigs","fasta"])),
        #           "-i",os.path.join(self.dir,".".join([self.basename,"edena_contigs","fasta"]))],os.path.join(self.dir,".".join([self.basename,"edena_contigs","fasta","blast"])))
    def soap(self):
        conf="""max_rd_len=150
[LIB]
avg_ins=400
reverse_seq=0
asm_flags=3
#rd_len_cutoff=200
#rank=1
#pair_num_cutoff=22
#map_len=100
"""
        conf += "q1=" + ".".join([self.basename+"_1","fastq"]) + "\n"
        conf += "q2=" + ".".join([self.basename+"_2","fastq"]) + "\n"
        of= open(".".join([self.basename,"soap","conf"]),'w')
        print(conf,file=of)
        of.close()
        run([soappath,'all',
             '-K',"31",'-p',str(args.cpu),'-d',"4",'-D',"4",
             '-s',".".join([self.basename,"soap","conf"]),
             '-o',os.path.join(self.dir,self.basename)],os.path.join(self.dir,".".join([self.basename,"contigs"])))
    def velvet(self):
        run(["".join(velvetpath,'h'),self.dir,"31",'-create_binary','-fastq',".".join([self.basename+"_1","fastq"]),".".join([self.basename+"_2","fastq"])],"xxx")
        run(["".join(velvetpath,'g'),self.dir,'-cov_cutoff',"auto",'-exp_cov',"auto"],"xxx")
    def SEQuel(self):
        print("Running prep on %s" % (self.basename),file=sys.stderr)
        run(["prep.pl",
             "-r1",".".join([self.basename+"_1","fastq"]),
             "-r2",".".join([self.basename+"_2","fastq"]),
             "-c",os.path.join(self.dir,".".join([self.basename,
                                                  "edena_contigs",
                                                  "fasta"])),
             "-o",os.path.join(self.dir,"prep"),
             "-t",str(args.cpu)],
            os.path.join(self.dir,"prep"))
        run(["rm","-f",os.path.join(self.dir,"prep","reads_aln.sam")],"xxx")

    def getsnps(self):
        run(["dnadiff",self.parent.path,os.path.join(self.dir,".".join([self.aln,"fasta"])),"-p",os.path.join(self.dir,self.aln)],os.path.join(self.dir,".".join([self.aln,"snps"])))

    def getGATKsnps(self):
        print("Calling SNPs round 2...")
        #if nonexist(".".join([self.aln,"snps","vcf"])) :
        if args.haploid:
        #Use haploid model:
            cmd=[javapath,"-jar",GATKpath,"-T","UnifiedGenotyper","-nt","1","-I",".".join([self.aln,"sorted.dedup.bam"]),"-R",self.parent.path,"-o",".".join([self.aln,"snps","vcf","tmp"]),"-ploidy","1","--genotype_likelihoods_model","SNP","--annotateNDA","-maxAltAlleles","3","--genotyping_mode","GENOTYPE_GIVEN_ALLELES","--alleles","snps.vcf","--output_mode","EMIT_ALL_SITES"]
            pGenotyper=sp.Popen(cmd,stderr=errorlog,stdout=devnull)
        else:
            #Use diploid model:
            cmd=[javapath,"-jar",GATKpath,"-T","UnifiedGenotyper","-nt","1","-I",".".join([self.aln,"sorted.dedup.bam"]),"-R",self.parent.path,"-o",".".join([self.aln,"snps","vcf","tmp"]),"-ploidy","2","--genotype_likelihoods_model","SNP","--annotateNDA","-maxAltAlleles","6","--genotyping_mode","GENOTYPE_GIVEN_ALLELES","--alleles",".".join([self.parent.basename,"snps","vcf"]),"--output_mode","EMIT_ALL_SITES"]
            pGenotyper=sp.Popen(cmd,stderr=errorlog,stdout=devnull)
                
        retVal=pGenotyper.wait()
        if retVal==0:
            untmp(".".join([self.aln,"snps","vcf"]))
        else:
            raise RunError(" ".join(cmd))

class SRA (Read):
    def __init__(self,parent,url):
        self.pairs=True
        self.url=list(map(os.path.abspath,[url]))
        self.parent=parent
        self.basename=os.path.splitext(os.path.basename(self.url[0]))[0]
        self.aln="_".join([self.parent.basename,self.basename])
        self.dir=os.path.join(self.parent.dir,self.basename)
        if nonexist(self.dir):
            os.mkdir(self.dir)
        self.path=os.path.join(
            self.dir,
            ".".join([self.basename,"sra"]))
    def getData(self):
        #wget(self.url,self.path)
        ## Convert NCBI SRA-file to fastQ format
        ## Extract reads to $tmpdir/$reads.fastq
        if args.offset:
            pass
            #run(["fastq-dump", "--split-spot","--skip-technical","--gzip","--Q",str(args.offset),self.path],
             #   ".".join([self.basename,"fastq","gz"]))
        else:
            pass
            #run(["fastq-dump", "--split-spot","--skip-technical","--gzip",self.path],
             #   ".".join([self.basename,"fastq","gz"]))
        ## Split reads into files $tmpdir/$read_[12].fastq
        cmd=["rndreads.py",
             "-T",str(args.trim),
             "-m",str(args.minlength),
             "-p",self.basename,".".join([self.basename,"fastq",".gz"])],
        if args.verbose:
            cmd.append("-v")
        if args.qtrim:
            cmd.append("-q")
            cmd.append(str(args.qtrim))
        if args.gzip:
            cmd+= ["-g"]
            self.A_OK=run(cmd,".".join([self.basename+"_1","fastq.gz"]))
        else:
            self.A_OK=run(cmd,".".join([self.basename+"_1","fastq"]))

# class TGZ (Read):
#     def __init__(self,parent,url):
#         self.url=list(map(os.path.abspath,url))
#         self.parent=parent
#         self.basename=os.path.splitext(os.path.basename(url))[0]
#         self.aln="_".join([self.parent.basename,self.basename])
#         self.dir=os.path.join(self.parent.dir,self.basename)
#         if nonexist(self.dir):
#             os.mkdir(self.dir)
#         self.path=os.path.join(
#             self.dir,
#             ".".join([self.basename,"sra"]))
#     def getData(self):
#         wget(self.url,self.path)
#         ## Extract reads from .tar.gz file
#         ## Extract reads to $tmpdir/$reads.fastq
#         os.chdir(self.parent.parent.workdir)
#         run(["fastq-dump", "--split-spot",self.path],
#             ".".join([self.basename,"fastq"]))
#         ## Split reads into files $tmpdir/$read_[12].fastq
#         run(["rndreads.py",
#              "-p",self.basename,".".join([self.basename,"fastq"])],
#             ".".join([self.basename+"_1","fastq"]))
#         os.chdir(self.parent.parent.curdir)

class PE (Read):
    def __init__(self,parent,urls):
        self.files=list()
        for url in urls:
            self.files.append(SeqUtils.Filename(url))
        self.pairs=True
        self.urls=list(map(os.path.abspath,urls))
        self.parent=parent
        self.basename=self.files[0].ID
        self.aln="_".join([self.parent.basename,self.basename])
        self.dir=os.path.join(self.parent.dir,self.basename)
        if nonexist(self.dir):
            os.mkdir(self.dir)
## Retrieve data
    def getData(self):
        cmd=["rndreads.py",
                 "-T",str(args.trim),
                 "-p",self.basename,
                 "-m",str(args.minlength)]+self.urls
        if args.offset:
            cmd += ["-Q",str(args.offset)]
        if args.qtrim:
            cmd += ["-q",str(args.qtrim)]
        if args.gzip:
            cmd+= ["-g"]
            self.A_OK=run(cmd,".".join([self.basename+"_1","fastq.gz"]))
        else:
            self.A_OK=run(cmd,".".join([self.basename+"_1","fastq"]))

class SE (Read):
    def __init__(self,parent,urls):
        self.pairs=False
        self.parent=parent
        self.files=SeqUtils.Filename(urls[0])
        self.basename=self.files.ID
        self.aln="_".join([self.parent.basename,self.basename])
        self.dir=os.path.join(self.parent.dir,self.basename)
        if nonexist(self.dir):
            os.mkdir(self.dir)
    def getData(self):
        readsfile=os.path
        cmd=["rndreads.py","-s",
                 "-T",str(args.trim),
                 "-p",self.basename,
                 "-m",str(args.minlength)]+[self.files.abspath]
        if args.offset:
            cmd += ["-Q",str(args.offset)]
        if args.qtrim:
            cmd += ["-q",str(args.qtrim)]
        self.A_OK=runStdout(cmd,".".join([self.basename+"_1","fastq"]))
        

class Reference:
    """ Holds methods for per reference genome analysis, which usually call pr read analysis for the associated reads. """
    def __init__(self,parent,url):
        self.parent=parent
        self.url=os.path.abspath(url)
        self.reads=list()
        self.basename=os.path.splitext(os.path.basename(self.url))[0]#was ref
        self.dir=os.path.abspath(self.basename) # was refdir
        if nonexist(self.dir):
            os.mkdir(self.dir)
        self.path=os.path.join(self.dir,".".join([self.basename,"fasta"]))

    # def addTGZ(self,url):
    #     self.reads.append()
    def addSRA(self,url):
        self.reads.append(SRA(self,url))
    def addPE(self,urls):
        self.reads.append(PE(self,urls))
    def addSE(self,urls):
        self.reads.append(SE(self,urls))
    def getData(self):
        #if not self.denovo:
        wget(self.url,self.path)
        os.chdir(self.parent.workdir)
        thrPool = ThreadPool(args.cpu)
        thrPool.map(lambda x: x.getData(),self.reads)
        os.chdir(self.parent.curdir)
    def index(self):
        ''' Index reference genome '''
        if args.long:
            run([bwapath,"index","-a","bwtsw",self.path],
            "xxx")
        else:
            run([bwapath,"index",self.path],
                os.path.join(self.dir,".".join([self.basename,"fasta","fai"])))
        run([samtoolspath,"faidx",self.path],"xxx")
        run([javapath,"-jar",os.path.join(picardpath,"CreateSequenceDictionary.jar"),"R=",self.path,"O=",self.path.replace(".fasta",".dict")],self.path.replace(".fasta",".dict"))
    def refalign(self):
        self.factor=2
        thrPool = ThreadPool(max(int(args.cpu/self.factor),1))
        if args.reallyold:
            thrPool.map(lambda x: x.runbwa(),self.reads)
        elif args.old:
            thrPool.map(lambda x: x.bwamem_picard_GATK(),self.reads)
        else:
            thrPool.map(lambda x: x.bwamem_elprep_GATK(),self.reads)

    def mergeSNPs(self):
        """Merges the SNPs called between reference and isolates to a single set."""
        runStdout(["egrep","^#",".".join([self.reads[0].aln,"vcf"])],".".join([self.basename,"snps","header"]))
        contigre=re.compile('##contig=<ID=(.+),length')
        contigs=list()
        ### Make sure we add contigs in the expected order by reading from the header file
        for line in open(".".join([self.basename,"snps","header"]),'r'):
            m=contigre.match(line)
            if(m):
                contigs.append(m.group(1).strip('"'))
        run(["rm","-f",".".join([self.basename,"snps","vcf","tmp"])],"xxx")
        for contig in contigs:
            contig=re.escape(contig)
            cmd="grep -P -h '^" + contig + "\s'"
            for read in self.reads:
                cmd += " \"" + ".".join([read.aln,"vcf"]) + "\""
            cmd += " | sort -k2,2n | gawk 'BEGIN{oldline=\"\"}{if($2 != oldline){print $0;oldline=$2}}' >> \""  + ".".join([self.basename,"snps","vcf","tmp"]) + "\""
            os.system(cmd)

        fo=open(".".join([self.basename,"snps","vcf"]),"w")
        sp.check_call(["cat",".".join([self.basename,"snps","header"]),".".join([self.basename,"snps","vcf","tmp"])], stdout=fo)
        fo.close()
        time.sleep(3)
        
    def recallSNPs(self):
        """Re-evaluate the SNPs in all potential SNP positions"""
        thrPool = ThreadPool(args.cpu)
        thrPool.map(lambda x: x.getGATKsnps(),self.reads)

    def makeSNPmat(self):
        """Collects the SNPs called in recallSNPs and saves it in various useful formats"""
        cmd = ["vcf2snpmat.py"] if args.ambiguous else ["vcf2snpmat.py","-c"]
        cmd.extend(["-s",".snps.vcf","-p","{}_".format(self.basename)])
        if args.mindepth:
            cmd.extend(["-d",str(args.mindepth)])
        for read in self.reads:
            cmd.append(".".join([read.aln,"snps","vcf"]))
        fo=open(".".join([self.basename,"snpmat"]),"w")
        fe=open(".".join([self.basename,"discarded_snps"]),"w")
        p=sp.check_call(cmd,stdout=fo,stderr=fe)
        fe.close()
        fo.close()
        fo=open(".".join([self.basename,"phylip"]),"w")
        p=sp.check_call(["snpmat2phylip.py",'-p',".".join([self.basename,"positions"]),'-s',".".join([self.basename,"rename","sed"]),".".join([self.basename,"snpmat"])],stdout=fo)
        fo.close()
        cmdP2F=["phylip2fasta.py",".".join([self.basename,"phylip"])]
        pP2F=sp.Popen(cmdP2F,
                      stdout=sp.PIPE)
        cmdRename=["sed","-f",".".join([self.basename,"rename","sed"])]
        pRename=sp.Popen(cmdRename,
                         stdin=pP2F.stdout,
                         stdout=open(".".join([self.basename,"snps","fsa","tmp"]),"w"))
        retVal=pRename.wait()
        if retVal==0:
            untmp(".".join([self.basename,"snps","fsa"]))
        else:
            raise RunError(" | ".join([" ".join(cmdP2F)," ".join(cmdRename)]))

        """ Clean recombination"""
        fo=open(".".join([self.basename,"clean","phylip"]),"w")
        fe=open(".".join([self.basename,"HTpositions"]),"w")
        p=sp.check_call(["cleanrecomb2.py",'-p',".".join([self.basename,"positions"]),".".join([self.basename,"phylip"])],stdout=fo,stderr=fe)
        fo.close()
        fe.close()
        cmdP2F=["phylip2fasta.py",".".join([self.basename,"clean","phylip"])]
        pP2F=sp.Popen(cmdP2F,
                      stdout=sp.PIPE)
        cmdRename=["sed","-f",".".join([self.basename,"rename","sed"])]
        pRename=sp.Popen(cmdRename,
                         stdin=pP2F.stdout,
                         stdout=open(".".join([self.basename,"clean","snps","fsa","tmp"]),"w"))
        pRename.wait()
        if retVal==0:
            untmp(".".join([self.basename,"clean","snps","fsa"]))
        else:
            raise RunError(" | ".join([" ".join(cmdP2F)," ".join(cmdRename)]))


    def snptree(self):
        """ Construct an MP tree """
        self.mergeSNPs()
        self.recallSNPs()
        self.makeSNPmat()
        print("Making tree")
        sp.check_call(["cp",".".join([self.basename,"clean","phylip"]),"infile"])
        run(["rm","-f","outfile","outtree"],"xxx")
        pEchoY= sp.Popen(["echo","Y"],stdout=sp.PIPE)
        pDnapars = sp.Popen(["phylip","dnapars"],stdin=pEchoY.stdout,stdout=devnull)
        retval=pDnapars.wait()
        if retval==0:
            sp.check_call(["sed","-f",".".join([self.basename,"rename","sed"]),"outfile"],
                          stdout=open(".".join([self.basename,"clean","dnapars"]),'w'))
            sp.check_call(["sed","-f",".".join([self.basename,"rename","sed"]),"outtree"],
                          stdout=open(".".join([self.basename,"clean","nwk"]),'w'))
        else:
            raiseRunError("echo Y | phylip dnapars")
        print("Done")
        
    def denovo(self):
        for read in self.reads:
            try:
                func=getattr(read,args.d)
            except AttributeError:
                print("Unknown denovo type: {}".format(args.d),file=sys.stderr)
            func()
    def SEQuel(self):
        for read in self.reads:
            read.SEQuel()

class Project:
    def __init__(self,args):
        if(args.workdir):
            self.workdir=args.workdir
            if not os.path.exists(self.workdir):
                mkdir=sp.Popen(["mkdir",self.workdir])
                mkdir.wait()
        else:
            self.workdir=tempfile.mkdtemp(dir=".")
        self.curdir=os.getcwd()
        fh=open(args.config)
        endl=fh.newlines
        self.reference=None
        for reads in fh:
            if reads.startswith('#'):
                continue
            if reads == endl:
                continue
            if not self.reference:
                self.reference=Reference(self,reads.strip(endl))
                continue
            if(reads.startswith("SRA=")):
                self.reference.addSRA(reads[4:].strip(endl))
            # if(reads.startswith("TGZ=")):
            #     self.reference.addTGZ(reads[4:].strip(endl))
            elif(reads.startswith("SE=")):
                self.reference.addSE(reads[3:].strip(endl))
            else:
                urls = list()
                for url in reads.strip(endl).split("\t"):
                    urls.extend(glob.glob(url))
                if len(urls)<3 and len(urls)>0:
                    self.reference.addPE(urls)
                else:
                    print("Didn't find two reads:\n{}".format(", ".join(urls)),file=sys.stderr)
        fh.close()

        self.reference.getData()
        if not args.denovo:
            self.reference.index()

    def refalign(self):
        os.chdir(self.workdir)
        self.reference.refalign()
        os.chdir(self.curdir)

    def snptree(self):
        os.chdir(self.workdir)
        self.reference.snptree()
        os.chdir(self.curdir)
        
    def denovo(self):
        os.chdir(self.workdir)
        self.reference.denovo()
        os.chdir(self.curdir)

    def SEQuel(self):
        os.chdir(self.workdir)
        self.reference.SEQuel()
        os.chdir(self.curdir)
        

def wget(url,destination):
    """ Fetches remote files """
    print("Fetching: %s" % (destination))
    if os.path.splitext(url)[1]==".gz":
        if nonexist(destination):
            run(["cp",url,destination+".gz"],destination+".gz")
            #urllib.request.urlretrieve(url,filename=destination+".gz")
            run(["gunzip",destination+".gz"],destination)
    else:
        if nonexist(destination):
            run(["cp",url,destination],destination)
            #urllib.request.urlretrieve(url,filename=destination)

def nonexist(file=""):
    """Checks for file nonexistance"""
    return(not os.access(file,os.F_OK))

def runStdout(cmd,outfile):
    """Run command and redirect output to outfile"""
    if nonexist(outfile):
        fh=open(outfile,"w")
        if args.verbose:
            sp.check_call(cmd,stdout=fh)
        else:
            sp.check_call(cmd,stdout=fh,stderr=open("/dev/null",'w'))
        fh.close()
        return False
    else:
        return True

def run(cmd,outfile):
    """Run command and check if outfile is nonexistant"""
    if nonexist(outfile):
        if args.verbose:
            sp.check_call(cmd)
        else:
            sp.check_call(cmd,stderr=open("/dev/null",'w'))
        return False
    else:
        return True

def rmwalk(dir):
    """ Remove directory """
    try:
        for (dirpath,dirlist,filelist) in os.walk(dir,topdown=False):
            for f in filelist:
                os.unlink(os.path.join(dirpath,f))
            os.rmdir(dirpath)

    except ValueError:
        print("Failed to delete %s" % (dir),file=sys.stderr)
        return


parser=argparse.ArgumentParser(description='NGS pipeline for reference assembly')
parser.add_argument("config", help="Config file")
parser.add_argument("--denovo",help="Run only denovo assembly",action="store_true")
parser.add_argument("-d",type=str,help="Type of denovo analysis",default="edena")
parser.add_argument("-m","--minlength",help="Discard reads shorter than readlength",type=int,default=145)
parser.add_argument("--mindepth",help="Minimum depth to call snps",type=int,default=10)
parser.add_argument("-f","--fraction",help="Fraction of reads to use",type=float,default=1.0)
parser.add_argument("-w","--workdir",help="Give a specific working directory")
parser.add_argument("-t","--cpu",help="Number of concurrent threads to allow",type=int,default=3)
parser.add_argument("-T","--trim",help="Number of initial bases to trim from reads",type=int,default=0)
parser.add_argument("-c","--clean",help="Clean up work directory",action="store_true")
parser.add_argument("-v","--verbose",help="Be more verbal about what is happening",action="store_true")
parser.add_argument("-Q","--offset",help="Quality score offset -- used to override the automatic quality detection",type=int,default=0)
parser.add_argument("-q","--qtrim",help="Trim ends by quality",type=int)
parser.add_argument("-g","--gzip",help="Gzip fastq files",action="store_true")
parser.add_argument("-a","--ambiguous",help="Allow ambiguous base calls as SNPS",action="store_true")
parser.add_argument("-l","--long",help="Long reads alignment -- use only with legacy bwa",action="store_true")
parser.add_argument("--reallyold",help="Run legacy bwa",action="store_true")
parser.add_argument("--old",help="Use picard and samtools instead of elprep",action="store_true")
parser.add_argument("-1","--haploid",help="Use GATK haploid model",action="store_true")
parser.add_argument("--sequel",help="Runs sequel on edena contigs",action="store_true")
args =parser.parse_args()

## Use only 50% of reads to save memory on desktop computer
frac=args.fraction


if len(sys.argv)<2:
    print("Usage: ./pipeline.py <configfile>\n")
    print("The config file starts with a reference file on the first line, and read files on subsequent lines\n")
    print("Example:\nNC_017178.fasta\nSRA=ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/ERX/ERX102/ERX102103/ERR125954/ERR125954.sra\n~/data/reads/SRR010965.sra\n")
    sys.exit()



prj=Project(args)
#errorlog=sys.stderr # Use to log output
errorlog=open(os.path.join(prj.workdir,"pipeline.log"),"w")
if( not args.denovo):
    prj.refalign()
    prj.snptree()
else:
    prj.denovo()
    if args.sequel:
        prj.SEQuel()
#if args.clean:
#    rmwalk(self.workdir)



