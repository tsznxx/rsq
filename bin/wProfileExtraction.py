#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Last-modified: 24 Mar 2015 03:35:47 PM
#
#         Module/Scripts Description
# 
# Copyright (c) 2014 Yunfei Wang <yfwang0405@gmail.com>
#
#   __     __           __     _  __          __               
#   \ \   / /          / _|   (_) \ \        / /               
#    \ \_/ /   _ _ __ | |_ ___ _   \ \  /\  / /_ _ _ __   __ _ 
#     \   / | | | '_ \|  _/ _ \ |   \ \/  \/ / _` | '_ \ / _` |
#      | || |_| | | | | ||  __/ |    \  /\  / (_| | | | | (_| |
#      |_| \__,_|_| |_|_| \___|_|     \/  \/ \__,_|_| |_|\__, |
#                                                         __/ |
#                                                        |___/
#   
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @status:  experimental
# @version: 1.0.0
# @author:  Yunfei Wang
# @contact: yfwang0405@gmail.com

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import time
import numpy
import pysam
import argparse
import itertools
from bisect import bisect_left,bisect_right
from ngslib import IO,BigWigFile,Utils

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Convert BAM file to Wiggle file. This is for PARS data analysis. The BAM file is sorted and indexed. See http://rsqwiki.appspot.com/mapping for details. Contact Yunfei Wang to report any bugs (yfwang0405@gmail.com).',epilog='dependency ngslib')
    p._optionals.title = "Options"
    p.add_argument("-i","--input",dest="input",type=str,metavar="GM12878_S1.bam", required=True, nargs='+',help="Sorted and indexed BAM file(s).")
    p.add_argument("-a","--anno",dest="anno",type=str,metavar="hg19_RefSeq.genepred",required=True,help="Gene annotation file in genepred format.")
    p.add_argument("-c","--chrsizes",dest="chrsizes",type=str,metavar='hg19.sizes',required=True,help="Chrom sizes file. A file with chromosome names and sizes in the first two columns.")
    p.add_argument("-s","--shift",dest="shift",type=int,metavar='-1',default=-1,help="Read shift to enzyme recognition site. [default=-1 for PARS technology]")
    p.add_argument("-o","--outdir",dest="outdir",type=str,metavar="outdir",default=".",help="Directory to put BigWig file(s). [default= \".\"].")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

def groupGeneByChroms(anno,chroms):
    '''
    Group transcripts by chromosome names.
    '''
    genes = {}
    gene_in_chroms = {chrom:[] for chrom in chroms}
    for gene in IO.BioReader(anno,ftype='genepred'):
        offsets = [0]
        starts = []
        for exon in gene.exons():
            starts.append(exon.start)
            offsets.append(offsets[-1]+exon.length())
        genes[gene.id] = (offsets,starts,gene.strand)
        gene_in_chroms[gene.chrom].append(gene.id)
    return genes,gene_in_chroms

def BAMToWig(bfile,chroms,genes,gene_in_chroms,shift,outdir):
    ''' 
    Convert BAM file to Wiggle then BigWig file.
    '''
    sname=os.path.basename(bfile)
    sname=os.path.splitext(sname)[0]
    sam = pysam.Samfile(bfile,'rb')
    bwfs = {'+':open(outdir+"/"+sname+"_plus.wig",'w'),'-':open(outdir+"/"+sname+"_minus.wig",'w')}
    for strand in bwfs:
        bwfs[strand].write("track type=wiggle_0 name={0}_{1}\n".format(sname,strand == '+' and 'plus' or 'minus'))
    maxchrlen = max(chroms.values()) + abs(shift)
    depth = {'+':numpy.zeros(maxchrlen,dtype=numpy.uint16),'-':numpy.zeros(maxchrlen,dtype=numpy.uint16)}
    for chrom in sorted(chroms):
        chrlen = chroms[chrom]
        # fetch reads mapped to chrom
        for read in sam.fetch(chrom):
            if not read.is_paired or read.is_read1: # only consider the first read.
                strand = read.is_reverse and "-" or "+"
                pos = read.pos + (-shift if read.is_reverse else shift)
                depth[strand][pos] += 1
        for gid in gene_in_chroms[chrom]:
            offsets,starts,strand = genes[gid]
            try:
                for read in sam.fetch(gid):
                    if not read.is_paired or read.is_read1:
                        idx = bisect_right(offsets,read.pos)-1
                        pos = starts[idx]+read.pos-offsets[idx]
                        gstrand = (strand == '+')^read.is_reverse and '+' or '-'# gstrand = '+' if strand == read.strand else '-'
                        pos += gstrand=='+' and shift or -shift # plus strand -1 , minus strand +1
                        depth[gstrand][pos] += 1
            except:
                pass # genes have no read mapped to
        # report wig file
        for strand,bwf in bwfs.iteritems():
            # depth[strand][chrlen:] = 0 # set values to zero if exceed the chrom length range.
            bwf.write("variableStep chrom={0}\n".format(chrom))
            for i in numpy.nonzero(depth[strand][:chrlen])[0]:
                bwf.write("{0}\t{1}\n".format(i+1,depth[strand][i]))
                depth[strand][i] = 0 # reset to zero
    sam.close()
    for bwf in bwfs.values():
        bwf.close()
    return sname

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    args=argParser()

    # Read chrom sizes
    chroms = Utils.genomeSize(args.chrsizes)
    # Read gene annotation
    print >>sys.stderr, "Grouping transcripts by chromosome names:", args.anno
    timestart = time.time()
    genes,gene_in_chroms = groupGeneByChroms(args.anno,chroms)
    print >>sys.stderr, "Finished in {0}s.\n".format(time.time()-timestart)
        
    # read BAM file
    print >>sys.stderr, "BAM->Wiggle conversion ..."
    for bfile in args.input:
        print >> sys.stderr, "Processing BAM file: {0}".format(bfile),
        timestart = time.time()
        sname = BAMToWig(bfile,chroms,genes,gene_in_chroms,args.shift,args.outdir)
        print >> sys.stderr, "Wiggle file created in {0}s.".format(time.time()-timestart)
        # Wiggle -> BigWig
        print >>sys.stderr, "Covert Wiggle files to BigWig files ..."
        for ext in ['_plus.wig','_minus.wig']:
            timestart = time.time()
            bwfname = sname +ext[:-3]+'bw'
            BigWigFile.wigToBigWig(sname+ext,args.chrsizes,bwfname)
            print >>sys.stderr, "Finished {0} -> {1} in {2}s.".format(fname,bwfname,time.time()-timestart)

    
