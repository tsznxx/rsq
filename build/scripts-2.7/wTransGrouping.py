#!/net/uu/nm/bi/yxw120430/local/bin/python
# -*- coding: utf-8 -*-
#Last-modified: 25 Mar 2015 02:34:10 PM
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
import rsq
import time
import argparse
from ngslib import IO,BigWigFile,DB,Utils,GeneBed

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
    p.add_argument("-q","--query",dest="qbed",type=str,metavar="query.bed",required=True,help="Query genomic regions in BED format.")
    p.add_argument("-a","--anno",dest="anno",type=str,metavar="hg19_RefSeq.genepred",required=True,help="Gene annotation file in genepred format.")
    p.add_argument("-m","--method",dest="method",type=str,choices=['extend','truncate'],required=False,default='truncate',help="Method to report transcripts overlapped with the given region. 'truncate' will trucate all the overlapped transcripts, while 'extend' will extend the region to cover all the overlaping transcripts.")
    p.add_argument("-S","--bwfS",dest="bwfS",type=str,metavar='prefix',required=False,default=None,help="Prefix of loop (S) BigWig file(s). For example, for GM12878_S1_plus.bw and GM12878_S1_minu.bw, the prefix is 'GM12878_S1'.")
    p.add_argument("-V","--bwfV",dest="bwfV",type=str,metavar='prefix',required=False,default=None,help="Prefix of stem (V) BigWig file(s). For example, for GM12878_V1_plus.bw and GM12878_V1_minu.bw, the prefix is 'GM12878_V1'.")
    p.add_argument("-o","--outfile",dest="outfile",type=str,metavar="outfile",default="stdout",help="Grouped transcripts.")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

def transGroup(bedfile,annofile,bwfS_prefix,bwfV_prefix,method):
    '''
    Fetch transcripts in given BED regions.
    '''
    # Open files
    if bwfS_prefix:
        Utils.mustexist(bwfS_prefix+"_plus.bw")
        Utils.mustexist(bwfS_prefix+"_minus.bw")
        bwfS= {'+':BigWigFile(bwfS_prefix+"_plus.bw"),'-':BigWigFile(bwfS_prefix+"_minus.bw")}
    if bwfV_prefix:
        Utils.mustexist(bwfV_prefix+"_plus.bw")
        Utils.mustexist(bwfV_prefix+"_minus.bw")
        bwfV= {'+':BigWigFile(bwfV_prefix+"_plus.bw"),'-':BigWigFile(bwfV_prefix+"_minus.bw")}
    Utils.mustexist(annofile)
    annodb  = DB(annofile,'genepred')
    # Read Bed file and group transcripts    
    for bed in IO.BioReader(bedfile,'bed'):
        trans = {}
        start = bed.start
        end = bed.stop
        for gene in annodb.fetch(bed.chrom,bed.start,bed.stop,bed.strand,converter=GeneBed):
            trans[gene.id] = gene
            start = min(start,gene.start)
            end   = max(end,gene.stop)
        if method == 'extend': # extend
            cnt = len(trans)
            while True: # extend
                for gene in annodb.fetch(bed.chrom,start,end,bed.strand,converter=GeneBed):
                    if not trans.has_key(gene.id):
                        trans[gene.id] = gene
                        start = min(start,gene.start)
                        end   = max(end,gene.stop)
                if len(trans) > cnt:
                    cnt = len(trans)
                else:
                    break
        elif method == 'truncate': # truncate
            for tr in trans.values():
                if tr.start < bed.start or tr.stop > bed.stop:
                    tr.id += '_truncated'
                starts = []
                stops  = []
                for exon in tr.exons():
                    if bed.overlapLength(exon)>0:
                        starts.append(max(exon.start,bed.start))
                        stops.append(min(exon.stop,bed.stop))
                # update tr
                starts.sort()
                stops.sort()
                tr.exonstarts = starts
                tr.exonstops  = stops
                tr.start = starts[0]
                tr.stop  = stops[-1]
                tr.txstart = max(tr.txstart,tr.start)
                tr.txstop  = min(tr.txstop,tr.stop)
                tr.exoncount = len(tr.exonstarts)
        else:
            raise ValueError("ERROR: Method '{0}' is not recognized. Choices are 'e','extend' or 't','truncate'.".format(method))
        # Report transcripts
        addinfo = (bed.id,start,end,len(trans))
        for tr in sorted(trans.values()):
            tr.otherfields.extend(addinfo)
            yield tr

    if bwfS_prefix:
        for bwf in bwfS.values():
            bwf.close()
    if bwfV_prefix:
        for bwf in bwfV.values():
            bwf.close()
    annodb.close()

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    args = argParser()
    if not ( args.bwfS or args.bwfV):
        sys.exit("ERROR: Either loop (S) or stem (V) bigwig files should be provided.")
    with IO.mopen(args.outfile,'w') as ofh:
        for tr in transGroup(args.qbed,args.anno,args.bwfS,args.bwfV,args.method):
            print >>ofh, tr.allFields()
