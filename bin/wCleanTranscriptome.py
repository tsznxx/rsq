#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Last-modified: 20 Mar 2015 03:57:25 PM
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

import os,sys
from ngslib import IO,BedList

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        sys.exit("Example: python "+sys.argv[0]+" hg19_RefSeq.genepred >hg19_RefSeq_relabel.genepred")
    genes = BedList()
    gids = {}
    # Reading transcriptome annotation file.
    for gene in IO.BioReader(sys.argv[1],ftype='genepred'):
        if not '_' in gene.chrom:
            genes.append(gene)
            if gids.has_key(gene.id):
                gids[gene.id][1] += 1
            else:
                gids[gene.id] = [1,1]
    # Relabeling and printing
    genes.sort()
    for gene in genes:
        if gids[gene.id][1] >1:
            gid = gene.id
            gene.id += ":{0}/{1}".format(*gids[gene.id])
            gids[gid][0] += 1
        print gene

