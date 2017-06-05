#!/bin/sh
#Last-modified: 02 Mar 2015 11:31:25 AM

####################### Module/Scripts Description ######################
#  
#  Copyright (c) 2014 Yunfei Wang <yfwang0405@gmail.com>
#  
#  This code is free software; you can redistribute it and/or modify it
#  under the terms of the BSD License (see the file COPYING included with
#  the distribution).
#  
#  @status:  experimental
#  @version: $Revision$
#  @author:  Yunfei Wang
#  @contact: yfwang0405@gmail.com
#
#########################################################################


USAGE=" Usage: $0 sce_genes_demo.fasta sce_S1_demo.tab sce_V1_demo.tab\n\tFor demo example, run\n\t\t>$0 demo\n\tTo clean demo results, run\n\t\t>$0 clean\n\tPut in real data for your own case.\n"

case $# in
    0)  echo -en $USAGE
        exit;;
    1)  if [ "$1" == "demo" ]; then
            FA=sce_genes_demo.fasta
			S=sce_S1_demo.tab
			V=sce_V1_demo.tab
        else if [ "$1" == "clean" ]; then
                rm -rf workdir *fd *fc *fs *isf *txt *sizes *bw *fa *fai *ps
                exit
             else
                 echo -en $USAGE
                 exit
             fi
        fi
        ;;
    2)  echo -en $USAGE
        exit;;
    *)  FA=$1
        S=$2
        V=$3
        ;;
esac

prefix=${FA%%.tab}

# Dataset:
# FastA: sequences of transcripts
# S: RNase S1 depth
# V: RNase V1 depth

# 1. generating FastD file from preprocessed data
echo -e "\n###################################\nStep 1: Generating FastD."
rsq FastD generate -f $FA -S $S -V $V -o ${prefix}.fd

# 2. calculating constraints using 'exclusive' method
echo -e "\n###################################\nStep 2: Generating FastC."
rsq FastC generate -i ${prefix}.fd -m e -o ${prefix}.fc

# 3. fold RNA secondary structure using sfold with 5 CPUs
echo -e "\n###################################\nStep 3: Fold secondary structures."
rsq fold -i ${prefix}.fc -p 5 -o ${prefix}.fs

# 4. Quantification of secondary structures
echo -e "\n###################################\nStep 4: Quantification of structures."
rsq quantify -d ${prefix}.fd -s ${prefix}.fs -o ${prefix}.txt

# 5. Visualization of RNA secondary structures
echo -e "\n###################################\nStep 5: Visualization of structures."
rsq draw -i ${prefix}.fs -s $prefix
