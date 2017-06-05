#!/bin/sh
#Last-modified: 02 Mar 2015 11:31:53 AM

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


USAGE=" Usage: $0 demo\n\tFor demo example, run\n\t\t>$0 demo\n\tTo clean demo results, run\n\t\t>$0 clean\n\tPut in real data for your own case.\n"

case $# in
    0)  echo -en $USAGE
        exit;;
    *)  if [ "$1" == "demo" ]; then
			prefix="demo"
        else if [ "$1" == "clean" ]; then
                rm -rf workdir *merged.fd *merged.fc *merged.fs *fitness
                exit
             else
                 echo -en $USAGE
                 exit
             fi
        fi
        ;;
esac



# Dataset:
# FastD: demo01.fd demo02.fd
# FastC: demo01.fc demo02.fc
# FastS: demo01.fs demo02.fs known.fs

# 1. Calculate fitness scores given known strucutres.
echo -e "\n###################################\nStep 1: FastD fitness to know structures."
rsq fitness -d ${prefix}01.fd -s known.fs -o ${prefix}01.fitness
rsq fitness -d ${prefix}02.fd -s known.fs -o ${prefix}02.fitness

# 2. Merge FastD with weights from fitness analysis.
echo -e "\n###################################\nStep 2: Merge FastD based on fitness score to known structures."
w1=$(awk 'BEGIN{FS=OFS="\t";cnt=0;score=0}{if(NR != 1) {cnt+=1;score+=$4}}END{print score/cnt}' ${prefix}01.fitness)
w2=$(awk 'BEGIN{FS=OFS="\t";cnt=0;score=0}{if(NR != 1) {cnt+=1;score+=$4}}END{print score/cnt}' ${prefix}02.fitness)
rsq FastD merge -i ${prefix}01.fd ${prefix}02.fd -w $w1 $w2 -o ${prefix}_merged.fd

# 3. Merge FastC files
echo -e "\n###################################\nStep 3: Merge FastC files by intesecting constraints."
rsq FastC merge -i ${prefix}01.fc ${prefix}02.fc -m i -o ${prefix}_merged.fc

# 4. Merge FastS files
echo -e "\n###################################\nStep 4: Merge FastS files."
rsq FastS merge -i ${prefix}01.fs ${prefix}02.fs -o ${prefix}_merged.fs

