#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Last-modified: 11 Dec 2014 03:56:18 PM


        #############################################################################
        #   ███╗   ███╗██╗   ██╗██╗  ████████╗██╗███████╗ ██████╗ ██╗     ██████╗   #
        #   ████╗ ████║██║   ██║██║  ╚══██╔══╝██║██╔════╝██╔═══██╗██║     ██╔══██╗  #
        #   ██╔████╔██║██║   ██║██║     ██║   ██║█████╗  ██║   ██║██║     ██║  ██║  #
        #   ██║╚██╔╝██║██║   ██║██║     ██║   ██║██╔══╝  ██║   ██║██║     ██║  ██║  #
        #   ██║ ╚═╝ ██║╚██████╔╝███████╗██║   ██║██║     ╚██████╔╝███████╗██████╔╝  #
        #   ╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝   ╚═╝╚═╝      ╚═════╝ ╚══════╝╚═════╝   #
        #############################################################################


#         Module/Scripts Description
# 
# Copyright (c) 2014 Yunfei Wang <yfwang0405@gmail.com>
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

# python packages
import os
import re
import sys
import copy
import time
import bisect
import random
import tempfile
from itertools import izip
from subprocess import Popen, PIPE, call
from multiprocessing import Pool,Lock

# external python packages
import numpy
from fisher import pvalue,pvalue_npy

# custom packages
import wRNA
import ngslib

# ------------------------------------
# constants
# ------------------------------------

HOME = os.path.expanduser('~')+"/"
debug = False
rctable = {'A':'U','U':'A','C':'G','G':'C'}

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

class FastA(object):
    '''
    FastA format
        >YIL140W
        AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
    '''
    def __init__(self,name,seq):
        ''' Initiation. '''
        self.name = name
        self.seq = seq.upper()
    def __len__(self):
        ''' length of sequence. '''
        return len(self.seq)
    def __str__(self):
        ''' FastA in string. '''
        return ">{0}\n{1}".format(self.name,self.seq)
    def rc(self):
        '''
        Do reverse complement of self.seq. No returns. 
        '''
        self.seq = ''.join([rctable[c] for c in self.seq])[::-1]

class FastD(FastA):
    '''
    FastD format for storing sequencing depth.
    Example file:
        >tE-UUC-I
        UCCGAUAUAGUGUAACGGCUAUCACAUCACGCUUUCACCGUGGAGACCGGGGUUCGACUCCCCGUAUCGGAG
        27023.0;99.0;16.0;2.0;2.0;56.0;25.0;58.0;27.0;28.0;41.0;51.0;138.0;220.0;37.0;50.0;52.0;12.0;8.0;43.0;23.0;47.0;74.0;34.0;19.0;31.0;17.0;27.0;33.0;36.0;13.0;7.0;19.0;371.0;12.0;146.0;45.0;22.0;31.0;22.0;34.0;46.0;34.0;16.0;4.0;1.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0
        3187.0;13.0;27.0;9.0;0.0;18.0;17.0;15.0;13.0;15.0;16.0;13.0;43.0;29.0;17.0;34.0;21.0;5.0;5.0;12.0;5.0;9.0;25.0;7.0;9.0;13.0;7.0;533.0;55.0;22.0;11.0;4.0;5.0;79.0;2.0;26.0;17.0;18.0;19.0;20.0;15.0;49.0;33.0;20.0;5.0;1.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0
    '''
    def __init__(self, name, seq, loops = None, stems = None):
        ''' Initiation. '''
        super(FastD,self).__init__(name,seq)
        self.seq = self.seq.replace('T','U')
        self.loops = numpy.array(loops,dtype=float) if loops is not None else numpy.zeros(len(seq),dtype=float)
        self.stems = numpy.array(stems,dtype=float) if stems is not None else numpy.zeros(len(seq),dtype=float)
    def __str__(self):
        ''' string of FastD format. '''
        return ">{0}\n{1}\n{2}\n{3}".format(self.name,
                                            self.seq, 
                                            ";".join([str(round(i,2)) for i in self.loops]) if numpy.sum(self.loops)>0 else "+",
                                            ";".join([str(round(i,2)) for i in self.stems]) if numpy.sum(self.stems)>0 else "+") 
    def depth(self,mode = 'b'):
        '''
        Calculate the depth of FastD. 
        Parameters:
            mode: char
                'b': both
                's': stem only, depth 
                'l': loop only
        Note:
            For stem or loop only cases, we assume that the average stem/loop raito is 1. Then the depth for stems and loops are sum(stems)/(length/2) and sum(loops)/(length/2), respectively.
        '''
        if mode == 'b':
            depth = sum(self.loops+self.stems)/len(self.seq)
        elif mode == 's':
            depth = 2.*sum(self.stems)/len(self.seq)
        else:
            depth = 2.*sum(self.loops)/len(self.seq)
        return depth
    def normalize(self, sratio = 1, vratio = 1, ntrim = 5):
        ''' Normalization. '''
        self.loops *= sratio
        self.stems *= vratio
        if ntrim == 0:
            return
        self.loops[0:ntrim] = 0
        self.stems[0:ntrim] = 0
        self.loops[-ntrim:] = 0
        self.stems[-ntrim:] = 0
        return
    def rc(self):
        ''' Reverse complement. '''
        super(FastD,self).rc()
        self.loops = self.loops[::-1]
        self.stems = self.stems[::-1]
    def toFastC(self, method='exclusive',sthreshold = 0.05, vthreshold = 0.05):
        '''
        Calculate constraints from FastD depth. 
        Parameters:
            method: string
                'p':'percentile', top p depths are considered as constraints.
                'f':'fisher', fisher exact test
                'l':'logS/V', log(loop+1)/(stem+1)
                'e':'exclusive', fisher exact test and arbitrary cut on the other data.
                'n':'none', no constraints
                Note: If either S or V data is not provided, and 'n' is not specified, 'p' is used.
            sthreshold: float
                Threshold for loop data to calculate constraints
            vthreshold: float
                Threshold for stem data to calculate constraints
        Returns:
            tfc: FastC object
                Calculated constraints.
        '''
        def pval(d):
            p = pvalue(*d)
            return p.left_tail,p.right_tail
        constraints = numpy.repeat('.', len(self.seq))
        if (sum(self.loops) == 0 or sum(self.stems) == 0) and method not in ('n','none'):
            method = 'p'
            #sthreshold = vthreshold = 0.1
        if method in ('p','percentile'):
            # loops
            if sum(self.loops)>0:
                try:
                    idx = int(len(self.loops)*(1-float(sthreshold)) - 0.00000001)
                    constraints[self.loops > numpy.sort(self.loops)[idx]] = 'x'
                except:
                    print self.name,idx
            # stems
            if sum(self.stems)>0:
                idx = int(len(self.stems)*(1-float(vthreshold)) - 0.00000001)
                constraints[self.stems > numpy.sort(self.stems)[idx]] = '|'
            return FastC(self.name,self.seq,constraints.tostring())
        if method in ['f','fisher']:
            Ssum = sum(self.loops)
            Vsum = sum(self.stems)
            S,V = zip(pvalue_npy(self.loops,Ssum-self.loops,self.stems,Vsum-self.stems)[0:2])
            constraints[numpy.where(S<sthreshold)] = 'x'
            constraints[numpy.where(V<vthreshold)] = '|'
            return FastC(self.name, self.seq, constraints.tostring())
        if method == 'e' or method == 'exclusive':
            # Noise boundary calculation
            nz= sorted(self.loops[numpy.nonzero(self.loops)[0]])
            Sbound = int(round(max(nz[int(len(nz)*0.05+0.5)],5)))
            nz= sorted(self.stems[numpy.nonzero(self.stems)[0]])
            Vbound = int(round(max(nz[int(len(nz)*0.05+0.5)],5)))
            # sum of depth
            Ssum = int(round(sum(self.loops)))
            Vsum = int(round(sum(self.stems)))
            # calcuate boundary values
            fishS = {} # S < V
            for i in range(Sbound+1):
                for j in xrange(Vsum):
                    # p = fisher_exact(((i,Ssum),(j,Vsum)),alternative='less')[1]
                    p = pvalue(i,Ssum,j,Vsum).left_tail
                    if p <= vthreshold:
                        break
                fishS[i] = j
            fishV = {} # S > V
            for i in range(Vbound+1):
                for j in xrange(Ssum):
                    # p = fisher_exact(((i,Vsum),(j,Ssum)),alternative='less')[1]
                    p = pvalue(i,Vsum,j,Ssum).left_tail
                    if p <= sthreshold:
                        break
                fishV[i] = j
            if debug: print >> sys.stderr, self.name,Sbound, Vbound, max(fishS),max(fishV)
            # calculate constraints
            pos = numpy.where(numpy.logical_and(self.loops>Sbound,self.stems<=Vbound))[0]
            constraints[pos] = ['x' if self.loops[i]>=fishV[self.stems[i]] else '.' for i in pos]
            pos = numpy.where(numpy.logical_and(self.stems>Vbound,self.loops<=Sbound))[0]
            constraints[pos] = ['|' if self.stems[i]>=fishS[self.loops[i]] else '.' for i in pos]
            return FastC(self.name, self.seq, constraints.tostring())
        if method == 'l' or method == 'logS/V': # log(S+1)/(V+1)
            pars = numpy.log2((self.loops+1)/(self.stems+1))
            constraints[numpy.where(pars >  sthreshold)] = 'x'
            constraints[numpy.where(pars <  -vthreshold)] = '|'
            return FastC(self.name, self.seq, constraints.tostring())
        return FastC(self.name, self.seq, constraints.tostring())

class FastC(FastA):
    '''
    FastC format for RNA structure prediction.
    Example file:
        >YIL140W
        AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
        .|xx.|.x..||.|...xx.x.xx...||.x|..||...|..........
    Note:
        '.' means free to fold
        'x' means loops
        '|' means stems

    '''
    def __init__(self, name, seq, constraints=None):
        super(FastC,self).__init__(name,seq)
        self.seq = self.seq.replace('T','U')
        self.constraints = constraints
        if constraints is None:
            self.constraints = '.'*len(self.seq)
    def __str__(self):
        ''' FastC string. '''
        return ">{0}\n{1}\n{2}".format(self.name,self.seq,self.constraints)
    def rc(self):
        ''' Reverse complementary. '''
        super(FastC,self).rc()
        self.constraints = self.constraints[::-1]

class FastS(FastA):
    ''' 
    RNA structure.
    Example: Multiple structures are supported.
        >YIL140W
        AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
        .......((((((((............)))).))))..............  (-6.0)
        .......(((((((..............))).))))..............  (-4.0)
    '''
    def __init__(self, name, seq, structures, scores=None):
        '''
        Initiation.
        Parameters:
            name: string
                Name of sequence.
            seq: string
                RNA sequences.
            structures: string or list of strings
                RNA structure in Dot-Bracket format.
            scores: None, float or list of float
                A score of each structure. Could be the MFE score, the percentage or any other scores. If it is None, each structure will be set to equal percentages.
        '''
        super(FastS,self).__init__(name,seq)
        self.seq = self.seq.replace('T','U')
        if isinstance(structures,basestring):
            self.structures = [structures]
        else:
            self.structures = [st for st in structures]
        if scores is None:
            self.scores = [1./len(self.structures)] * len(self.structures)
        elif isinstance(scores,int) or isinstance(scores,float):
            self.scores = [float(scores)]
        else:
            self.scores = [float(sc) for sc in scores]
    def __str__(self):
        ''' FastS string. '''
        lstr = ">{0}\n{1}".format(self.name,self.seq)
        for st, sc in izip(self.structures, self.scores):
            lstr += "\n{0}\t({1})".format(st, round(sc,3))
        return lstr
    def rc(self):
        ''' Reverse complementary. '''
        super(FastS,self).rc()
        for i in range(len(self.structures)):
            self.structures[i] = self.structures[i][::-1]

class IO(object):
    ''' Process raw data to FastC format. '''
    def FastAReader(infile):
        '''Read sequence files.'''
        # Read lines
        with open(infile) as fh:
            line = fh.next()
            if line[0] != ">":
                raise ValueError("Records in FastA files should start with '>' character")
            line = line.lstrip('>').split()
            name = line[0]
            seq = ''
            while True:
                try:
                    line = fh.next()
                except:
                    if seq != '':
                        yield FastA(name, seq)
                    raise StopIteration
                if line[0] != ">":
                    seq += line.rstrip()
                else:
                    yield FastA(name, seq)
                    line = line.lstrip('>').split()
                    name = line[0]
                    seq = ''
    FastAReader=staticmethod(FastAReader)
    def FastDReader(infile):
        ''' Read FastD file (sequences, S1 depth and V2 depth) and calculate the constraints with threshold. '''
        with open(infile) as fh:
            for line in fh:
                name = line.lstrip(">").rstrip()
                seq  = fh.next().rstrip()
                line = fh.next().rstrip()
                if line[0].startswith('+'):
                    S = numpy.zeros(len(seq))
                else:
                    S = numpy.array([float(i) for i in line.split(";")])
                line = fh.next().rstrip()
                if line[0].startswith('+'):
                    V = numpy.zeros(len(seq))
                else:
                    V = numpy.array([float(i) for i in line.split(";")])
                yield FastD(name, seq, S, V)
    FastDReader=staticmethod(FastDReader)
    def FastCReader(infile):
        ''' Read FastC file. '''
        with open(infile) as fh:
            for line in fh:
                name = line.lstrip(">").rstrip()
                seq  = fh.next().rstrip()
                constraints = fh.next().rstrip()
                yield FastC(name, seq, constraints)
    FastCReader=staticmethod(FastCReader)
    def FastSReader(infile):
        ''' Read FastS file. '''
        with open(infile) as fh:
            name = ''
            seq = ''
            structures = []
            scores = []
            for line in fh:
                if line.startswith('>'):
                    if name != '':
                        yield FastS(name, seq, structures, scores)
                    name = line.lstrip('>').rstrip()
                    seq  = fh.next().rstrip()
                    structures = []
                    scores = []
                else:
                    lstr = line.split()
                    structures.append(lstr[0])
                    try:
                        scores.append(float(lstr[1].rstrip(")").lstrip("(")))
                    except:
                        scores.append(1.0)
            if name != '':
                yield FastS(name, seq, structures, scores)
        raise StopIteration
    FastSReader=staticmethod(FastSReader)
    def Reader(infile, ftype = 'guess'):
        ''' 
        Read multifold format files. 
        Parameters:
            infile: string
                Input file name.Supported file types: FastA, FastD, FastC, FastS and EFastS.
            ftype: string
                By default, the file type can be guessed from the extensions.
                Supported file types: FastA, FastD, FastC, FastS and EFastS. ftype is not case sensitive.
        '''
        ftype = ftype.lower()
        rdict = {'fs':IO.FastSReader,'fasts':IO.FastSReader,'fa':IO.FastAReader,'fasta':IO.FastAReader,'fd':IO.FastDReader,'fastd':IO.FastDReader,'fc':IO.FastCReader,'fastc':IO.FastCReader,'efs':IO.FastSReader,'efasts':IO.FastSReader}
        if ftype == "guess":
            ftype = os.path.splitext(infile)[1][1:]
        if not rdict.has_key(ftype):
            raise IOError("ERROR: File type cannot be guessed from extension.")
        reader = rdict[ftype]
        for item in reader(infile):
            yield item
    Reader=staticmethod(Reader)

class Predictor(object):
    '''
    RNA structure predictor.
    Usage:
        structure, score = Predictor.RNAfold(tFastC, Temperature, threshold)
        print score, structure
        
    '''
    def FoldMerge(tfastc, T = 37, threshold = 0., predictors=['RNAfold','Fold','pknots','UNAFold','sfold'],**kwargs):
        ''' Merge the folded structures from multiple software. '''
        predictor_dict = {'RNAfold': Predictor.RNAfold, 'Fold': Predictor.Fold, 'pknots':  Predictor.pknots, 'mfold':  Predictor.mfold, 'UNAFold':  Predictor.UNAFold, 'sfold':Predictor.sfoldext,'pknotsRG':Predictor.pknotsRG, 'ipknot':Predictor.ipknot}
        structures = []
        scores = []
        for predictor in predictors:
            if predictor_dict.has_key(predictor):
                tfasts = predictor_dict[predictor](tfastc,T,**kwargs)
                structures.extend(tfasts.structures)
                scores.extend(tfasts.scores)
        # unique structure
        struct_dict= {}
        for st, sc in izip(structures,scores):
            if struct_dict.has_key(st):
                struct_dict[st] = min(sc, struct_dict[st]) # energy the lower the better
            else:
                struct_dict[st] = sc
        return FastS(tfastc.name, tfastc.seq, struct_dict.keys(), struct_dict.values())
    FoldMerge=staticmethod(FoldMerge)
    def RNAfold(tfastc, T = 37, threshold = 0,**kwargs):
        ''' Call RNAfold to predict the structure. '''
        # Run RNAfold with contraints
        p = Popen(['RNAfold','-C', '-T', str(T)], stdin = PIPE, stdout = PIPE, stderr = PIPE)
        p.stdin.write(str(tfastc))
        pstdout, pstderr = p.communicate()
        
        # Check the threshold here
        # to do 

        # Extract structure and score
        pre = re.compile(r'(\S+)\s+\((\s*-*\d+\.\d+)\)')
        s = pstdout.rstrip().split('\n')[-1]
        structure , score = pre.search(s).groups()
        if os.path.isfile(tfastc.name+"_ss.ps"):
            os.remove(tfastc.name+"_ss.ps")
        return FastS(tfastc.name, tfastc.seq, [structure], [float(score)])
    RNAfold=staticmethod(RNAfold)
    def Fold(tfastc, T = 37, threshold = 0,**kwargs):
        '''
        Call Fold (RNAStructure package) to predict structures.
        Note: Fold won't give structures if it cannot satisfy all the constraints. [Fix this in the future].
        '''
        # check DATAPATH environment variable
        if not os.environ.has_key('DATAPATH'):
            raise KeyError("Fold (RNAStructure) Error: Please set environment variable $DATAPATH to the location of the data_tables.")
        T = Utils.TempConverter(T,'C','K') # from C to K
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa')
        (fdct, ctname) = tempfile.mkstemp(suffix='.ct')
        (fdcn, cnname) = tempfile.mkstemp(suffix='.CON')
        with os.fdopen(fdfa, 'w') as fh:
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        with os.fdopen(fdcn, 'w') as fh:
            fh.write(Utils.FastCToFold(tfastc.constraints))
        os.close(fdct)
        structures = []
        scores = []
        try:
            # Run Fold
            # Fold 
            p = Popen(['Fold', faname, ctname, '-T', str(T), '-C', cnname], stdin = None, stdout = PIPE, stderr = PIPE)
            pstdout, pstderr = p.communicate()
            if pstderr:
                raise ValueError("Fold (RNA Structure package) run error: {0}".format(pstderr))
            # Parse the *.ct file
            tfs = Utils.ct2dot(ctname)
        finally:
            os.remove(faname)
            os.remove(ctname)
            os.remove(cnname)
        if len(structures) == 0:
            print >> sys.stderr, "No structures predicted by Fold."
        return tfs # FastS(tfastc.name, tfastc.seq, structures, scores)
    Fold=staticmethod(Fold)
    def pknots(tfastc, T = 37, threshold = 0,**kwargs):
        '''
        Fold by pknots. 
        pknots doesn't have temperature parameter.
        pknots has -k parameter for pseudoknots folding.
        '''
        if T != 37:
            print >> sys.stderr, "WARNING: pknots doesn't accept temperature != 37 !"
            return FastS(tfastc.name, tfastc.seq, [], [])
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa')
        (fdct, ctname) = tempfile.mkstemp(suffix='.ct')
        with os.fdopen(fdfa, 'w') as fh:
            fh.write(">{0}\n{1}".format(tfastc.name,tfastc.seq))
        os.close(fdct)
        structures = []
        scores = []
        try:
            # Run pknots
            p = Popen(['pknots', '-k', faname, ctname], stdin = None, stdout = PIPE, stderr = PIPE)
            pstdout, pstderr = p.communicate()
            if pstderr:
                raise ValueError("ERROR: pknots run error: {0}".format(pstderr))
            # Parse the *.ct file
            #structures, scores = Utils.ct2dot(ctname)
            structures, scores = Utils.pknots2dot(ctname)
        finally:
            os.remove(faname)
            os.remove(ctname)
        if len(structures) == 0:
            print >> sys.stderr, "No structures predicted by pknots."
        return FastS(tfastc.name, tfastc.seq, structures, scores)
    pknots=staticmethod(pknots)
    def pknotsRG(tfastc, T = 37, threshold = 0,**kwargs):
        ''' Fold by pknotsRG. '''
        if T != 37:
            print >> sys.stderr, "WARNING: pknotsRG doesn't accept temperature != 37 !"
            return FastS(tfastc.name, tfastc.seq, [], [])
        structures = []
        scores = []
        try:
            p = Popen(['pknotsRG', '-s'], stdin = PIPE, stdout = PIPE, stderr = PIPE)
            p.stdin.write(">{0}\n{1}".format(tfastc.name, tfastc.seq))
            pstdout, pstderr = p.communicate()
        except OSError:
            raise OSError('pknotsRG is not installed.')
        if pstderr:
            raise IOError('ERROR: pknotsRG run error: '+pstderr)
        # parse the output
        for line in pstdout.split('\n')[6:-1]: # first 5 lines are descriptions
            structure, score = line.split()
            score = float(score[1:-2]) ## (-5.80)
            if score < 0 and structure not in structures:
                structures.append(structure)
                scores.append(score) 
        return FastS(tfastc.name, tfastc.seq, structures, scores)
    pknotsRG=staticmethod(pknotsRG)
    def mfold(tfastc, T = 37, threshold = 0,**kwargs):
        '''
        Fold by mfold.
        mfold SEQ=test.fa AUX=test.con T=37 RUN_TYPE=html
        '''
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa')
        (fdcn, cnname) = tempfile.mkstemp(suffix='.con')
        with os.fdopen(fdfa, 'w') as fh: # with statement will close the fh and fd in the end
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        with os.fdopen(fdcn, 'w') as fh:
            fh.write(Utils.FastCToMfold(tfastc.constraints))
        structures = []
        scores = []
        try:
            # mfold 
            p = Popen(['mfold', 'SEQ='+faname, 'AUX='+cnname, 'T={0}'.format(T), 'RUN_TYPE=html'], stdin = None, stdout = PIPE, stderr = PIPE) 
            pstdout, pstderr = p.communicate()
            if pstderr:
                raise ValueError("ERROR: mfold run error: {0}".format(pstderr))
        finally:
            os.remove(faname)
            os.remove(cnname)
            # Parse the *ct file
            idx = 1
            bname = os.path.basename(faname)
            fn = "{0}_{1}.ct".format(bname,idx)
            while os.path.isfile(fn):
                tfs = Utils.ct2dot(fn)
                structures.extend(tfs.structures)
                scores.extend(tfs.scores)
                idx += 1
                fn = "{0}_{1}.ct".format(bname,idx)
            for fn in os.listdir("."):
                if fn.startswith(bname):
                    os.remove(fn)
        if len(structures) == 0:
            print >> sys.stderr, "No structures predicted by mfold."
        if os.path.isfile('date.test'):
            os.remove('date.test')
        return FastS(tfastc.name, tfastc.seq, structures, scores) 
    mfold=staticmethod(mfold)
    def UNAFold(tfastc, T = 37, threshold = 0,**kwargs):
        ''' Fold by UNAFold. '''
        # Create temp file for fa and constraints
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa',dir=".")
        (fdcn, cnname) = tempfile.mkstemp(suffix='.aux',dir=".")
        faname = os.path.basename(faname)
        cnname = os.path.basename(cnname)
        with os.fdopen(fdfa, 'w') as fh:
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        with open(cnname, 'w') as fh:
            fh.write(Utils.FastCToMfold(tfastc.constraints))
        structures = []
        scores = []
        try:
            # mfold 
            p = Popen(['UNAFold.pl', '-t', str(T), '-c', cnname, '--run-type=html', faname], stdin = None, stdout = PIPE, stderr = PIPE) 
            pstdout, pstderr = p.communicate()
            if pstderr:
                print >> sys.stderr, "ERROR: UNAFold run error: {0}".format(pstderr)
        finally:
            os.remove(faname)
            os.remove(cnname)
            # Parse the *ct file
            bname = os.path.basename(faname)
            for fn in os.listdir("."):
                if fn.startswith(bname):
                    if fn.endswith(".ct"):
                        tfs = Utils.ct2dot(fn)
                        for st,sc in izip(tfs.structures,tfs.scores):
                            if st not in structures:
                                structures.append(st)
                                scores.append(sc)
                    os.remove(fn)
        if len(structures) == 0:
            print >> sys.stderr, "No structures predicted by UNAFold."
        return FastS(tfastc.name, tfastc.seq, structures, scores) 
    UNAFold=staticmethod(UNAFold)
    def ipknot(tfastc, T = 37, threshold = 0, **kwargs):
        ''' ipknot fold. '''
        if T != 37:
            print >> sys.stderr, "WARNING: pknotsRG doesn't accept temperature != 37 !"
            return FastS(tfastc.name, tfastc.seq, [], [])
        (fdfa, faname) = tempfile(suffix='.fa')
        with os.fdopen(fdfa) as fh:
            fh.write(">{0}\n{1}".format(tfastc.name, tfastc.seq))
        try:
            p = Popen(['ipknot', '-r', '10', faname], stdout = PIPE, stderr = PIPE)
            pstdout, pstderr = p.communicate()
        except OSError:
            raise OSError('ERROR: ipknot is not installed.')
        finally:
            os.remove(faname)
        # Parse the result
        if pstderr:
            raise IOError('ERROR: ipknot run error: '+pstderr)
        structure =  pstdout.split('\n')[-2].split()
        if structure.strip(".") != "":
            return FastS(tfastc.name, tfastc.seq, [structure], [0.0])
        return FastS(tfastc.name, tfastc.seq, [], [])
    ipknot=staticmethod(ipknot)
    def sfold(tfastc, T = 37, threshold = 0,**kwargs):
        ''' Fold by sfold. '''
        # Test if sfold exists
        try:
            p = Popen(['sfold'],stdin=None,stdout=PIPE,stderr=PIPE)
        except:
            raise ValueError("ERROR: sfold cannot be found. Download it from http://sfold.wadsworth.org/License_info.html\
                              and make it public accessible by\n   >export PATH=$PATH:SFOLD_HOME/bin\n")
        # Create temp file for fa and constraints
        prefix = tempfile.mkdtemp(dir=".")
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa',dir='.')
        (fdcn, cnname) = tempfile.mkstemp(suffix='.aux',dir='.')
        with os.fdopen(fdfa, 'w') as fh:
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        # constraints
        cns = False
        if '|' in tfastc.constraints or 'x' in tfastc.constraints:
            cns = True
            with os.fdopen(fdcn, 'w') as fh:
                fh.write(Utils.FastCToMfold(tfastc.constraints))
        structures = []
        scores = []
        try:
            # sfold
            if cns:
                p = Popen(['sfold', '-f', cnname, '-o', prefix, faname], stdin = None, stdout = PIPE, stderr = PIPE)
            else:
                p = Popen(['sfold', '-o', prefix, faname], stdin = None, stdout = PIPE, stderr = PIPE)
            pstdout, pstderr = p.communicate()
            if pstderr:
                print >> sys.stderr, "ERROR: sfold run error: {0}".format(pstderr)
            if os.path.isdir(prefix+"/clusters"): # clusters exists
                for f in os.listdir(prefix+"/clusters"):
                    if f.endswith(".ct"):
                        tfs = Utils.ct2dot(prefix+"/clusters/"+f)
                        for st,sc in izip(tfs.structures,tfs.scores):
                            structures.append(st)
                            scores.append(sc)
        finally:
            os.remove(faname)
            os.remove(cnname)
            os.system('rm -rf '+prefix)
        return FastS(tfastc.name,tfastc.seq, structures, scores)
    sfold=staticmethod(sfold)
    def sfoldext(tfastc,T=37,threshold = 0, **kwargs):
        '''
        Extended verision of sfold, in which the number of structures can be specified. 
            Parameters:
                tfastc: FastC object.
                    Fastc object to be folded. Constraints are set to '.' if not applicable.
                T: float
                    Temperature for fold. Not applicable here.
                threshold: float
                    threshold of energy of structures.
                kwargs: dict
                    dictionary of additional parameters. 
                    "N": integer, number of sample size. Default is 1000.
                    "workdir": string, directory for temporary files during fold. Default is "workdir"
                    "keep": bool, keep the workdir or not. Default is False to save space.
            Returns:
                fs: FastS object
                    Centroid of structures.
       ''' 
        # Test if sfold exists
        p = Popen(['which','sfold'],stdin = None, stdout = PIPE, stderr = PIPE)
        pstdout, pstderr = p.communicate()
        if pstdout == '':
            raise IOError("ERROR: sfold cannot be found. Download it from http://sfold.wadsworth.org/License_info.html and make it public accessible by\n   >export PATH=$PATH:SFOLD_HOME/bin")
        # set parameters        
        envfile = pstdout[:-10]+"/sfoldenv"
        if not os.path.isfile(envfile):
            #os.system('cd /net/uu/nm/bi/yxw120430/software/RNAStructure/sfold-2.2 && ./configure')
            raise IOError("ERROR: sfoldenv cannot be found. Please go to sfold home directory and run >./configure")
        with open(envfile) as fh:
            for line in fh:
                if "=" in line:
                    var, val = line.rstrip().split('=')
                    os.environ[var] = os.path.expandvars(val)

        # Create file for fa and constraints
        workdir = os.path.expanduser(kwargs.get('workdir','./workdir'))
        keep = kwargs.get('keep',False)

        prefix = workdir+"/"+tfastc.name
        faname = prefix+".fa"
        cnname = prefix+".aux"
        with open(faname, 'w') as fh:
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        # constraints
        cns = False
        if '|' in tfastc.constraints or 'x' in tfastc.constraints:
            cns = True
            with open(cnname, 'w') as fh:
                fh.write(Utils.FastCToMfold(tfastc.constraints))
        structures = []
        scores = []
        N = int(kwargs.get("N",1000))
        try:
            # sfold
            if debug: start = time.time()
            if cns:
                cmd = [os.environ['SFOLDEXE'], '-n', str(N), '-p', os.environ['SFOLDPAR'], '-f', cnname, '-o', prefix, faname]
                p = Popen(cmd, stdin = None, stdout = PIPE, stderr = PIPE)
            else:
                cmd = [os.environ['SFOLDEXE'], '-n', str(N), '-p', os.environ['SFOLDPAR'], '-o', prefix, faname]
                p = Popen(cmd, stdin = None, stdout = PIPE, stderr = PIPE)
            # print ' '.join(cmd)
            pstdout, pstderr = p.communicate()
            if not os.path.isfile(prefix+"/sample_1000.out"): # samples exist
                raise IOError("ERROR: {0} sampling failed for {1}.".format(os.environ['SFOLDEXE'],tfastc.name))
            if debug: stop1 = time.time()
            # SCLASS
            os.system('mkdir '+prefix+'/clusters')
            cmd = [os.environ['SCLASSEXE'], '-b', prefix+'/bp.out', '-c', prefix+'/clusters', '-e', prefix, '-f', prefix+'/fe.out', '-p', prefix+'/bprob.out', '-s', faname]
            #print ' '.join(cmd)
            p = Popen(cmd,stdin = None, stdout = PIPE, stderr = PIPE)
            pstdout, pstderr = p.communicate()
            #if not pstderr.startswith('Read'):
                #sys.stderr.write(pstderr)
            if not os.path.isfile(prefix+"/clusters/c01.ccentroid.ct"):
                raise IOError("ERROR: {0} clustering failed for {1}.".format(os.environ['SCLASSEXE'],tfastc.name))
            if debug: 
                stop2 = time.time()
                print >> sys.stderr, tfastc.name,stop1-start,stop2-stop1
            # read centroids
            if os.path.isdir(prefix+"/clusters"): # clusters exists
                for f in os.listdir(prefix+"/clusters"):
                    if f.endswith(".ct"):
                        tfs = Utils.ct2dot(prefix+"/clusters/"+f)
                        for st,sc in izip(tfs.structures,tfs.scores):
                            structures.append(st)
                            scores.append(sc)
        finally:
            if not keep:
                os.system('rm -rf {0}*'.format(prefix))
        if len(structures) == 0:
            raise ValueError("ERROR: no structure was predicted for {0}.".format(tfastc.name))
        fs = FastS(tfastc.name,tfastc.seq, structures, scores)
        return fs
    sfoldext=staticmethod(sfoldext)
    def snoGPS(tfc, T = 37, threshold = 0,**kwargs):
        '''
        Predict H/ACA box snoRNAs.
        '''
        pass
        return tfs
    snoGPS=staticmethod(snoGPS)
    def sfold_star(fc_fh):
        '''
        star wrapper of sfold.
        '''
        fc,fh = fc_fh
        print >>fh, sfold(fc)
        return

class Algorithm(object):
    ''' Algorithm for RNA structure prediction. '''
    def simulation(fs,numR,pS=0.5,noise_ratio=0.):
        '''
        Simulate reads for structures with a given percentages.
            fs: FastS object
                Contains multiple structures and percentage of each structures
            numR: int
                number of reads mapped to the RNA
            pS: float,range in [0,1]
                Percentage of reads for S. pV = 1 - pS
            noise_ratio: float
                Percentage of noise reads
        Returns:
            fp: FastD object
                Randomly attribute reads based on pS, pV and noise_raito.
        Note:
            numR is attributed to loops and stems by nS:nV.
            pS and pV are attributed to structures by scores of structures.
        '''        
        # initiation
        L = len(fs)
        M = len(fs.structures)
        S = numpy.zeros(L)
        V = numpy.zeros(L)
        # ranges of each structure
        scores = numpy.array(fs.scores,dtype=float)
        # Number of noise reads
        numN = int(round(numR*noise_ratio))
        numR -= numN
        if debug: print numR,numN
        
        # Number of RNA molecular for S and V
        pV = 1-pS
       
        # I matrix
        IV = numpy.zeros((L,M))+1
        IS = numpy.zeros((L,M))
        for j in range(M):
            s = numpy.array(list(fs.structures[j]))
            IS[:,j] += s == '.'
            IV[:,j] -= s == '-'
        IV -= IS
        IN = IS+IV
        # Effective length
        LS = numpy.sum(IS,0) # loop effective length
        LV = numpy.sum(IV,0) # stem effectvie length
        LN = numpy.sum(IN,0) # noise effective length
        # theta
        theta_S = scores*LS
        theta_S /= sum(theta_S)
        print theta_S
        theta_V = scores*LV
        theta_V /= sum(theta_V)
        theta_N = LN/sum(LN)
        print theta_V
                    
        # generate random reads for loops and stems
        numpy.random.seed(int(time.time()*100)%100)
        nS = numpy.round(numR*pS*theta_S)
        nV = numpy.round(numR*pV*theta_V)
        nNS = numpy.round(numN*pS*theta_N)
        nNV = numpy.round(numN*pV*theta_N)
        for j in xrange(M):
            idx = numpy.nonzero(IS[:,j])[0] # loop reads
            for i in numpy.random.randint(LS[j],size=int(nS[j])):
                S[idx[i]] += 1
            idx = numpy.nonzero(IV[:,j])[0] # stem reads
            for i in numpy.random.randint(LV[j],size=int(nV[j])):
                V[idx[i]] += 1
            idx = numpy.nonzero(IN[:,j])[0] 
            for i in numpy.random.randint(LN[j],size=int(nNS[j])): # noise loop reads
                S[idx[i]] += 1
            for i in numpy.random.randint(LN[j],size=int(nNV[j])): # noise stem reads
                V[idx[i]] += 1
        return FastD(fs.name,fs.seq,S,V)
    simulation=staticmethod(simulation)
    def EM(fd,fs,centroids=None,exprs=None,threshold = 1e-4, maxiter=100):
        ''' 
        EM algorithm to calculate percentages and expression levels for all centroid structures.. 
        Parameters:
            fd: FastD object
                RNA footprinting data
            fs: FastS object
                Predicted structures
            exprs: tuple or None
                None means no exprs
                exprss in format: ((E1,n1),(E2,n2),...,(El,nl)).Ei is expression level and ni is number of structures for the ith exprs.
            threshold: float
                Threshold for convergence
            maxiter: int
                Maximum iteration number
        Returns:
            Pi: list of float
                Percentages of each structure.
            uS: list of float
                Estimated expression level from loop data for each structure
            uV: list of float
                Estimated expression level from stem data for each structure
        '''
        def normalize(x,p=1):
            s = sum(x*p)
            if s == 0:
                return x
            else:
                return x*p/s
        # check selected centroids
        if centroids is None:
            centroids = range(len(fs.structures))
        # check exprs
        if not exprs is None:
            Es,Ns = zip(*exprs)
            Es = normalize(Es)
        M = len(centroids) # Number of clusters
        if M == 0:
            return [],[],[]
        L = len(fd.seq)  # Length of RNA
        
        S = numpy.array(fd.loops,dtype=float)
        V = numpy.array(fd.stems,dtype=float)
       
        # Iij, effective positions
        IV = numpy.zeros((L,M))+1
        IS = numpy.zeros((L,M))
        for j in range(M):
            s = numpy.array(list(fs.structures[centroids[j]]))
            IS[:,j] += s == "."
            IV[:,j] -= s == "-"
        IV -= IS
        # effective length
        LS = numpy.sum(IS,0)
        LV = numpy.sum(IV,0)
        
        if M == 1:
            return [1.0],[numpy.mean(S[IS[:,0]==1])],[numpy.mean(V[IV[:,0]==1])]
        # Initial percentage
        Pi = numpy.zeros(M)+1.0/M
        
        # iteration
        cnt = 0
        while (True):
            # E step: E[P(X,Z|theta)]
            # normalize the possibilities by row
            ES = numpy.apply_along_axis(lambda x:normalize(x,Pi),1,IS)
            EV = numpy.apply_along_axis(lambda x:normalize(x,Pi),1,IV)
            # expected reads for each structure
            muS = numpy.sum(numpy.apply_along_axis(lambda x: x*S,0,ES),0)
            muV = numpy.sum(numpy.apply_along_axis(lambda x: x*V,0,EV),0)
            # test Q(j)
            mu = muS/LS + muV/LV
            # M step: argmax(E[P(X,Z|theta)])
            if exprs is None:
                mu /= sum(mu)
            else:
                cnt = 0
                for e,n in zip(Es,Ns):
                    mu[cnt:cnt+n] *= e/sum(mu[cnt:cnt+n])
                    cnt +=n            
            
            # check if converge
            if cnt > maxiter or numpy.sum(numpy.abs(mu-Pi)) < threshold:
                break
            Pi = mu
            cnt += 1
        if debug: print >> sys.stderr, "EM iterates {0} times.".format(cnt)
        return Pi, muS/LS,muV/LV
    EM=staticmethod(EM)
    def fitness(fd,fs,threshold=0):
        '''
        Calculate fitness, the percentage of compatible reads, given all reads and structures. 
        Parameters: 
            fd: FastD object
                FastD object. Flanking bases and outliers are trimmed.
            fs: FastS object
                FastS with known or predicted structures.
            treshold: float
                Only structures with percentages higher or equal to percentages are considered. 
                Set threshold to 0 if all structures are wanted. This is also used when negative MFE scores are used.
        Output:
            fitS: float
                fitness score of loops (S reads)
            fitV: float
                fitness score of stems (V reads)
            fitness: float
                overall fitness score
        '''
        # indice of wanted structures
        idx = numpy.arange(len(fs.structures))
        if threshold > 0:
            idx = idx[numpy.array(fs.scores)>=threshold]
        # Count loops and stems
        L = len(fs)
        IS = numpy.zeros(L) != 0
        IV = numpy.zeros(L) != 0
        for i in idx:
            s = numpy.array(list(fs.structures[i]))
            IS |= s == "."
            IV |= (s !='.') & (s !='-')
        # fitness scores
        loops = numpy.array(fd.loops,dtype=float)
        stems = numpy.array(fd.stems,dtype=float)
        S = sum(loops)
        V = sum(stems)
        fS = sum(loops[IS])
        fV = sum(stems[IV])
        return S>0 and fS/S or 0, V>0 and fV/V or 0, S+V>0 and (fS+fV)/(S+V) or 0
    fitness=staticmethod(fitness)
    def bpReader(prefix):
        '''
        Read bp.out and centroids from Boltzmann sampling. The centroids are at the end of the structures.
        Parameters:
            prefix: string
                Folder name of the predicted structures, which contains the 'clusters' folder.
        Returns:
            fs: FastS
                FastS object with structures in bp.out
        '''    
        fafh = ngslib.IO.BioReader(prefix+".fa",'fasta')
        fa = fafh.next()
        fafh.close()
        structures = []
        l = len(fa)
        pairs = numpy.zeros(l,dtype=int)
        # read structures from bp.out
        for line in open(prefix+"/bp.out"):
            if not line.startswith(' Structure'):
                s,e = [int(i) for i in line.split()]
                pairs[s-1] = e
                pairs[e-1] = s
            else:
                if max(pairs) != 0:
                    structures.append(Utils._pairsToDot(pairs))
                pairs = numpy.zeros(l,dtype=int)        
        if max(pairs) != 0:
            structures.append(Utils._pairsToDot(pairs))
        # read structures from centroids
        for fn in os.listdir(prefix+"/clusters"):
            if fn.endswith('.ct'):
                tfs = Utils.ct2dot(prefix+"/clusters/"+fn)
                structures.append(tfs.structures[0])
        return FastS(fa.id,str(fa.seq),structures)
    bpReader=staticmethod(bpReader)    
    def boltzman_sampling(tfastc,T,N=1000):
        # Test if sfold exists
        p = Popen(['which','sfold'],stdin = None, stdout = PIPE, stderr = PIPE)
        pstdout, pstderr = p.communicate()
        if pstdout == '':
            raise IOError("ERROR: sfold cannot be found. Download it from http://sfold.wadsworth.org/License_info.html and make it public accessible by\n   >export PATH=$PATH:SFOLD_HOME/bin")
        # set parameters        
        envfile = pstdout[:-10]+"/sfoldenv"
        if not os.path.isfile(envfile):
            #os.system('cd /net/uu/nm/bi/yxw120430/software/RNAStructure/sfold-2.2 && ./configure')
            raise IOError("ERROR: sfoldenv cannot be found. Please go to sfold home directory and run >./configure")
        with open(envfile) as fh:
            for line in fh:
                if "=" in line:
                    var, val = line.rstrip().split('=')
                    os.environ[var] = os.path.expandvars(val)

        # Create temp file for fa and constraints
        prefix = tempfile.mkdtemp(dir=".")
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa',dir='.')
        (fdcn, cnname) = tempfile.mkstemp(suffix='.aux',dir='.')
        with os.fdopen(fdfa, 'w') as fh:
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        # constraints
        cns = False
        if '|' in tfastc.constraints or 'x' in tfastc.constraints:
            cns = True
            with os.fdopen(fdcn, 'w') as fh:
                fh.write(Utils.FastCToMfold(tfastc.constraints))
        structures = []
        scores = []
        try:
            # sfold
            if cns:
                cmd = [os.environ['SFOLDEXE'], '-n', str(N), '-p', os.environ['SFOLDPAR'], '-f', cnname, '-o', prefix, faname]
                p = Popen(cmd, stdin = None, stdout = PIPE, stderr = PIPE)
            else:
                cmd = [os.environ['SFOLDEXE'], '-n', str(N), '-p', os.environ['SFOLDPAR'], '-o', prefix, faname]
                p = Popen(cmd, stdin = None, stdout = PIPE, stderr = PIPE)
            print ' '.join(cmd)
            pstdout, pstderr = p.communicate()
            if not os.path.isfile(prefix+"/sample_1000.out"): # samples exist
                raise IOError("ERROR: {0} sampling failed.".format(os.environ['SFOLDEXE']))
                if os.path.isfile(prefix+"/sample_1000.out"): # samples exist
                    with open(prefix+"/sample_1000.out") as fh:
                        B = [0]*len(tfastc)
                        for line in fh:
                            fields = line.split()
                            if len(fields) == 0:
                                pass
                            elif fields[0] == 'Structure':
                                if fields[1] != '1':
                                    structures.append(Utils._pairsToDot(B))
                                    B = [0]*len(tfastc)
                            elif len(fields) == 3:
                                a = int(fields[0])
                                b = int(fields[1])
                                for i in range(int(fields[2])):
                                    B[a+i-1]=b-i
                        structures.append(Utils._pairsToDot(B))
        finally:
            os.remove(faname)
            os.remove(cnname)
            os.system('rm -rf '+prefix)
        scores = numpy.ones(N)
        return FastS(tfastc.name,tfastc.seq, structures, scores)
    boltzman_sampling=staticmethod(boltzman_sampling)
    def fdr_bh(pvalues):
        n = len(pvalues)
        if n == 0:
            return []
        # sort pvlaues
        args, pv = zip(*sorted(enumerate(pvalues),key= lambda x:x[1]))
        bh_values = [0] * n
        pre_v = 0
        for i, pv in enumerate(pv):
            cur_v = min(1, pv * n /(i+1)) # correction
            pre_v = max(cur_v, pre_v)
            bh_values[args[i]] = pre_v
        return bh_values
    fdr_bh=staticmethod(fdr_bh)
    def rank(x):
        ''' rank of a list. '''
        array=numpy.array(x)
        temp = array.argsort()
        ranks = numpy.arange(len(array))[temp.argsort()]
        return ranks
    rank=staticmethod(rank)
    def trimed(x, lowbound = 0.05, upbound = 0.95):
        ''' Trim the top ones and the bottom ones, and set them to bound values. '''
        r = Algorithm.rank(x) # type(r) is numpy.ndarray
        l = len(r) - 1
        ub = x[numpy.where(r == round(upbound * l))[0]]
        lb = x[numpy.where(r == round(lowbound *l))[0]]
        if type(x) == numpy.ndarray:
            nx = numpy.array(x)
            nx[nx>ub] = ub
            nx[nx<lb] = lb
        else:
            nx = [ ub if i > ub else i for i in x]
            nx = [ lb if i < lb else i for i in nx]
        return nx
    trimed=staticmethod(trimed)
    def distanceMatrix(fs,method="bp"):
        '''
        Calculate the distance matrix, given a FastS object with multiple structures.
        Options:
            method: 'top' means calculate the topological distance using wRNA.distance from Vienna RNA package [1].
                    'rbp' means calculate the relaxed base-pair distance [2].
                    'bp'  means calculate the base-pair distance [3].
        References:
            [1] Hofacker, I. L., et al. (1994). "Fast Folding and Comparison of Rna Secondary Structures." Monatshefte Fur Chemie 125(2): 167-188.
            [2] Agius, P., et al. (2010). "Comparing RNA secondary structures using a relaxed base-pair score." RNA 16(5): 865-878.
            [3] Ding, Y., et al. (2005). "RNA secondary structure prediction by centroids in a Boltzmann weighted ensemble." RNA 11(8): 1157-1166.
        '''
        if method == 'top':
            dist = wRNA.distance
        elif method == 'bp':
            dist = Algorithm.bp_distance2 
        elif method == 'rbp':
            dist = Algorithm.rbp_distance
        else:
            raise ValueError('ERROR: method not recognized. Should be in ["rbp","pos","bp","top"].')
        if len(fs.structures)==1:
            return numpy.matrix([0])
        dim = len(fs.structures)
        dismatrix = numpy.zeros((dim,dim))
        for i in xrange(dim):
            for j in xrange(i):
                dismatrix[i,j] = dist(fs.structures[i], fs.structures[j])
                dismatrix[j,i] = dismatrix[i,j]
        return dismatrix
    distanceMatrix=staticmethod(distanceMatrix)
    def rbp_distance(S1,S2, t = 1):
        '''
        Relaxed base-pair distance between two structures [1].
        Options:
            t: relaxed parameter.
        References:
            [1] Agius, P., et al. (2010). "Comparing RNA secondary structures using a relaxed base-pair score." RNA 16(5): 865-878.
        '''
        # basepairs for S1
        bps1 = []
        m = re.search('[\(\)]+',S1)
        pos= 0 
        while m is not None :            
            bps1.append((pos+m.start()+1,pos+m.end()))
            pos+=m.end()
            m = re.search('[\(\)]+',S1[pos:])
        # basepairs for S2
        bps2 = []
        m = re.search('[\(\)]+',S2)
        pos= 0 
        while m is not None :            
            bps2.append((pos+m.start()+1,pos+m.end()))
            pos+=m.end()
            m = re.search('[\(\)]+',S2[pos:])
        # pairwise distances
        ds = []
        for i,j in bps1:
            ds.append( min([max(abs(i-i1),abs(j-j1)) for i1,j1 in bps2]))
        for i,j in bps2:
            ds.append( min([max(abs(i-i1),abs(j-j1)) for i1,j1 in bps1]))
        # sort ds
        ds = numpy.array(sorted(ds,reverse=True))
        # calcluate min m
        minm = len(ds)
        for m in range(len(ds)+1):
            if numpy.all(ds[m:] <= t*m):
                minm = m
                break
        return minm
    rbp_distance=staticmethod(rbp_distance)
    def pos_distance(s1, s2):
        ''' structure distance by position. euclidean distance between two structures.'''
        m1 = numpy.array(list(s1))=='.'
        m2 = numpy.array(list(s2))=='.'
        return numpy.linalg.norm(m1-m2)
    pos_distance=staticmethod(pos_distance)
    def bp_distance2(s1,s2):
        ''' base pair distance defined by Ding, Y., et al. (2005). '''
        return numpy.sum(numpy.abs(Algorithm.dot2mat(s1)-Algorithm.dot2mat(s2)))
    bp_distance2=staticmethod(bp_distance2)
    def dot2mat(s):
        ''' no pseudoknots. '''
        l = len(s)
        m = numpy.zeros((l,l),dtype=bool)
        stack1 = []
        for i,c in enumerate(s):
            if c == '(':
                stack1.append(i)            
            elif c == ')':
                m[stack1.pop(),i] = 1
        return m
    dot2mat=staticmethod(dot2mat)
    def accessibility(fs,start,end):
        ''' 
        Calculate RNA accessibility. 
        Parameters:
            fs: FastS object
                RNA structures with percentages as scores.
            start: int
                Start position, zerobased, included.
            end:   int
                End position, zerobased, not included.
        '''
        acc = 0.
        for st,sc in zip(fs.structures,fs.scores):
            if st[start:end].count('.') == end-start:
                acc += sc
        return acc
    accessibility=staticmethod(accessibility)

class Utils(object):
    ''' Utils for RNA structure prediction. '''
    DegreeAlphabet = ['F','K','C', 'R']
    AnyToC = { 'F': lambda x: (x-32.)/1.8, 'K': lambda x: x-273.16, 'R': lambda x: x*1.25 }
    CToAny = { 'F': lambda x: x*1.8+32., 'K': lambda x: x+273.16, 'R': lambda x: x* 0.8 }
    def TempConverter(temp, _from, _to):
        '''
        Convert temperature among different systems.
            C = Celsius
            F = Fahrenheit
            K = Kelvin (absolute temperature)
            R = Reaumur Equals
        '''
        _from = _from.upper()
        _to = _to.upper()
        if _from not in Utils.DegreeAlphabet or _to not in Utils.DegreeAlphabet:
            raise ValueError("Temperature format not supported. Please choose from ['F','C','K', 'R'].")
        if _from == 'C':
            return Utils.CToAny[_to](temp)
        if _to == 'C':
            return Utils.AnyToC[_from](temp)
        t = Utils.AnyToC[_from](temp)
        return Utils.CToAny[_to](t)
    TempConverter=staticmethod(TempConverter)
    def FastCToMfold(constraints):
        ''' Convert FastC constraints to mfold constraints. '''
        mconstraints = []
        for match in re.finditer(r'(\|+)', constraints):
            mconstraints.append('F {0} 0 {1}'.format(match.start()+1, match.end() - match.start()))
        for match in re.finditer(r'(x+)', constraints):
            mconstraints.append('P {0} 0 {1}'.format(match.start()+1, match.end() - match.start()))
        return '\n'.join(mconstraints)
    FastCToMfold=staticmethod(FastCToMfold)
    def FastCToFold(constraints):
        '''
        Convert FastC constraints to Fold (RNAStructure package) constraints.
        Example of Fold contraints:
            DS:
            2
            3
            -1
            SS:
            5
            7
            -1
        '''
        fconstraints = "DS:"
        for match in re.finditer(r'(\|+)', constraints):
            for i in range(match.start(),match.end()):
                fconstraints += "\n{0}".format(i+1)
        fconstraints += "\n-1\nSS:"
        for match in re.finditer(r'(x+)', constraints):
            for i in range(match.start(),match.end()):
                fconstraints += "\n{0}".format(i+1)
        fconstraints +="\n-1\nMod:\n-1\nPairs:\n-1 -1\nFMN:\n-1\nForbids:\n-1 -1"
        return fconstraints
    FastCToFold=staticmethod(FastCToFold)
    def ct2dot(ctname):
        ''' 
        CT format to DOT format.
        Input can be CT file names or a CT file string.
            50      dG = -6.2       seqname
            1       A       0       2       0       1       0       0
            2       U       1       3       0       2       0       0
            3       G       2       4       0       3       0       0
            4       A       3       5       0       4       0       0
            5       C       4       6       0       5       0       0
            6       A       5       7       0       6       0       0
            7       C       6       8       0       7       0       8
            8       A       7       9       36      8       7       9
            9       G       8       10      35      9       8       10
            10      C       9       11      34      10      9       11
            11      U       10      12      33      11      10      12
            12      U       11      13      31      12      11      13
            13      C       12      14      30      13      12      14
            14      A       13      15      29      14      13      15
        Note:
            (1) dG can be ENERGY, energy or dG.
            (2) the first 6 fields are required.
            (3) for RNA without structure predicted, no dG field provided.
        For pseudoknots like this:
            '(((..[[[...)))...]]]'
            The first seen pairs are labelled with '()', and the following conflict ones are labelled as '[]'
        '''
        structures = []
        scores = []
        try: # file name
            fh = open(ctname)
            cts = fh.readlines()
            fh.close()
        except: # string of CTs
            cts = ctname.split('\n')
        pairs = []
        score = 0.
        seq = ''
        name = ''
        # suppose there might be multiple structures in the ct file.
        idx = 0
        while idx < len(cts):
            line = cts[idx].split()
            if len(line)>=6 and line[0] == line[5]:
                if line[5] == '1': # first line
                    if len(pairs) > 0:
                        s = Utils._pairsToDot(pairs)
                        if '(' in s or '[' in s:
                            structures.append(s)
                            scores.append(score)
                    pairs = [int(line[4])]
                    seq += line[1]
                    # Parse header                
                    name = cts[idx-1].split()[-1]
                    if name.strip('.') == "":
                        name = "NONAME"
                    m = re.findall('-*\d+\.*\d*',cts[idx-1])
                    if len(m) >= 1:
                        ctlen = int(m[0])
                    if len(m) >=2:
                        score = float(m[1])
                else:
                    pairs.append(int(line[4]))
                    seq += line[1]
            idx +=1
        if len(pairs) > 0:
            s = Utils._pairsToDot(pairs)
            if '(' in s or '[' in s:
                structures.append(s)
                scores.append(score)
        return FastS(name,seq,structures,scores) 
    ct2dot=staticmethod(ct2dot)
    def _pairsToDot(pairs,pseudo=False):
        '''
        Generate DOT structure from a list of pairs.
        Parameters:
            pairs: list of stem indice
                Indice in pairs are 1-based. For example, if i is paired with j, then pairs[i-1] = j and pairs[j-1] = i
            pseudo: bool
                Structure has pseudoknot or not. Defalt is False.
        Returns:
            structure: string
                dot format structure.
        '''
        structure = numpy.repeat('.', len(pairs))
        if max(pairs) == 0:
            return structure.tostring()
        if not pseudo:
            for a,b in enumerate(pairs):
                if a < b:
                    structure[a] = '('
                    structure[b-1] = ')'
            return structure.tostring()

        bs = [len(pairs) + 1]
        for a, b in enumerate(pairs):
            a += 1
            if a > b:
                try:
                    bs.remove(a)
                except:
                    pass
            elif b > min(bs):
                structure[a-1] = '['
                structure[b-1] = ']'
            else:
                structure[a-1] = '('
                structure[b-1] = ')'
                bs.append(b)
        return structure.tostring()
    _pairsToDot=staticmethod(_pairsToDot)
    def pknots2dot(pkname):
        ''' Convert pknots output to dot format. '''
        structures = []
        scores = []
        try:
            with open(pkname) as fh:
                pks = fh.readlines()
        except:
            pks = pkname.split('\n')
        idx = 0
        while idx < len(pks):
            # Find SEQ
            while idx < len(pks) and not pks[idx].startswith('SEQ'):
                idx += 1
            if idx >= len(pks):
                break
            # Read pairs
            pairs = []
            idx += 1
            while not pks[idx].startswith('---'):
                pairs.extend([int(i) for i in pks[idx+2].replace('.','0').split()])
                idx += 4
            # Parse pairs
            structures.append(Utils._pairsToDot(pairs))
            # Find energy line
            while not pks[idx].startswith('energy'):
                idx += 1
            scores.append(float(pks[idx].split()[-1]))
        return (structures, scores)
    pknots2dot=staticmethod(pknots2dot)    
    def dot2ct(tfasts): # 
        ''' Convert DOT format to CT format. '''
        ctstring = []
        for st, sc  in izip(tfasts.structures, tfasts.scores):
            # print header
            ctstring.append ("%5d  ENERGY = %-3.2f  %s" % (len(tfasts), sc, tfasts.name))
            stack1=[]
            stack2=[]
            pairs={}
            for i,c in enumerate(st):
                if c == '(':
                    stack1.append(i+1)
                elif c == '[':
                    stack2.append(i+1)
                elif c == ')':
                    pairs[i+1] = stack1.pop()
                    pairs[pairs[i+1]] = i+1
                elif c == ']':
                    pairs[i+1] = stack2.pop()
                    pairs[pairs[i+1]] = i+1
            for i in xrange(1,len(tfasts)+1): # ###24#A######23###25####0###24
                ctstring.append( " %8d %s %8d %8d %8d %8d" % (i, tfasts.seq[i-1], i-1, i+1, pairs.get(i,0),i) )
        return '\n'.join(ctstring)
    dot2ct=staticmethod(dot2ct)
    def draw(fs, ftype = 'ps'):
        '''
        Draw structure in Postscript or SVG format.
        Input: FastS object
            >YIL140W
            AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
            .(((.......)))............(((...)))...............    (-3.1)
            ..............(((.........(((...)))......)))......    (-3.0)
            ...........((((............))))...................    (-3.0)
        Usage:
            Utils.draw(fs, 'ps')
        Parameters:
            fs: FastS object
                FastS for drawing
            ftype: string
                Options in 'ps', 'svg' and most common figure types.
        Output:
            YIL140W-1.ps  YIL140W-2.ps YIL140W-3.ps
        '''
        ftype = ftype.lower()
        cts = Utils.dot2ct(fs)
        for i in range(len(fs.structures)):
            fprefix = "{0}-{1}".format(fs.name, i+1)
            if ftype == 'svg':
                wRNA.plot(cts,fprefix+".svg", 1, i+1)
            elif ftype == 'ps':
                wRNA.plot(cts,fprefix+".ps", 0, i+1)
            else:
                print >> sys.stderr, "ERROR: file format '{0}' is not supported.".format(ftype)
        return
    draw=staticmethod(draw)
    def mcpu(ncpu,func,paras):
        '''
        Run program in multiple process.
        Parameters:
            ncpu: int
                number of cpus
            func: function
                function to run
            paras: list
                list of parameter list
        Returns:
        '''
        pool = Pool(processes=ncpu)
        for rst in pool.map(func,paras):
            yield rst
        return
    mcpu=staticmethod(mcpu)

class Test(object):
    ''' Test Module. '''
    def testPredictor():
        ''' Test Predictor. '''
        # test FastC class
        T = 37
        threshold = 0
        tfastc = FastC('YIL140W','AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU','..xx......||..........xx...|......||..............')
        print "RNAfold (Vienna RNA package) Test:"
        print "Sequence for prediction:"
        print tfastc
        
        # test RNAfold
        print "\nRNAfold predicted structure:"
        print Predictor.RNAfold(tfastc, T, threshold)
        print 

        # test Fold
        print "\nFold (RNAStructure package) predicted structure:"
        print Predictor.Fold(tfastc,T, threshold)
        print 

        # test mfold
        print "\nmfold predicted structure:"
        print Predictor.mfold(tfastc,T, threshold)
        print 

        # test UNAfold
        print "\nUNAFold predicted structure:"
        print Predictor.UNAFold(tfastc, T, threshold)
        print 

        # test pknots
        print "\npknots predicted structure:"
        print Predictor.pknots(tfastc, T, threshold)
        print

        # test pknotsRG
        print "\npknotsRG predicted structure:"
        print Predictor.pknotsRG(tfastc, T, threshold)
        print
        
        # test ipknot
        print "\nipknot predicted structure:"
        print Predictor.ipknot(tfastc, T, threshold)
        print 

        # test FoldMerge
        print "\n Merge all the predicted structures together."
        print Predictor.FoldMerge(tfastc, T, predictors = ['RNAfold','pknots','Fold','mfold'])
        print
    testPredictor=staticmethod(testPredictor)
    def testUtils():
        ''' Test Utils. '''
        tfastc = FastC('YIL140W','AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU','.|xx.|.x..||.|...xx.x.xx...||.x|..||...|..........')
        tfasts = Predictor.RNAfold(tfastc, 37, 0)
        print "DOT format to CT format:"
        cts = Utils.dot2ct(tfasts)
        print cts
        print
        print "CT format back to DOT format:"
        print Utils.ct2dot(cts)
    testUtils=staticmethod(testUtils)

class ThreadSafeFile(object):
    ''' Lock file while written from thread. '''
    def __init__(self,f):
        ''' Initiation. '''
        self.fh = f
        self.lock = Lock()
        self.locked = True
    def write(self,lstr):
        ''' Write string to file. '''
        with self.lock:
            self.fh.write(lstr)
            self.fh.flush()
    def __enter__(self):
        ''' On enter. '''
        return self
    def __exit__(self,etype,value,traceback):
        ''' On exit. '''
        pass

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv) == 1:
        sys.exit('''
        Module description:
            RNA Structure Prediction Module.
        Test:
            python wRNA.py test
        
        Usage:
            import wRNA
            tfastc = wRNA.FastC('YIL140W','AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU','.|xx.|.x..||.|...xx.x.xx...||.x|..||...|..........')
            print tfastc
            print wRNA.Predictor.RNAfold(tfastc)

        Output:
            >YIL140W
            AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
            .|xx.|.x..||.|...xx.x.xx...||.x|..||...|..........
            >YIL140W
            AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
            ((...(..((((.((............)))..))))...)).........    (9.70)
        ''')
    else:
        Test.testPredictor()
        #Test.testUtils()

