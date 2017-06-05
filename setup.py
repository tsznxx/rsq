#!/usr/bin/python
#Last-modified: 24 Mar 2015 01:56:04 PM

#         Module/Scripts Description
# 
# Copyright (c) 2008 Yunfei Wang <Yunfei.Wang1@utdallas.edu>
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @status:  experimental
# @version: 1.0.0
# @author:  Yunfei Wang
# @contact: Yunfei.Wang1@utdallas.edu

# ------------------------------------
# python modules
# ------------------------------------

import os,sys
import string
from setuptools import setup, find_packages, Extension

# ------------------------------------
# constants
# ------------------------------------

EXTERNAL = ['RNAlib']

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__ == '__main__':
    if sys.version[:3] != '2.7':
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)

    # includepy = "%s/include/python%s" % (sys.prefix, sys.version[:3])
    with open("README",'r') as fh:
        long_description = fh.read()

    # PROG and VERSION
    with open('VERSION','r') as fh:
        PROG, VERSION = fh.next().split()

    # Compile external lib
    if 'clean' in sys.argv:
        for elib in EXTERNAL:
            print >>sys.stderr, "Clean {0} Package ...".format(elib)
            os.system('cd external/{0} && make clean && cd ../..'.format(elib))
        print >>sys.stderr, "Clean dist and egg info ..."
        os.system('if [ -d dist ]; then rm -rf dist; fi')
        os.system('if [ -f {0}.egg-info ]; then rm {1}.egg-info; fi'.format(PROG,PROG))
        os.system('if [ -d {0}.egg-info ]; then rm -rf {1}.egg-info; fi'.format(PROG,PROG))
    else:
        for elib in EXTERNAL:
            print >>sys.stderr, "Compile {0} Package ...".format(elib)
            os.system('cd external/{0} && make && cd ../..'.format(elib))
    
    # install requirement
    install_requires = [ "fisher >= 0.1.4",
                         "ngslib>=1.1.10",
                         "numpy >= 1.4.1"]

    setup(name=PROG,
          version=VERSION,
          author='Yunfei Wang',
          author_email='yfwang0405@gmail.com',
          url='http://http://rsqwiki.appspot.com',
          license="GNU General Public License (GPL)",
          keywords = "Python, RNA structure prediction, fold",
          description = ("Python Package for RNA structurome quantification analysis."),
          long_description = long_description,
          package_dir={PROG:'src'},
          packages = [PROG],
          scripts=['bin/'+PROG,
                   'bin/wProfileExtraction.py',
                   'bin/wCleanTranscriptome.py',
                   'bin/wTransGrouping.py'],
          ext_modules=[Extension('wRNA',['external/RNAlib/wRNA/wRNA.cpp'],extra_link_args=['-lm','external/RNAlib/libRNA.a'],extra_compile_args='-w -shared -fPIC -p -Iexternal/RNAlib/fold -Iexternal/RNAlib/plot -Iexternal/RNAlib/wRNA'.split(' '))],
          classifiers=['Environment :: Console',
                       'Development Status :: 3 - Alpha',
                       'Intended Audience :: Developers',
                       'License :: OSI Approved :: GNU General Public License (GPL)',
                       'License :: Free for non-commercial use',
                       'Operating System :: Unix',
                       'Programming Language :: Python :: 2.7',
                       'Topic :: Scientific/Engineering :: Bio-Informatics'],
          install_requires=install_requires)

