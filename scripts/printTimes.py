#!/usr/bin/env python

import sys
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from   astropy.time import Time

import kaipy.kaiH5 as kaiH5
import kaipy.remix.remix as remix

ftag = 'wsa_0004_0002_0004_0000_0000_0000'

MainS = """Description"""

parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")

args = parser.parse_args()

wsaFile = args.id+'.gam.h5'

nsteps,sIds=kaiH5.cntSteps(wsaFile)

T=kaiH5.getTs(wsaFile,sIds,aID='MJD')

for i,tt in enumerate(T):
	print('Step#%06d: '%sorted(sIds)[i],Time(tt,format='mjd').iso)
sys.exit(0)
