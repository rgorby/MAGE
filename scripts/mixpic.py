#!/usr/bin/env python

################ first figure out the time ################

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

# Defaults
nStp = 0
ftag = "msphere"
doNflux = False	

MainS = """Creates simple multi-panel REMIX figure for a GAMERA magnetosphere run.
If run without arguments, prints all steps and UT labels found in the file.
Top Row - TBD
Bottom Row - TBD
"""

parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
parser.add_argument('-n' ,type=int,metavar="step" ,default=nStp,help="Time slice to plot, similar to msphpic.py (default: %(default)s)")
parser.add_argument('-nflux', action='store_true', default=doNflux,help="Show number flux instead of energy flux (default: %(default)s)")

#Finalize parsing
args = parser.parse_args()

remixFile = args.id+'.mix.h5'
 
nsteps,sIds=kaiH5.cntSteps(remixFile)
T=kaiH5.getTs(remixFile,sIds,aID='MJD')

if not(args.n):
	for i,tt in enumerate(T):
		print('Step#%06d: '%sorted(sIds)[i],Time(tt,format='mjd').iso)
	sys.exit(0)
else:
	if args.n not in sIds:
			sys.exit("Time step not in the h5 file.")

	print('Found time:',Time(T[args.n],format='mjd').iso)
################################################################

# Now plotting

mpl.use('Agg')
mpl.rc('mathtext',fontset='stixsans',default='regular')
mpl.rc('font',size=10)

# Initialize the remix class
ion = remix.remix(remixFile,args.n)

# if only plotting one variable, could just do this:
# ion.init_vars(h)
# ion.plot('potential')

for h in ['NORTH','SOUTH']:
	fig = plt.figure(figsize=(12,7.5))
	plt.figtext(0.5,0.94,'MIX ('+h+')\n'+Time(T[args.n],format='mjd').iso,
			fontsize=12,multialignment='center',horizontalalignment='center')

	ion.init_vars(h)
#	ion.efield()
	
	gs = gridspec.GridSpec(2,3,figure=fig,left=0.03,right=0.97, top=0.9,bottom=0.03)

	ion.plot('potential',gs=gs[0,0])
	ion.plot('current'  ,gs=gs[1,0])
	ion.plot('sigmap'   ,gs=gs[0,1])
	ion.plot('sigmah'   ,gs=gs[1,1])
	ion.plot('energy'   ,gs=gs[0,2])
	if args.nflux:
		ion.plot('flux' ,gs=gs[1,2])
	else:
		ion.plot('eflux',gs=gs[1,2])

	if (h.lower()=='north'):
		plt.savefig('remix_n.png',dpi=300)
	else:
		plt.savefig('remix_s.png',dpi=300)


