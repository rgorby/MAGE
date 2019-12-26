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
nStp = -1
ftag = "msphere"
doNflux = False
printAll = False

MainS = """Creates simple multi-panel REMIX figure for a GAMERA magnetosphere run.
Top Row - FAC (with potential contours overplotted), Pedersen and Hall Conductances
Bottom Row - Joule heating rate, particle energy and energy flux
"""

parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
parser.add_argument('-n' ,type=int,metavar="step" ,default=nStp,help="Time slice to plot, similar to msphpic.py (default: %(default)s)")
parser.add_argument('-print', action='store_true', default=printAll,help="Print list of all steps and time labels (default: %(default)s)")
parser.add_argument('-nflux', action='store_true', default=doNflux,help="Show number flux instead of energy flux (default: %(default)s)")

# also, optional min/max values for plotting
# it's ugly to specify them in the command line
# but I don't know of a better way to expose them to the user
# config file wouldn't work because we'd need to provide a sample somewhere
# on second thought, this is too ugly -- don't want to do it. Leaving here just in case.
# parser.add_argument('-curMin' ,type=float,metavar="current minimum",help="FAC minimum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-curMax' ,type=float,metavar="current maximum",help="FAC maximum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-sigpMin' ,type=float,metavar="sigmaP minimum",help="Pedersen conductance minimum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-sigpMax' ,type=float,metavar="sigmaP maximum",help="Pedersen conductance maximum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-sighMin' ,type=float,metavar="sigmaH minimum",help="Hall conductance minimum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-sighMax' ,type=float,metavar="sigmaH maximum",help="Hall conductance maximum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-jouleMin' ,type=float,metavar="joule minimum",help="Joule heating minimum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-jouleMax' ,type=float,metavar="joule maximum",help="Joule heating maximum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-energyMin' ,type=float,metavar="energy minimum",help="Particle energy minimum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-energyMax' ,type=float,metavar="energy maximum",help="Particle energy maximum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-efluxMin' ,type=float,metavar="eflux minimum",help="Particle energy flux minimum for plotting. If not specified, defaults in the remix class are used.")
# parser.add_argument('-efluxMax' ,type=float,metavar="eflux maximum",help="Particle energy flux maximum for plotting. If not specified, defaults in the remix class are used.")

#Finalize parsing
args = parser.parse_args()

remixFile = args.id+'.mix.h5'
 
nsteps,sIds=kaiH5.cntSteps(remixFile)
T=kaiH5.getTs(remixFile,sIds,aID='MJD')

if args.print:    # if no arguments provided
	for i,tt in enumerate(T):
		print('Step#%06d: '%sorted(sIds)[i],Time(tt,format='mjd').iso)
	sys.exit(0)
else:
	if args.n == -1: args.n = sorted(sIds)[-1]    # take last step by default
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
	
	gs = gridspec.GridSpec(2,3,figure=fig,left=0.03,right=0.97, top=0.9,bottom=0.03)

	ion.plot('current'  ,gs=gs[0,0])
	ion.plot('sigmap'   ,gs=gs[0,1])
	ion.plot('sigmah'   ,gs=gs[0,2])		
	ion.plot('joule'    ,gs=gs[1,0])
	ion.plot('energy'   ,gs=gs[1,1])
	if args.nflux:
		ion.plot('flux' ,gs=gs[1,2])
	else:
		ion.plot('eflux',gs=gs[1,2])

	if (h.lower()=='north'):
		plt.savefig('remix_n.png',dpi=300)
	else:
		plt.savefig('remix_s.png',dpi=300)


