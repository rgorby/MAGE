#!/usr/bin/env python

################ first figure out the time ################

import sys
import argparse
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from   astropy.time import Time

import kaipy.kaiH5 as kaiH5
import kaipy.remix.remix as remix

parser = argparse.ArgumentParser()
parser.add_argument('remixFile',help='REMIX file to use')
parser.add_argument('-n',type=int,help="Time slice to plot (similar to msphpic.py)")

#Finalize parsing
args = parser.parse_args()
 
nsteps,sIds=kaiH5.cntSteps(args.remixFile)
T=kaiH5.getTs(args.remixFile,sIds,aID='MJD')

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
ion = remix.remix(args.remixFile,args.n)

for h in ['NORTH','SOUTH']:
	fig = plt.figure(figsize=(12,7.5))
	plt.figtext(0.5,0.94,'MIX ('+h+')\n'+Time(T[args.n],format='mjd').iso,
			fontsize=14,multialignment='center',horizontalalignment='center')

	ion.init_vars(h)
	
	gs = gridspec.GridSpec(2,3,figure=fig,left=0.03,right=0.97, top=0.9,bottom=0.03)

	ax1=fig.add_subplot(gs[0,0],polar=True)
	ion.plot('potential',ax=ax1)
	
	ax2=fig.add_subplot(gs[1,0],polar=True)
	ion.plot('current',ax=ax2)

	ax3=fig.add_subplot(gs[0,1],polar=True)
	ion.plot('sigmap',ax=ax3)

	ax4=fig.add_subplot(gs[1,1],polar=True)
	ion.plot('sigmah',ax=ax4)

	ax5=fig.add_subplot(gs[0,2],polar=True)
	ion.plot('energy',ax=ax5)

	ax6=fig.add_subplot(gs[1,2],polar=True)
	ion.plot('flux',ax=ax6)

	if (h.lower()=='north'):
		plt.savefig('remix_n.png',dpi=300)
	else:
		plt.savefig('remix_s.png',dpi=300)


