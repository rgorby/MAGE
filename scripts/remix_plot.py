#!/usr/bin/env python

# DEFINE DATA LIMITS
variables = { 'potential' : {'min':-100,
							'max': 100},
			 'current'   : {'min':-1,
							'max':1},
			 'sigmap'    : {'min':1,
							'max':10},
			 'sigmah'    : {'min':2,
							'max':20},
			 'energy'    : {'min':0,
							'max':100},
			 'flux'      : {'min':0,
							'max':3.e8},
			 'efield'    : {'min':-1,
							'max':1},
			  'joule'    : {'min':0,
							'max':50},
			 }
#ncontours = 101
#nticks = 11  # how many ticks on the colorbar

## ====== SHOULD NOT NEED TO CHANGE ANYTHING BELOW ==============================

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

ion = remix.get_data(args.remixFile,args.n)

x = ion['X']
y = ion['Y']

for h in ['NORTH','SOUTH']:
	fig = plt.figure(figsize=(12,7.5))
	plt.figtext(0.5,0.94,'MIX ('+h+')\n'+Time(T[args.n],format='mjd').iso,
			fontsize=14,multialignment='center',horizontalalignment='center')

	# note, because of how the grid is set up in the h5 file,
	# the code below produces the first theta that's just shy of 2pi.
	# this is because the original grid is staggered at half-cells from data.
	# empirically, this is OK for pcolormesh plots under remix.plot.
	# however, contour plots have issues across the periodic boundary. 
	# still working on it as of 24 Dec 2019
	theta=np.arctan2(y,x)
	theta[theta<0]=theta[theta<0]+2*np.pi
	r=np.sqrt(x**2+y**2)

	if (h.lower()=='north'):
		variables['potential']['data'] = ion['Potential '+h]
		variables['current']['data']   =-ion['Field-aligned current '+h]  # note, converting to common convention (upward=positive)
		variables['sigmap']['data']    = ion['Pedersen conductance '+h]
		variables['sigmah']['data']    = ion['Hall conductance '+h]
		variables['energy']['data']    = ion['Average energy '+h]
		variables['flux']['data']      = ion['Number flux '+h]
		# variables['efield']['data']    = efield_n*1.e6
		# variables['joule']['data']     = sigmap_n*efield_n**2*1.e-3
	else:  # note flipping the y(phi)-axis
		variables['potential']['data'] = ion['Potential '+h][:,::-1]
		variables['current']['data']   = ion['Field-aligned current '+h][:,::-1]
		variables['sigmap']['data']    = ion['Pedersen conductance '+h][:,::-1]
		variables['sigmah']['data']    = ion['Hall conductance '+h][:,::-1]
		variables['energy']['data']    = ion['Average energy '+h][:,::-1]
		variables['flux']['data']      = ion['Number flux '+h][:,::-1]
	
	gs = gridspec.GridSpec(2,3,figure=fig,left=0.03,right=0.97, top=0.9,bottom=0.03)

	ax1=fig.add_subplot(gs[0,0],polar=True)
	remix.plot(theta,r,variables,'potential',ax=ax1)
	
	ax2=fig.add_subplot(gs[1,0],polar=True)
	remix.plot(theta,r,variables,'current',ax=ax2,x=x,y=y)

	ax3=fig.add_subplot(gs[0,1],polar=True)
	remix.plot(theta,r,variables,'sigmap',ax=ax3)

	ax4=fig.add_subplot(gs[1,1],polar=True)
	remix.plot(theta,r,variables,'sigmah',ax=ax4)

	ax5=fig.add_subplot(gs[0,2],polar=True)
	remix.plot(theta,r,variables,'energy',ax=ax5)

	ax6=fig.add_subplot(gs[1,2],polar=True)
	remix.plot(theta,r,variables,'flux',ax=ax6)

	if (h.lower()=='north'):
		plt.savefig('remix_n.png',dpi=300)
	else:
		plt.savefig('remix_s.png',dpi=300)


