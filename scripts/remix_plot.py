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
from astropy.time import Time
import kaipy.kaiH5 as kaiH5

parser = argparse.ArgumentParser()
parser.add_argument('remixFile',help='REMIX file to use')
parser.add_argument('-UT',"--UniversalTime",help="UT to plot in the format YYYY:MM:DDThh:mm:ss'")
#Finalize parsing
args = parser.parse_args()
 
nsteps,sIds=kaiH5.cntSteps(args.remixFile)
T=kaiH5.getTs(args.remixFile,sIds,aID='MJD')

if not(args.UniversalTime):
    for i,tt in enumerate(T):
        print('Step#%06d: '%sorted(sIds)[i],Time(tt,format='mjd').iso)
    sys.exit(0)
else:
    t0 = Time(args.UniversalTime)
    if (t0.mjd<T.min()) or (t0.mjd>T.max()):
        sys.exit('Time outside bounds. Stopping. ')
        
    # find closest time
    imin = np.argmin(np.abs(t0.mjd-T))
    print('Found closest time:',Time(T[imin],format='mjd').iso)
################################################################

# now plotting

import kaipy.remix.remix as remix
import sys
from matplotlib.pyplot import rc,figure,figtext,subplot,show,pcolormesh
from numpy import arctan2,sqrt,pi

rc('mathtext',fontset='stixsans',default='regular')
rc('font',size=11)

step = sorted(sIds)[imin]
ion = remix.get_data(args.remixFile,step)

x = ion['X']
y = ion['Y']

for h in ['NORTH','SOUTH']:
    figure(figsize=(10,6))
    figtext(0.5,0.92,'MIX ('+h+')\n'+Time(T[imin],format='mjd').iso,
            fontsize=14,multialignment='center')

    theta=arctan2(y,x)
    theta[theta<0]=theta[theta<0]+2*pi
    r=sqrt(x**2+y**2)

    if (h.lower()=='north'):
        variables['potential']['data'] = ion['Potential '+h]
        variables['current']['data']   = ion['Field-aligned current '+h]
        variables['sigmap']['data']    = ion['Pedersen conductance '+h]
        variables['sigmah']['data']    = ion['Hall conductance '+h]
        variables['energy']['data']    = ion['Average energy '+h]
        variables['flux']['data']      = ion['Number flux '+h]
        # variables['efield']['data']    = efield_n*1.e6
        # variables['joule']['data']     = sigmap_n*efield_n**2*1.e-3
    else:
        variables['potential']['data'] = ion['Potential '+h][:,::-1]
        variables['current']['data']   = ion['Field-aligned current '+h][:,::-1]
        variables['sigmap']['data']    = ion['Pedersen conductance '+h][:,::-1]
        variables['sigmah']['data']    = ion['Hall conductance '+h][:,::-1]
        variables['energy']['data']    = ion['Average energy '+h][:,::-1]
        variables['flux']['data']      = ion['Number flux '+h][:,::-1]
    
    subplot(231,polar=True)
    remix.plot(theta,r,variables,'potential')
    
    subplot(234,polar=True)
    remix.plot(theta,r,variables,'current')

    subplot(232,polar=True)
    remix.plot(theta,r,variables,'sigmap')

    subplot(235,polar=True)
    remix.plot(theta,r,variables,'sigmah')

    subplot(233,polar=True)
    remix.plot(theta,r,variables,'energy')

    subplot(236,polar=True)
    remix.plot(theta,r,variables,'flux')

show()


