#!/usr/bin/env python

# DEFINE DATA LIMITS
#hemisphere = 'north'
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
import kaipy.remix.remix as remix
import sys
from matplotlib.pyplot import rc,figure,figtext,subplot,show,pcolormesh
from numpy import arctan2,sqrt,pi

rc('mathtext',fontset='stixsans',default='regular')
rc('font',size=11)

ion = remix.get_data(sys.argv[1])

for h in ion.keys():
    if h=='Sim time': continue

    figure(figsize=(10,6))
    figtext(0.5,0.92,'MIX ('+h+')\n%4d:%02d:%02d  %02d:%02d:%02d' % 
            (ion['Sim time'][0],ion['Sim time'][1],ion['Sim time'][2],ion['Sim time'][3],ion['Sim time'][4],ion['Sim time'][5]),
            fontsize=14,multialignment='center')

    x = ion[h]['x']
    y = ion[h]['y']
    theta=arctan2(y,x)
    theta[theta<0]=theta[theta<0]+2*pi
    r=sqrt(x**2+y**2)

    if (h.lower()=='north'):
        variables['potential']['data'] = ion[h]['psi']
        variables['current']['data']   = ion[h]['fac']
        variables['sigmap']['data']    = ion[h]['sigmap']
        variables['sigmah']['data']    = ion[h]['sigmah']
        variables['energy']['data']    = ion[h]['energy']
        variables['flux']['data']      = ion[h]['flux']
        # variables['efield']['data']    = efield_n*1.e6
        # variables['joule']['data']     = sigmap_n*efield_n**2*1.e-3
    else:
        variables['potential']['data'] = ion[h]['psi'][:,::-1]
        variables['current']['data']   = ion[h]['fac'][:,::-1]
        variables['sigmap']['data']    = ion[h]['sigmap'][:,::-1]
        variables['sigmah']['data']    = ion[h]['sigmah'][:,::-1]
        variables['energy']['data']    = ion[h]['energy'][:,::-1]
        variables['flux']['data']      = ion[h]['flux'][:,::-1]
    
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

