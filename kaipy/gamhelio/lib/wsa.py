#!/usr/bin/env python

from astropy.io import fits

from numpy import linspace,pi,meshgrid,sin,cos,zeros,ones,dstack,diff,sqrt,array,savetxt,flatnonzero,insert,asarray,zeros_like,argmin,unravel_index

def read(wsa_file,densTempInfile,normalized=False,verbose=True):
    if verbose: info(wsa_file)

    hdul = fits.open(wsa_file)  
    jd_c = hdul[0].header['JULDATE'] #Julian date at the center of the map  
    n_phi_wsa_v = hdul[0].header['NAXIS1']+1  # number of cell vertices
    n_phi_wsa_c = hdul[0].header['NAXIS1']    # number of cell centers
    phi_wsa_v     = linspace(0,360,n_phi_wsa_v)/180.*pi
    phi_wsa_c     = 0.5*(phi_wsa_v[:-1]+phi_wsa_v[1:])

    n_theta_wsa_v = hdul[0].header['NAXIS2']+1  # number of cell vertices
    n_theta_wsa_c = hdul[0].header['NAXIS2']    # number of cell centers
    theta_wsa_v     = linspace(0,180,n_theta_wsa_v)/180.*pi
    theta_wsa_c     = 0.5*(theta_wsa_v[:-1]+theta_wsa_v[1:])

    bi_wsa = hdul[0].data[0,::-1,:]  # note the theta reversal to convert from wsa theta to gamera theta definition
    v_wsa  = hdul[0].data[1,::-1,:]  # note the theta reversal to convert from wsa theta to gamera theta definition
    if densTempInfile:
        n_wsa  = hdul[0].data[2,::-1,:]
        T_wsa  = hdul[0].data[3,::-1,:]
    else:
        n_wsa = 112.64+9.49e7/v_wsa**2
        T0 = 1.44e6 #8e5
        n0 = 300.  
        B0 = bi_wsa[unravel_index(argmin(abs(n_wsa-n0)),n_wsa.shape)] # this is in nT
        n0 = n_wsa[unravel_index(argmin(abs(n_wsa-n0)),n_wsa.shape)]  # this is in cm^-3

#        T_wsa = n0*T0/n_wsa

# The code below allows total pressure conservation in the non-radial
# direction: nT+B**2/8pi = n0*T0+B0**2/8pi. The unit conversion assumes
# B and B0 in nT (which they should be by now) and botzmann constant in the
# denominator.
        T_wsa = n0*T0/n_wsa + (B0**2-bi_wsa**2)/1.38/8./pi*1.e6/n_wsa
# units: B in nT = * 10^-5 G; Kb in Gauss = 1.38*10^-16; n in cm^-3
        P_tot_wsa = n_wsa*1.38e-16 * T_wsa + (bi_wsa**2 * 1.e-10)/8./pi
    hdul.close()
    
    if normalized:
        return(jd_c,phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa,v_wsa,n_wsa*1.67e-24,T_wsa)
    else:
        return(jd_c,phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa*1.e-5,v_wsa*1.e5,n_wsa*1.67e-24,T_wsa)

def info(wsa_file):
    hdul = fits.open(wsa_file)
    print(repr(hdul[0].header))
    hdul.close()

def plot(wsa_file,savefig):
    hdul = fits.open(wsa_file)    
    bi_wsa = hdul[0].data[0,:,:]  
    v_wsa  = hdul[0].data[1,:,:]
    hdul.close()

    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(16,12))
    ax1 = plt.subplot(311)
    p1=ax1.pcolormesh(bi_wsa,cmap='RdBu_r',vmin=bi_wsa.min(),vmax=-bi_wsa.min())
    plt.colorbar(p1,ax=ax1).set_label('Br')
    ax1.set_xlim((0,bi_wsa.shape[1]))
    ax1.set_ylim((0,bi_wsa.shape[0]))

    ax2 = plt.subplot(312,sharex=ax1)
    p2=ax2.pcolormesh(v_wsa)
    plt.colorbar(p2,ax=ax2).set_label('V')
    ax2.set_xlim((0,bi_wsa.shape[1]))
    ax2.set_ylim((0,bi_wsa.shape[0]))

    fig.suptitle(wsa_file)
    if not savefig:
        plt.show()
    else:
        plt.savefig(wsa_file[:-4]+'png')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('wsaFile',help='WSA file to use')
    parser.add_argument('--savefig',help='Drop figure?',default=False,action='store_true')
    args = parser.parse_args()
    
    info(args.wsaFile)
    plot(args.wsaFile,args.savefig)
