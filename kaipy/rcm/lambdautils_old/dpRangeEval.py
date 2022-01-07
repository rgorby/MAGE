import numpy as np
import argparse
from argparse import RawTextHelpFormatter

from kaipy.rcm.lambdautils.AlamData import AlamData
import kaipy.rcm.lambdautils.dpetadp as dpetadp
import kaipy.rcm.lambdautils.plotter as plotter

#Constants
boltz = 1.38e-23  # boltzmann constant
ev = 1.6e-19  # electron volt
nt = 1.0E-9
m_e = 9.109E-31
m_p = 1.673E-27
radius_earth_m = 6.371E6
pressure_factor = 2./3.*ev/radius_earth_m*nt
density_factor = nt/radius_earth_m

#Input defaults
f_in = "rcmconfig.h5"
tag = ""
tiote = 4.0
L_in = 10
kap_in = -1
doShowPlot = False

def getEtas(alamData, dens, press, vm, kap):
    if kap == kap_in:
        etas = dpetadp.dp2eta(dens, press, vm, tiote, alamData)
    else:
        etas = dpetadp.dp2eta_kap_runov(dens, press, vm, tiote, kap, alamData)
    return etas

"""
    Calculate how many channels contribute (1-epsilon)*P pressure out of total pressure P
    Highest channels contributing epsilon pressure can be ignored by RCM advector

    **This implementation keeps all electron channels and only removes ion channels
      Need to check how okay this is to do
"""
def getFracChannelUse(alamData, dens_arr, press_arr, vm, kap, epsilon=1E-6):
    fracChannelUse = np.zeros((len(press_arr)))  # Fraction of total proton channel use

    alams_p = alamData.alams['spec2']

    for i in range(len(press_arr)):
        etas = getEtas(alamData, dens_arr[i], press_arr[i], vm, kap)
        dTotal, pTotal = dpetadp.eta2dp(etas, alamData, vm)

        numK = len(etas['spec2'])
        k_cutoff = numK-1
        eps_pressure = 0  # Epsilon pressure
        doContinue = True
        while k_cutoff > 0:
            press_new = pressure_factor*np.abs(alams_p[k_cutoff])*etas['spec2'][k_cutoff]*vm**2.5
            if press_new + eps_pressure < epsilon*pTotal:  # Safe to include this channel to those that will be cut off
                k_cutoff -= 1
                eps_pressure += press_new
            else:
                break

        fracChannelUse[i] = (k_cutoff+1)/numK
    return fracChannelUse


def getRatios1D(alamData, dens_arr, press_arr, vm, kap):
    dRatios = np.zeros((len(press_arr)))
    pRatios = np.zeros((len(press_arr)))

    for i in range(len(press_arr)):
        etas = getEtas(alamData, dens_arr[i], press_arr[i], vm, kap)
        dRatios[i], pRatios[i] = dpetadp.eta2dp(etas, alamData, vm)
        dRatios[i] /= dens_arr[i]
        pRatios[i] /= press_arr[i]
    return dRatios, pRatios

def getRatios2D(alamData, dens_arr, press_arr, vm):

    dRatios = np.zeros((len(press_arr), len(dens_arr)))
    pRatios = np.zeros((len(press_arr), len(dens_arr)))
    
    for i in range(len(press_arr)):
        for j in range(len(dens_arr)):
            etas = dpetadp.dp2eta(dens_arr[j], press_arr[i], vm, tiote, alamData)
            dRatios[i,j], pRatios[i,j] = dpetadp.eta2dp(etas, alamData, vm)
            dRatios[i,j] /= dens_arr[j]
            pRatios[i,j] /= press_arr[i]
    return dRatios, pRatios

if __name__=="__main__":
    
    MainS = """Evaluate DP->eta->DP accuracy for a given alam distribution"""

    
    parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
    parser.add_argument('--file', type=str, metavar="<eta file>",default=f_in,
            help="Name of file containing alam channels (default: %(default)s)")
    parser.add_argument('--tiote', type=float,default=tiote,
            help="T_ion/T_electron (default: %(default)s)")
    parser.add_argument('--L', type=float, metavar="[Re]",default=L_in,
            help="L shell in [nPa] (default: %(default)s)")
    parser.add_argument('--tag', type=str,default=tag,
            help="Tag for plot title (default: %(default)s)")
    parser.add_argument('--kap', type=float, metavar="K", default=kap_in,
            help="Kappa value. If unspecified, Maxwellian is used.")
    parser.add_argument("--show", default=False, action='store_true',
            help="Show superplot after saving (default: %(default)s)")

    args = parser.parse_args()
    fname = args.file
    tiote = args.tiote
    tag = args.tag
    Lshell = args.L
    kap = args.kap
    doShowPlot = args.show

    bVol = 0.34 #* 6.371E15 # [Re/nT] * 1E9 * 6.371E6
    vm = bVol**(-2./3.)
    bVol = dpetadp.L_to_bVol(Lshell)
    vm = bVol**(-2./3.)

    alamData = AlamData(fname)

    nSamples = 100
    dMin = 0.01  * 1E6
    dMax = 500   * 1E6
    pMin = 0.01  * 1E-9
    pMax = 1000   * 1E-9
    doLog = True
    
    if(doLog):
        dens_arr = np.logspace(np.log10(dMin), np.log10(dMax), nSamples)
        press_arr = np.logspace(np.log10(pMin), np.log10(pMax), nSamples)
    else:
        dens_arr = np.linspace(dMin, dMax, nSamples)
        press_arr = np.linspace(pMin, pMax, nSamples)


    #dRatios, pRatios = getRatios2D(alamData, dens_arr, press_arr, vm)
    #plotter.plotRatios_dpRange2D(dens_arr, press_arr, dRatios, pRatios, tag)

    # Do range of temperatures instead
    ktmin = 1E-3 # [keV]
    ktmax = 1000 # [keV]
    dconst = 10 # [1/cc]
    kt_arr = np.logspace(np.log10(ktmin), np.log10(ktmax), nSamples)
    dens_arr = np.ones((nSamples))*dconst*1E6 # [1/m^3]
    press_arr = np.array([kt*dconst/6.25 * 1E-9 for kt in kt_arr]) # [Pa]

    #print(press_arr)

    dRatios, pRatios = getRatios1D(alamData, dens_arr, press_arr, vm, kap)
    fracChannelUse = getFracChannelUse(alamData, dens_arr, press_arr, vm, kap)

    title = tag
    if kap == kap_in:
        title += "    Maxwellian"
    else:
        title += "    Kappa = " + str(kap)
    title += "\nL = {}  bVol = {:4.2e}".format(Lshell, bVol)
    title += "\nd = {} [1/cc], p = ({:1.2e}, {:1.2e}) [nPa]".format(dconst, press_arr[0]*1E9, press_arr[-1]*1E9)
    plotter.plotRatios_tempRange1D(kt_arr, dRatios, pRatios, title, fracChannelUse, doShow=doShowPlot)
