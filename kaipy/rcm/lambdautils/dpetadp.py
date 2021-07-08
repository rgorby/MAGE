import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from math import gamma

from kaipy.rcm.lambdautils.AlamData import AlamData
from kaipy.rcm.lambdautils.helperdefs import L_to_bVol
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

isMain = False

def erfexpdiff(A, xplus, xminus):
    
    #Difference of erf's using Abramowitz & Stegun, 7.1.26
    p  =  0.3275911  
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429

    tp = 1.0/(1+p*xplus)
    tm = 1.0/(1+p*xminus)
    ep = np.exp(-(xplus**2))
    em = np.exp(-(xminus**2))

    erfdiff = - (a1*tp + a2*(tp**2.0) + a3*(tp**3.0) + a4*(tp**4.0) + a5*(tp**5.0))*ep \
              + (a1*tm + a2*(tm**2.0) + a3*(tm**3.0) + a4*(tm**4.0) + a5*(tm**5.0))*em
    expdiff = 2.0*(xplus*ep - xminus*em)/np.sqrt(np.pi)
    #print("efrfdiff:{:1.2e}\texpdiff:{:1.2e}".format(erfdiff, expdiff))
    eta = A*(erfdiff-expdiff)
    return eta

def dp2eta(dens, press, vm, tiote, alamData): # dens = [1/m^3], press = [Pa]
    pion = press*tiote/(1.0+tiote)
    pele = press*  1.0/(1.0+tiote)

    ti = pion/dens/boltz
    te = pele/dens/boltz

    if isMain:
        print("Ti={:2.2e} [eV]\tTe={:2.2e} [eV]".format(ti*boltz/1.6E-19,te*boltz/1.6E-19))

    A0 = (dens/density_factor)/vm**1.5

    etadict = {}

    for s in alamData.alams.keys():  # s = species
        alams = alamData.alams[s]
        amins = alamData.amins[s]
        amaxs = alamData.amaxs[s]

        tk = ti if alams[1] > 0 else te
        
        eta_arr = np.array([])
        for k in range(len(alams)):
            xplus = np.sqrt(ev*np.abs(amaxs[k])*vm/boltz/tk)
            xminus = np.sqrt(ev*np.abs(amins[k])*vm/boltz/tk)

            eta = erfexpdiff(A0, xplus, xminus)
            eta_arr = np.append(eta_arr, eta)
        etadict[s] = eta_arr

    return etadict

#  Python implementation of Kareem's code (which is an implementation of original kappa notes)
def dp2eta_kap_kfort(dens, press, vm, tiote, kap, alamData):
    pion = press*tiote/(1.0+tiote)
    pele = press*  1.0/(1.0+tiote)

    tk_i = pion/dens
    tk_e = pele/dens

    if isMain:
        print("Ti={:2.2e} [eV]\tTe={:2.2e} [eV]".format(tk_i/1.6E-19,tk_e/1.6E-19))

    kap15 = kap-1.5
    kapgam = gamma(kap+1.0)/gamma(kap-0.5)
    kapbar = kap15

    etadict = {}

    for s in alamData.alams.keys():
        alams = alamData.alams[s]
        amins = alamData.amins[s]
        amaxs = alamData.amaxs[s]

        tk = tk_i if alams[1] > 0 else tk_e

        Tev = tk/1.6E-19
        E0_ev = Tev*kap15/kap
        A0 = (2.0/np.sqrt(np.pi))*(dens/density_factor)/vm**1.5

        eta_arr = np.array([])
        for k in range(len(alams)):
            E_ev = abs(alams[k])*vm  # [eV]
            kArg = 1 + (E_ev/E0_ev)/kapbar
            delscl = (amaxs[k]-amins[k])/E0_ev
            etak = A0*kapgam/kapbar**1.5 * np.sqrt(E_ev/E0_ev)*delscl*kArg**(-kap-1)
            
            eta_arr = np.append(eta_arr, etak)
        etadict[s] = eta_arr

    return etadict

# Runov et al. (2015) pg. 4377 10.1002/2015JA021166
def dp2eta_kap_runov(dens, press, vm, tiote, kap, alamData):
    pion = press*tiote/(1.0+tiote)
    pele = press*  1.0/(1.0+tiote)

    tk_i = pion/dens
    tk_e = pele/dens

    if isMain:
        print("Ti={:2.2e} [eV]\tTe={:2.2e} [eV]".format(tk_i/1.6E-19,tk_e/1.6E-19))

    kap15 = kap-1.5
    kapgam = gamma(kap+1.0)/gamma(kap-0.5)

    etadict = {}

    for s in alamData.alams.keys():
        alams = alamData.alams[s]
        amins = alamData.amins[s]
        amaxs = alamData.amaxs[s]

        tk = tk_i if alams[1] > 0 else tk_e
        m_s = m_p if alams[1] > 0 else m_e
        Tev = tk/1.6E-19
        E0_ev = Tev*kap15/kap
        #E0_ev = Tev
        A0 = (2.0/np.sqrt(np.pi))*(dens/density_factor)/vm**1.5
        #A0 = (dens/density_factor)/vm**1.5

        eta_arr = np.array([])
        
        for k in range(len(alams)):
            E_ev = abs(alams[k])*vm  # [eV]
            kArg = 1 + E_ev/(kap*E0_ev)
            delscl = (amaxs[k]-amins[k])*vm/E0_ev

            etak = A0*np.sqrt(E_ev/E0_ev) * delscl / (kap)**1.5 * kapgam * kArg**(-kap-1)

            #etak = A0*kapgam/kapbar**1.5 * np.sqrt(E_ev/E0_ev)*delscl*kArg**(-kap-1)
            
            eta_arr = np.append(eta_arr, etak)
        """
        for k in range(len(alams)):
            E_ev_min = abs(amins[k])*vm
            E_ev_max = abs(amaxs[k])*vm
            kArg_min = (1+E_ev_min/(kap*E0_ev))
            kArg_max = (1+E_ev_max/(kap*E0_ev))

            consts = A0 /(kap*E0_ev)**1.5 * kapgam

            etak = (amaxs[k]-amins[k])*vm/2 * consts*((np.sqrt(E_ev_max)*kArg_max**(-kap-1)) + np.sqrt(E_ev_min)*kArg_min**(-kap-1))

            eta_arr = np.append(eta_arr, etak)
        """
        etadict[s] = eta_arr

    return etadict

def etaLinearRescaler(alamData, etaData, targetN, targetP, vBol):
    vm = vBol**(-2/3)
    tN = targetN/density_factor#/vm**1.5
    tP = targetP/pressure_factor#/vm**2.5

    newEtaData = etaData

    for s in alamData.alams.keys():
        alams = alamData.alams[s]
        etas = etaData[s]

        #Do the summations
        
        kMax = len(alams)-1
        doIterate = True
        while doIterate:
            S1 = 0  # Sum eta_k
            S2 = 0  # Sum eta_k*abs(alam_k)
            S3 = 0  # Sum eta_k*abs(alam_k)^2
            for k in range(kMax+1):
                S1 += etas[k]
                S2 += etas[k]*np.abs(alams[k])
                S3 += etas[k]*alams[k]**2

            denom = 1/(S1*S3 - S2**2)
            a = denom*(S3*tN*bVol - S2*1.5*tP*bVol**(5/3))
            b = denom*(-S2*tN*bVol + S1*1.5*tP*bVol**(5/3))
            print("a, b = {},{}".format(a, b))
            #TODO: set new kMax and iterate if necessary
            if (a+b*np.abs(alams[kMax]))*etas[kMax] < 0:
                etas[kMax] = 0
                kMax -= 1
            else:
                doIterate = False
        newEtaData[s] = (a+b*np.abs(alams))*etas
    
    return newEtaData

def eta2dp(etaData, alamData, vm):

    dTot = 0
    pTot = 0

    for s in alamData.alams.keys():  # s = species
        alams = alamData.alams[s]
        etas = etaData[s]

        dSpec, pSpec = eta2dp_s(etas, alams, vm)
        if isMain:
            print("{}: d = {:1.2e} [1/cc]   p = {:1.2e} [nPa]".format(s, dSpec*1E-6, pSpec*1E9))
        pTot += pSpec
        if "1" not in s:
            dTot += dSpec
    if isMain:
        print("Total: d = {:1.2e} [1/cc]   p = {:1.2e} [nPa]".format(dSpec*1E-6, pSpec*1E9))

    return dTot, pTot

def eta2dp_s(etas_s, alams_s, vm, k1=-1, k2=-1):  # etas_s, alams_s are of a specific species
    if k1 == -1:
        k1 = 0
    if k2 == -1:
        k2 = len(etas_s)-1

    dens = 0  # [1/m^2]
    press = 0  # [Pa]

    for k in range(k1,k2+1):
        dens += density_factor*etas_s[k]*vm**1.5
        press += pressure_factor*np.abs(alams_s[k])*etas_s[k]*vm**2.5
    return dens, press


if __name__=="__main__":
    isMain = True
    MainS = """Calculate DP => RCM etas => DP for a given alam distribution"""

    f_in = "rcmconfig.h5"
    dens_in = 0.5
    press_in = 0.85
    tiote_in = 4.0
    L_in = 10
    kap_in = -1
    tag = ""
    parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
    parser.add_argument('--file', type=str, metavar="<eta file>",default=f_in,
            help="Name of file containing alam channels (default: %(default)s)")
    parser.add_argument('--d', type=float, metavar="[1/cc]",default=dens_in,
            help="Density in [1/cc] (default: %(default)s)")
    parser.add_argument('--p', type=float, metavar="[nPa]",default=dens_in,
            help="Pressure in [nPa] (default: %(default)s)")
    parser.add_argument('--L', type=float, metavar="[Re]",default=L_in,
            help="L shell in [nPa] (default: %(default)s)")
    parser.add_argument('--tiote', type=float,default=tiote_in,
            help="T_ion/T_electron (default: %(default)s)")
    parser.add_argument('--tag', type=str,default=tag,
            help="Tag for plot title (default: %(default)s)")
    parser.add_argument('--kap', type=float, metavar="K", default=kap_in,
            help="Kappa value. If unspecified, Maxwellian is used.")
    parser.add_argument("--show", default=False, action='store_true',
            help="Show superplot after saving (default: %(default)s)")
    args = parser.parse_args()
    fname = args.file
    dens = args.d*1E6  # 1/cc to 1/m^3
    press = args.p*1E-9  # nPa to Pa
    Lshell = args.L
    tiote = args.tiote
    tag = args.tag
    kap = args.kap
    doShowPlot = args.show

    bVol = 0.34 #* 6.371E15 # [Re/nT] * 1E9 * 6.371E6
    bVol = L_to_bVol(Lshell)
    vm = bVol**(-2./3.)
    print("bVol = {}, vm = {}".format(bVol, vm))
    alamData = AlamData(filename=fname)
    if kap == -1:
        etaData = dp2eta(dens, press, vm, tiote, alamData)
    else:
        #etaData = dp2eta_kap_runov(dens, press, vm, tiote, kap, alamData)
        etaData = dp2eta_kap_runov(dens, press, vm, tiote, kap, alamData)

    #etaLinearRescaler(alamData, etaData, dens, press, bVol)

    dTot, pTot = eta2dp(etaData, alamData, vm)  # [1/m^3] and [Pa]
    dNew_o_dOrig = dTot/dens
    pNew_o_pOrig = pTot/press

    dens_cc = dens*1E-6
    press_npa = press*1E9
    kT = 6.25*press_npa/dens_cc  # [keV]

    title = tag
    if kap == kap_in:
        title += "    Maxwellian"
    else:
        title += "    Kappa = " + str(kap)
    title += "\nd = {:4.2f} [1/cc]  p = {:4.2f} [nPa]  kT = {:1.2e} [keV]".format(dens_cc, press_npa, kT)
    title += "\nL = {}  bVol = {:4.2f}".format(Lshell, bVol)
    title += "\nD'/D = {:1.2e}  P'/P = {:1.2e}".format(dNew_o_dOrig, pNew_o_pOrig)
    plotter.plotEtas(alamData, etaData, vm, title, doShow=doShowPlot)

    


    print("p_new/p_0: {:0.2e}".format(pNew_o_pOrig))
    print("d_new/d_0: {:0.2e}".format(dNew_o_dOrig))
