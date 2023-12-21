import os
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import kaipy.kdefs as kd
import kaipy.kaiTools as kt
import kaipy.kaiViz as kv
import kaipy.raiju.raijuutils as ru




# random values
rad2deg = 180/np.pi

# Settings
dir = '/Users/sciolam1/Workspace/runs/local/raijudev/test_evol'
raih5fname = 'raijuSA.raiju.h5'

def calc_bVol_ana(colat, b_surf):
    # bVol calc
    cSum =  (35.0      *np.cos(1.0*colat) -      7.0 *np.cos(3.0*colat) \
            +(7.0/5.0) *np.cos(5.0*colat) - (1.0/7.0)*np.cos(7.0*colat))/64.0
    S8 = np.sin(colat)**8.0
    bVol = 2*cSum/S8/b_surf  # [Rp/nT]
    return bVol


def calc_grad_bVol_ana(colat, b_surf, Ri_m):
    # Gradient of bVol calc
    cSum =  (35.0      *np.cos(1.0*colat) -      7.0 *np.cos(3.0*colat) \
            +(7.0/5.0) *np.cos(5.0*colat) - (1.0/7.0)*np.cos(7.0*colat))/64.0
    S8 = np.sin(colat)**8.0

    dSum = (-35.0*np.sin(1.0*colat) + 21.0*np.sin(3.0*colat) \
            - 7.0*np.sin(5.0*colat) +  1.0*np.sin(7.0*colat) )/64.0
    grad_bVol = (2.0/S8/b_surf) * (-8.0/np.tan(colat)*cSum + dSum) / Ri_m  # 1/Ri_m  w.r.t. arc length
    return grad_bVol


def vCorot(Rp_m, Ri_m, cPot, b_surf, r):
    """
    Calculate corotation velocity [m/s] in the ionosphere
    Rp: Planetary radius [m]
    cPot: Corotation potential [kV]
    b_surf: surface B-field [nT]
    r: eval equatorial radius [Rp]

    vCorot [rad/s] = cPot * b_surf / Rp^2
    vCorot [m/s] = cPot * b_surf * r [Rp] / Rp^2 = cPos * b_surf * r / Rp
    """
    v_eq = (cPot*1e3) / (b_surf * 1e-9) * (r*Rp_m) / Rp_m**2  # [m/s]

    sinColat = np.sqrt(1/r)
    r_iono = Ri_m * sinColat

    v_iono = v_eq * r_iono / (r*Rp_m)
    return v_iono


def vGradEq(Rp_m, Ri_m, signQ, b_surf, E, r):
    """
    Calculate gradient-drift velocity of equatorial particle [m/s] mapped to the ionosphere
    signQ: charge sign
    b_surf: surface B-field [nT]
    E: Particle equatorial energy [keV]
    r: radial distance [Rp]

    vd = 1/2* v_perp * r_g * (B x grad-B)/B^2

    1/2 * v_perp * r_g = E/qB
    grad-B/B = -3/r

    vd = -3*E*r^2 / (q * B-surf)
    
    If E in eV, then = -3*E*r^2/B-surf

    """
    sinColat = np.sqrt(1/r)
    #cosColat = np.cos(np.arcsin(sinColat))
    r_iono = Ri_m * sinColat

    #b_iono = b_surf * np.sqrt(1 + 3*sinColat**2)
    v_eq = signQ * -3 * (E * 1e3) * r**2 /(b_surf*1e-9 * rp_m)  # [m/s]
    
    v_iono = v_eq * r_iono / (r*Rp_m)
    return v_iono


def vGCiso(Rp_m: float, Ri_m: float, b_surf: float, alamc: float, r: np.ndarray):
    """
    Calculate gradient-curvature drift of an isotropic distribution in a dipole field
    eq. 8.1.40 from chapter 8 of Wolf's unpublished book :)

    Rp_m: [m]
    Ri_m: [m]
    b_surf: [nT]
    alamc: eV * (Rp/nT)^(2/3)
    r: Radial distance [Rp]

    v_d = B x grad-Wk / qB^2  --> grad-Wk / |B|

    lambda in eV so can drop the q

    grad-Wk = grad(lambda * bVol^-2/3) = lambda*(-2/3*bVol^(-5/3)) * grad(bVol)
    """

    sinColat = np.sqrt(1/r)
    r_iono = Ri_m * sinColat
    colat = np.arcsin(sinColat)

    bVol      = calc_bVol_ana     (colat, b_surf)  # [Rp/nT]
    grad_bVol = calc_grad_bVol_ana(colat, b_surf, Ri_m)  # [Rp/nT/m]
    

    grad_Wk = alamc*(-2./3.*bVol**(-5./3.)) * grad_bVol  # [eV/m]

    # 1/B = [m^2/V/s]
    #Beq = b_surf*1e-9/r**3  # [T]
    #v_drift = -1*grad_Wk / (Beq)  # [m/s]

    cosdip = 2.0*np.cos(colat)/np.sqrt(1.0 + 3.0*np.cos(colat)**2.0) 
    Bmag = (b_surf*1e-9) / (Ri_m/Rp_m)**3 * np.sqrt(1.0+3.0*np.cos(colat)**2.0)

    v_drift = -1*grad_Wk / cosdip / Bmag  # [m/s]
    #v_drift = -1*grad_Wk / cosdip / (b_surf*1e-9)  # [m/s]
    #v_drift = -1*grad_Wk / (b_surf*1e-9)  # [m/s]

    # Map to ionosphere
    #v_iono = v_drift * r_iono / (r*Rp_m)
    #return v_iono

    # We did gradient in iono and used iono b-field, so we already have stuff in ionosphere
    return v_drift

    


def checkBVol(raiInfo: ru.RAIJUInfo):

    with h5.File(raiInfo.fname, 'r') as f5:
        ri_m  = f5['Planet'].attrs['Rad_ionosphere']  # [m]
        rp_m  = f5['Planet'].attrs['Rad_surface'   ]  # [m]
        Bs_nt = f5['Planet'].attrs['Mag Moment'    ] * kd.G2nT  # [nT]
        s5 = f5[raiInfo.stepStrs[1]]
        bVol = ru.getVar(s5, 'bVol')  # [Rp/nT]
        gradVM = ru.getVar(s5, 'gradVM')
        rMin = ru.getVar(s5, 'xmin')

    bVol = bVol[:,0]
    gradVM = gradVM[:,0,0]
    rMin_cc = kt.to_center1D(rMin[:,0])
    sinColat = np.sqrt(1/rMin_cc)
    colat = np.arcsin(sinColat)

    bVol_ana = calc_bVol_ana(colat, Bs_nt)
    grad_bVol_ana = calc_grad_bVol_ana(colat, Bs_nt, ri_m)

    gradVM_ana = (-2/3*bVol**(-5/3)) * grad_bVol_ana

    plt.figure()
    plt.plot   (rMin_cc, bVol_ana , label='ana'  )
    plt.scatter(rMin_cc, bVol, label='model', c='orange')
    plt.title('bVol')
    plt.yscale('log')
    plt.legend()

    plt.figure()
    plt.plot   (rMin_cc, gradVM_ana , label='ana'  )
    plt.scatter(rMin_cc, gradVM, label='model', c='orange')
    plt.title('$\\nabla$ bVol')
    plt.yscale('log')
    plt.legend()


def checkCorot(raiInfo: ru.RAIJUInfo, n: int):
    """
    Check if corotation velocity is correct at timestep n
    For now, just using plasmasphere channel and assuming espot = 0
     If not, espot will affect the total velocity
     TODO: In the future, we can subtract its contribution
    """

    with h5.File(raiInfo.fname, 'r') as f5:
        ri_m  = f5['Planet'].attrs['Rad_ionosphere']  # [m]
        rp_m  = f5['Planet'].attrs['Rad_surface'   ]  # [m]
        Bs_nt = f5['Planet'].attrs['Mag Moment'    ] * kd.G2nT  # [nT]
        cPot  = f5['Planet'].attrs['Psi Corot']  # [kV]
        s5 = f5[raiInfo.stepStrs[n]]
        espot = ru.getVar(s5, 'espot') # [kV]

        if np.any(espot > 0):
            print("Error: espot must be zero for checkCorot to work")
            quit()
        
        rMin = ru.getVar(s5, 'xmin')[:,0]  # Just take all i's at j=0 (y=0)
        cVel_ph = ru.getVar(s5, 'cVel_ph')
        
    rMin_cc = kt.to_center1D(rMin)
    vCorot_ana = vCorot(rp_m, ri_m, cPot, Bs_nt, rMin_cc)
    # Just get velocity from plasmasphere channel
    kPsph = raiInfo.species[ru.flavs_s['PSPH']].kStart
    vCorot_model = cVel_ph[:,0,kPsph]
    plt.figure()
    plt.plot(rMin_cc, vCorot_ana, label='analytic')
    plt.scatter(rMin_cc, vCorot_model, s=10, c='orange', label='model')
    plt.title("Corotation")
    plt.legend()
 

def checkGC(raiInfo: ru.RAIJUInfo, n: int, k: int):
    """
    Check if corotation velocity is correct at timestep n
    For now, just using plasmasphere channel and assuming espot = 0
     If not, espot will affect the total velocity
     TODO: In the future, we can subtract its contribution
    """

    with h5.File(raiInfo.fname, 'r') as f5:
        ri_m  = f5['Planet'].attrs['Rad_ionosphere']  # [m]
        rp_m  = f5['Planet'].attrs['Rad_surface'   ]  # [m]
        Bs_nt = f5['Planet'].attrs['Mag Moment'    ] * kd.G2nT  # [nT]
        cPot  = f5['Planet'].attrs['Psi Corot']  # [kV]
        alamc = f5['alamc'][:]  # [eV *(Rp/nT)**2/3]
        s5 = f5[raiInfo.stepStrs[n]]
        espot = ru.getVar(s5, 'espot') # [kV]
        bVol = ru.getVar(s5, 'bVol')  # [Rp/nT]


        if np.any(espot > 0):
            print("Error: espot must be zero for checkCorot to work")
            quit()
        
        rMin = ru.getVar(s5, 'xmin')
        cVel_ph = ru.getVar(s5, 'cVel_ph')

    # Find species info
    iSpc = 0
    while raiInfo.species[iSpc].kEnd < k: iSpc += 1
    spc = raiInfo.species[iSpc]
    print("Species: ",spc.name)
    signQ = np.sign(spc.q)
    rMin_cc = kt.to_center2D(rMin)
    E = alamc[k]*bVol**(-2/3)*1e-3  # [keV]
    vGC_eq_ana = vGradEq(rp_m, ri_m, signQ, Bs_nt, E, rMin_cc)
    vGC_iso_ana = vGCiso(rp_m, ri_m, Bs_nt, alamc[k], rMin_cc)
    vCorot_ana  = vCorot(rp_m, ri_m, cPot , Bs_nt, rMin_cc)

    vPh_model = cVel_ph[:,:,k]
    vGC_model = vPh_model - vCorot_ana

    plt.figure()
    plt.plot(   rMin_cc[:,0], vGC_eq_ana   [:,0], label='ana (grad_eq)')
    plt.plot(   rMin_cc[:,0], vGC_iso_ana  [:,0], c='green', label='ana (GC_iso)')
    plt.scatter(rMin_cc[:,0], vGC_model[:,0], s=10, c='orange', label='model')
    #plt.plot(rMin_cc[:,0], 1 - vGC_model[:,0] / vD_ana[:,0], label='err')
    #plt.plot(rMin_cc[:,0], 1 - vGC_model[:,0] / vGC_iso_ana[:,0], label='err')
    plt.title("GC")
    plt.legend()
    plt.ylabel('V_iono [m/s]')
    #axE = plt.twinx()
    #axE.plot(rMin_cc[:,0], E[:,0], c='green')
    #axE.set_ylabel('Energy [kev]')


def plotStep(raijuInfo, n, k):
    
    with h5.File(raijuInfo.fname, 'r') as f5:
        alamc = f5['alamc'][k]
        s5 = f5[raijuInfo.stepStrs[n]]
        xmin = ru.getVar(s5, 'xmin')
        ymin = ru.getVar(s5, 'ymin')
        etas = ru.getVar(s5, 'eta')
        bVol = ru.getVar(s5, 'bVol')


    fig = plt.figure()
    ax = fig.add_subplot(111)
    a = ax.pcolormesh(xmin, ymin, etas[:,:,k])
    plt.colorbar(a)

    """
    energy = alamc*bVol**(-2/3)*1e-3  # [keV]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    a = ax.pcolormesh(xmin, ymin, energy)
    plt.colorbar(a)
    """
    plt.show()


if __name__=="__main__":
    # Open file
    f5name = os.path.join(dir, raih5fname)
    f5 = h5.File(f5name, 'r')

    # See if we have debug turned on
    for k in f5.attrs.keys():
        print(k,": ",f5.attrs[k])
    doDebug = f5.attrs['doDebugOutput']
    if not doDebug:
        print("ERROR: We assume we have all debug info. put KAIJU/RAIJU/debug/debugOutput=T in xml file")
        quit()

    # General file info
    raijuInfo = ru.RAIJUInfo.getInfo(f5name)
    for k in vars(raijuInfo).keys():
        attr = getattr(raijuInfo, k)
        if isinstance(attr, list):
            print(k,": list, len",len(attr))
        elif isinstance(attr, np.ndarray):
            print(k,": list, shape", attr.shape)
        else:
            print(k,": ",attr)

    # Get grid info
    xIono = ru.getVar(f5, 'X')
    yIono = ru.getVar(f5, 'Y')
    Ni, Nj = xIono.shape
    Ni -= 1; Nj -= 1
    print("Num cells Ni, Nj:", Ni, ",", Nj)
    spacing_th = xIono[1:,0] - xIono[:-1,0]
    spacing_ph = yIono[0,1:] - yIono[0,:-1]
    tol_spacing = 1e-4  # [deg]

    # Planet params
    ri_m = f5['Planet'].attrs['Rad_ionosphere']
    rp_m = f5['Planet'].attrs['Rad_surface']
    Bs_nt = f5['Planet'].attrs['Mag Moment'] * kd.G2nT

    #checkCorot(raijuInfo, 10)
    
    spc_p = raijuInfo.species[ru.flavs_s['HOTP']]
    spc_k = spc_p.kStart + int(spc_p.N//2)
    #checkGC(raijuInfo, 10, spc_k)


    #checkBVol(raijuInfo)
    #plt.show()

    plotStep(raijuInfo, 20, 150)