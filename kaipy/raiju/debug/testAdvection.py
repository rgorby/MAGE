import os
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import cmasher as cm

import kaipy.kaiH5 as kh5
import kaipy.kdefs as kd
import kaipy.kaiTools as kt
import kaipy.kaiViz as kv
import kaipy.raiju.raijuutils as ru




# random values
rad2deg = 180/np.pi

# Settings
#dir = '/Users/sciolam1/Workspace/runs/local/raijudev/fluxOverhaul'
dir = '.'
raih5fname = 'raijuSA_o2.raiju.h5'
#raih5fname = 'raijuOWD_g.raiju.h5'

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

    vCorot [rad/s] = cPot / b_surf / Rp^2
    vCorot [m/s] = cPot / b_surf * r [Rp] / Rp^2 = cPos / b_surf * r / Rp
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
        rMin = ru.getVar(s5, 'xmin')[:,0]

    bVol = bVol[:,0]
    gradVM = gradVM[:,0,1]
    #rMin_cc = kt.to_center1D(rMin[:,0])
    sinColat = np.sqrt(1/rMin)
    colat = np.arcsin(sinColat)

    bVol_ana = calc_bVol_ana(colat, Bs_nt)
    grad_bVol_ana = calc_grad_bVol_ana(colat, Bs_nt, ri_m)

    gradVM_ana = (-2/3*bVol**(-5/3)) * grad_bVol_ana

    plt.figure()
    plt.plot   (rMin, bVol_ana , label='ana'  )
    plt.scatter(rMin, bVol, label='model', c='orange')
    plt.title('bVol')
    plt.yscale('log')
    plt.legend()

    plt.figure()
    plt.plot   (rMin, gradVM_ana , label='ana'  )
    plt.scatter(rMin, gradVM, label='model', c='orange')
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
        espot = ru.getVar(s5, 'espot')  # [kV]

        if np.any(espot > 0):
            print("Error: espot must be zero for checkCorot to work")
            return
        
        rMin    = ru.getVar(s5, 'xmin')[:,0]  # Just take all i's at j=0 (y=0)
        iVel_ph = ru.getVar(s5, 'iVel')[:,:,:,1]
        
    #rMin_cc = kt.to_center1D(rMin)
    vCorot_ana = vCorot(rp_m, ri_m, cPot, Bs_nt, rMin)
    # Just get velocity from plasmasphere channel
    kPsph = raiInfo.species[ru.flavs_s['PSPH']].kStart
    vCorot_model = iVel_ph[:,0,kPsph]

    print("corot:")
    print(vCorot_model/vCorot_ana)
    plt.figure()
    plt.plot(rMin, vCorot_ana, label='analytic')
    plt.scatter(rMin, vCorot_model, s=10, c='orange', label='model')
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
            return
        
        rMin = ru.getVar(s5, 'xmin')
        iVel_ph = ru.getVar(s5, 'iVel')[:,:,:,1]

    # Find species info
    iSpc = 0
    while raiInfo.species[iSpc].kEnd < k: iSpc += 1
    spc = raiInfo.species[iSpc]
    print("Species: ",spc.name)
    signQ = np.sign(spc.q)
    #rMin_cc = kt.to_center2D(rMin)
    E = alamc[k]*bVol**(-2/3)*1e-3  # [keV]
    vGC_eq_ana = vGradEq(rp_m, ri_m, signQ, Bs_nt, E, rMin)
    vGC_iso_ana = vGCiso(rp_m, ri_m, Bs_nt, alamc[k], rMin)
    vCorot_ana  = vCorot(rp_m, ri_m, cPot , Bs_nt   , rMin)

    vPh_model = iVel_ph[:,:,k]
    vGC_model = vPh_model - vCorot_ana

    plt.figure()
    plt.plot(   rMin[:,0], vGC_eq_ana   [:,0], label='ana (grad_eq)')
    plt.plot(   rMin[:,0], vGC_iso_ana  [:,0], c='green', label='ana (GC_iso)')
    plt.scatter(rMin[:,0], vGC_model[:,0], s=10, c='orange', label='model')
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


def makeVid(raijuInfo, outdir, stride=1):

    fontsize_title=6
    vidOut = 'raijuVid_o2'
    kh5.CheckDirOrMake(os.path.join(outdir, vidOut))
    
    fig = plt.figure()

    gs = gridspec.GridSpec(3, 4, width_ratios = [1, 1, 1, 0.2], height_ratios=[1, 1.0, 1.0], hspace=0.2)

    Ax_eta_psph = fig.add_subplot(gs[0, 0])
    
    Ax_eta_pLow  = fig.add_subplot(gs[0, 1])
    Ax_eta_pMid  = fig.add_subplot(gs[1, 1])
    Ax_eta_pHigh = fig.add_subplot(gs[2, 1])

    Ax_eta_eLow  = fig.add_subplot(gs[0, 2])
    Ax_eta_eMid  = fig.add_subplot(gs[1, 2])
    Ax_eta_eHigh = fig.add_subplot(gs[2, 2])

    # Create the Colorbar Axes.
    AxC1 = fig.add_subplot(gs[0, 3])
    
    # Create the Colorbars``.
    #etaMax = 0.3*np.max(eta)
    norm_eta = kv.genNorm(1e-4,1e8, doLog=True)
    cm_eta = 'viridis'
    kv.genCB(AxC1, norm_eta, "eta [eV * (Rp/nT)^(2/3)]", cM=cm_eta, doVert=True)

    xBnds = [-15, 10]

    # Open file, get some values
    f5 = h5.File(raijuInfo.fname, 'r')
    alamc = f5['alamc'][:]
    spc_psph = raijuInfo.species[ru.flavs_s['PSPH']]
    k_psph = spc_psph.kStart
    spc_hotp = raijuInfo.species[ru.flavs_s['HOTP']]
    k_pLow  = spc_hotp.kStart
    k_pMid  = spc_hotp.kStart + int(spc_hotp.N//2)
    k_pHigh = spc_hotp.kEnd - 1
    spc_hote = raijuInfo.species[ru.flavs_s['HOTE']]
    k_eLow  = spc_hote.kStart
    k_eMid  = spc_hote.kStart + int(spc_hote.N//2)
    k_eHigh = spc_hote.kEnd - 1

    n_pad = int(np.log10(raijuInfo.Nt)) + 1
    nplt = 0
    for n in raijuInfo.steps[::stride]:
        filename = "vid.{:0>{npad}d}.png".format(nplt, npad=n_pad)
        filename = os.path.join(outdir, vidOut, filename)
        nplt += 1
        if os.path.exists(filename):
            continue
        if (np.mod(nplt,5) == 0):
            print("\tvid %d"%(n))

        s5 = f5[raijuInfo.stepStrs[n]]
        time = s5.attrs['time']
        xmin = ru.getVar(s5, 'xmin')
        ymin = ru.getVar(s5, 'ymin')
        eta  = ru.getVar(s5, 'eta' )
        bVol = ru.getVar(s5, 'bVol')
        active = ru.getVar(s5, 'active')
        topo = ru.getVar(s5, 'topo')

        active3d = np.broadcast_to(active[:,:,np.newaxis], eta.shape)
        eta = np.ma.masked_where(active3d != ru.domain["ACTIVE"], eta)
        bVol = np.ma.masked_where(topo != ru.topo['CLOSED'], bVol)

        energies = alamc[np.newaxis,np.newaxis,:] * bVol[:,:,np.newaxis]**(-2/3) * 1e-3  # [keV]
        energies = np.abs(energies)

        # Psph
        ru.plotXYMin(Ax_eta_psph , xmin, ymin, eta[:,:,k_psph ], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        Ax_eta_psph.set_title('Psph', fontsize=fontsize_title)
        kv.setBndsByAspect(Ax_eta_psph, xBnds)

        ru.plotXYMin(Ax_eta_pLow , xmin, ymin, eta[:,:,k_pLow ], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        ru.plotXYMin(Ax_eta_pMid , xmin, ymin, eta[:,:,k_pMid ], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        ru.plotXYMin(Ax_eta_pHigh, xmin, ymin, eta[:,:,k_pHigh], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        Ax_eta_pLow.xaxis.set_major_formatter(plt.NullFormatter())
        Ax_eta_pLow.yaxis.set_major_formatter(plt.NullFormatter())
        Ax_eta_pMid.xaxis.set_major_formatter(plt.NullFormatter())
        Ax_eta_pLow .set_title('H+ $\lambda$={}, E={:0.2f}-{:0.2f} keV'.format(k_pLow , np.min(energies[:,:,k_pLow ]), np.max(energies[:,:,k_pLow ])), fontsize=fontsize_title)
        Ax_eta_pMid .set_title('H+ $\lambda$={}, E={:0.2f}-{:0.2f} keV'.format(k_pMid , np.min(energies[:,:,k_pMid ]), np.max(energies[:,:,k_pMid ])), fontsize=fontsize_title)
        Ax_eta_pHigh.set_title('H+ $\lambda$={}, E={:0.2f}-{:0.2f} keV'.format(k_pHigh, np.min(energies[:,:,k_pHigh]), np.max(energies[:,:,k_pHigh])), fontsize=fontsize_title)
        kv.setBndsByAspect(Ax_eta_pLow, xBnds)
        kv.setBndsByAspect(Ax_eta_pMid, xBnds)
        kv.setBndsByAspect(Ax_eta_pHigh, xBnds)

        ru.plotXYMin(Ax_eta_eLow , xmin, ymin, eta[:,:,k_eLow ], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        ru.plotXYMin(Ax_eta_eMid , xmin, ymin, eta[:,:,k_eMid ], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        ru.plotXYMin(Ax_eta_eHigh, xmin, ymin, eta[:,:,k_eHigh], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        Ax_eta_eLow .xaxis.set_major_formatter(plt.NullFormatter())
        Ax_eta_eLow .yaxis.set_major_formatter(plt.NullFormatter())
        Ax_eta_eMid .xaxis.set_major_formatter(plt.NullFormatter())
        Ax_eta_eMid .yaxis.set_major_formatter(plt.NullFormatter())
        Ax_eta_eHigh.yaxis.set_major_formatter(plt.NullFormatter())
        Ax_eta_eLow .set_title('e- $\lambda$={}, E={:0.2f}-{:0.2f} keV'.format(k_eLow , np.min(energies[:,:,k_eLow ]), np.max(energies[:,:,k_eLow ])), fontsize=fontsize_title)
        Ax_eta_eMid .set_title('e- $\lambda$={}, E={:0.2f}-{:0.2f} keV'.format(k_eMid , np.min(energies[:,:,k_eMid ]), np.max(energies[:,:,k_eMid ])), fontsize=fontsize_title)
        Ax_eta_eHigh.set_title('e- $\lambda$={}, E={:0.2f}-{:0.2f} keV'.format(k_eHigh, np.min(energies[:,:,k_eHigh]), np.max(energies[:,:,k_eHigh])), fontsize=fontsize_title)
        kv.setBndsByAspect(Ax_eta_eLow, xBnds)
        kv.setBndsByAspect(Ax_eta_eMid, xBnds)
        kv.setBndsByAspect(Ax_eta_eHigh, xBnds)

        fig.suptitle("step# {}; t = {:0.2f} s".format(n, time))
        kv.savePic(filename)

        Ax_eta_psph.clear()
        Ax_eta_pLow .clear()
        Ax_eta_pMid .clear()
        Ax_eta_pHigh.clear()
        Ax_eta_eLow .clear()
        Ax_eta_eMid .clear()
        Ax_eta_eHigh.clear()



def makeVid2(raijuInfo, outdir, stride=1):

    fontsize_title=6
    vidOut = 'raijuVid2'
    kh5.CheckDirOrMake(os.path.join(outdir, vidOut))

    fig = plt.figure(figsize=(20,10))
    cmW = 0.2
    gs = gridspec.GridSpec(2, 8, width_ratios = [1, cmW, 1, cmW, 1, cmW, 1, cmW], height_ratios=[1, 1], hspace=0.2)

    # Top row: drivers
    #  MHD bVol | MHD den | MHD press | None

    # Botom row: RAIJU
    # Psph | Mid H+ | Mid e- | total Pressure

    norm_den = kv.genNorm(0.1, 1000, doLog=True)
    norm_press = kv.genNorm(0, 70, doLog=False)
    norm_bVol = kv.genNorm(3e-3, 0.3, doLog=True)
    norm_eta = kv.genNorm(1e-1,1e10, doLog=True)
    cm_den = 'viridis'
    cm_press = 'plasma'
    cm_bVol = cm.prinsenvlag
    cm_eta = cm.voltage

    axBvol  = fig.add_subplot(gs[0,0])
    axDen   = fig.add_subplot(gs[0,2])    
    axPress = fig.add_subplot(gs[0,4])
    axCM_bvol  = fig.add_subplot(gs[0,1])
    axCM_den   = fig.add_subplot(gs[0,3])
    axCM_press = fig.add_subplot(gs[0,5])

    axPsph   = fig.add_subplot(gs[1,0])
    axMidH   = fig.add_subplot(gs[1,2])
    axMidE   = fig.add_subplot(gs[1,4])
    axRPress = fig.add_subplot(gs[1,6])
    #axCM_psph   = fig.add_subplot(gs[1,1])
    #axCM_midH   = fig.add_subplot(gs[3,1])
    #axCM_midE   = fig.add_subplot(gs[5,1])
    axCM_eta   = fig.add_subplot(gs[1,1])
    #axCM_RPress = fig.add_subplot(gs[7,1])

    kv.genCB(axCM_bvol , norm_bVol , 'Flux tube volume', cm_bVol, doVert=True)
    kv.genCB(axCM_den  , norm_den  , 'Density [#/cc]', cm_den, doVert=True)
    kv.genCB(axCM_press, norm_press, 'Pressure [nPa]', cm_press, doVert=True)
    
    kv.genCB(axCM_eta, norm_eta, 'eta', cm_eta, doVert=True)

    xBnds = [-15, 10]

    f5 = h5.File(raijuInfo.fname, 'r')
    alamc = f5['alamc'][:]
    spc_psph = raijuInfo.species[ru.flavs_s['PSPH']]
    k_psph = spc_psph.kStart
    spc_hotp = raijuInfo.species[ru.flavs_s['HOTP']]
    k_pMid  = spc_hotp.kStart + int(spc_hotp.N//2)
    spc_hote = raijuInfo.species[ru.flavs_s['HOTE']]
    k_eMid  = spc_hote.kStart + int(spc_hote.N//2)

    n_pad = int(np.log10(raijuInfo.Nt)) + 1
    nplt = 0
    for n in raijuInfo.steps[::stride]:
        filename = "vid.{:0>{npad}d}.png".format(nplt, npad=n_pad)
        filename = os.path.join(outdir, vidOut, filename)
        nplt += 1
        if os.path.exists(filename):
            continue
        if (np.mod(nplt,5) == 0):
            print("\tvid %d"%(n))

        s5 = f5[raijuInfo.stepStrs[n]]
        time = s5.attrs['time']
        xmin = ru.getVar(s5, 'xmin')
        ymin = ru.getVar(s5, 'ymin')
        active = ru.getVar(s5, 'active')
        topo = ru.getVar(s5, 'topo')
        eta  = ru.getVar(s5, 'eta' )
        bVol = ru.getVar(s5, 'bVol')
        denIn = ru.getVar(s5,'Davg_in')[:,:,2]
        pressIn = ru.getVar(s5,'Pavg_in')[:,:,2]
        pressOut = ru.getVar(s5,'Pressure')[:,:,2]

        bVol = kt.to_center2D(bVol)

        active3d = np.broadcast_to(active[:,:,np.newaxis], eta.shape)
        eta      = np.ma.masked_where(active3d != ru.domain["ACTIVE"], eta     )
        #bVol    = np.ma.masked_where(topo != ru.topo['CLOSED'], bVol   )
        bVol = np.ma.masked_where(active != ru.domain['ACTIVE'], bVol)
        denIn   = np.ma.masked_where(active != ru.domain['ACTIVE'], denIn  )
        pressIn = np.ma.masked_where(active != ru.domain['ACTIVE'], pressIn)
        pressOut = np.ma.masked_where(active != ru.domain["ACTIVE"], pressOut)

        energies = alamc[np.newaxis,np.newaxis,:] * bVol[:,:,np.newaxis]**(-2/3) * 1e-3  # [keV]
        energies = np.abs(energies)

        # energy eval loc
        evalLoc = [-4, 0]
        distSq = (xmin-evalLoc[0])**2 + (ymin-evalLoc[1])**2  # Each slice point's distance from particle
        aMinFlattened = distSq.argmin()  # Index of min distance, if it was a flattened array
        i, j = np.unravel_index(aMinFlattened, xmin.shape)  # Convert back into i,j location
        eEval_H = energies[i, j, k_pMid]
        eEval_E = energies[i, j, k_eMid]

        # MHD bVol
        ru.plotXYMin(axBvol , xmin, ymin, bVol, norm=norm_bVol, cmap=cm_bVol, shading='auto', lblsize=None)
        kv.setBndsByAspect(axBvol, xBnds)
        # MHD den
        ru.plotXYMin(axDen , xmin, ymin, denIn, norm=norm_den, cmap=cm_den, shading='auto', lblsize=None)
        kv.setBndsByAspect(axDen, xBnds)
        # MHD press
        ru.plotXYMin(axPress , xmin, ymin, pressIn, norm=norm_press, cmap=cm_press, shading='auto', lblsize=None)
        kv.setBndsByAspect(axPress, xBnds)

        # RAIJU psph
        ru.plotXYMin(axPsph , xmin, ymin, eta[:,:,k_psph ], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        axPsph.set_title('Plasmasphere', fontsize=fontsize_title)
        kv.setBndsByAspect(axPsph, xBnds)
        # RAIJU mid H+
        ru.plotXYMin(axMidH , xmin, ymin, eta[:,:,k_pMid ], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        axMidH.set_title('H+ $\lambda$={}, {:0.2f} keV @ {} R$_E$'.format(k_pMid, eEval_H, evalLoc[0]), fontsize=fontsize_title)
        kv.setBndsByAspect(axMidH, xBnds)
        # RAIJU mid e-
        ru.plotXYMin(axMidE , xmin, ymin, eta[:,:,k_eMid ], norm=norm_eta, cmap=cm_eta, shading='auto', lblsize=None)
        axMidE.set_title('e- $\lambda$={}, {:0.2f} keV @ {} R$_E$'.format(k_eMid, eEval_E, evalLoc[0]), fontsize=fontsize_title)
        kv.setBndsByAspect(axMidE, xBnds)
        # RAIJU total Pressure
        ru.plotXYMin(axRPress , xmin, ymin, pressOut, norm=norm_press, cmap=cm_press, shading='auto', lblsize=None)
        axRPress.set_title('Total pressure', fontsize=fontsize_title)
        kv.setBndsByAspect(axRPress, xBnds)

        axBvol.xaxis.set_major_formatter(plt.NullFormatter())
        axDen.xaxis.set_major_formatter(plt.NullFormatter())
        axPress.xaxis.set_major_formatter(plt.NullFormatter())
        axBvol.yaxis.set_major_formatter(plt.NullFormatter())
        axDen.yaxis.set_major_formatter(plt.NullFormatter())

        axPsph.yaxis.set_major_formatter(plt.NullFormatter())
        axMidH.yaxis.set_major_formatter(plt.NullFormatter())
        axMidE.yaxis.set_major_formatter(plt.NullFormatter())


        fig.suptitle("step# {}; t = {:0.2f} s".format(n, time))
        kv.savePic(filename)

        axBvol.clear()
        axDen.clear()
        axPress .clear()
        axPsph.clear()
        axMidH .clear()
        axMidE .clear()
        axRPress.clear()


def makeVid3(raijuInfo, outdir, stride=1):

    fontsize_title=18
    vidOut = 'raijuVid3'
    kh5.CheckDirOrMake(os.path.join(outdir, vidOut))

    fig = plt.figure(figsize=(24,10))
    cmW = 0.1
    gs = gridspec.GridSpec(4, 3, 
                           width_ratios = [1, 1, 1], 
                           height_ratios=[cmW, 1, 1, cmW], 
                           hspace=0.12, wspace=0.1)

    # Top row: drivers
    # MHD den | MHD press | MHD bVol

    # Botom row: RAIJU
    # Psph | Mid H+ | Mid e-

    norm_den = kv.genNorm(0.1, 1e4, doLog=True)
    norm_press = kv.genNorm(0, 50, doLog=False)
    norm_bVol = kv.genNorm(3e-3, 0.3, doLog=True)
    norm_psph = kv.genNorm(1e-1,1e10, doLog=True)
    norm_midh = kv.genNorm(5e2,1e8, doLog=True)
    norm_mide = kv.genNorm(5e2,1e8, doLog=True)
    cm_den = 'viridis'
    cm_press = 'plasma'
    cm_bVol = cm.prinsenvlag
    cm_eta = cm.voltage

    axCM_den   = fig.add_subplot(gs[0,0])
    axCM_press = fig.add_subplot(gs[0,1])
    axCM_bvol  = fig.add_subplot(gs[0,2])
    axDen   = fig.add_subplot(gs[1,0])    
    axPress = fig.add_subplot(gs[1,1])
    axBvol  = fig.add_subplot(gs[1,2])
    
    axPsph   = fig.add_subplot(gs[2,0])
    axMidH   = fig.add_subplot(gs[2,1])
    axMidE   = fig.add_subplot(gs[2,2])
    axCM_Psph   = fig.add_subplot(gs[3,0])
    axCM_MidH   = fig.add_subplot(gs[3,1])
    axCM_MidE   = fig.add_subplot(gs[3,2])

    kv.genCB(axCM_bvol , norm_bVol , 'Flux tube volume', cm_bVol , doVert=False, cbSz='large')
    kv.genCB(axCM_den  , norm_den  , 'Density [#/cc]'  , cm_den  , doVert=False, cbSz='large')
    kv.genCB(axCM_press, norm_press, 'Pressure [nPa]'  , cm_press, doVert=False, cbSz='large')
    axCM_bvol .xaxis.set_ticks_position('top')
    axCM_den  .xaxis.set_ticks_position('top')
    axCM_press.xaxis.set_ticks_position('top')
    axCM_bvol .xaxis.set_label_position('top')
    axCM_den  .xaxis.set_label_position('top')
    axCM_press.xaxis.set_label_position('top')

    
    kv.genCB(axCM_Psph, norm_psph, 'eta', cm_eta, doVert=False, cbSz='large')
    kv.genCB(axCM_MidH, norm_midh, 'eta', cm_eta, doVert=False, cbSz='large')
    kv.genCB(axCM_MidE, norm_mide, 'eta', cm_eta, doVert=False, cbSz='large')

    xBnds = [-15, 10]

    f5 = h5.File(raijuInfo.fname, 'r')
    alamc = f5['alamc'][:]
    spc_psph = raijuInfo.species[ru.flavs_s['PSPH']]
    k_psph = spc_psph.kStart
    spc_hotp = raijuInfo.species[ru.flavs_s['HOTP']]
    k_pMid  = spc_hotp.kStart + int(spc_hotp.N//4)
    spc_hote = raijuInfo.species[ru.flavs_s['HOTE']]
    k_eMid  = spc_hote.kStart + int(spc_hote.N//3)

    n_pad = int(np.log10(raijuInfo.Nt)) + 1
    nplt = 0
    for n in raijuInfo.steps[::stride]:
        filename = "vid.{:0>{npad}d}.png".format(nplt, npad=n_pad)
        filename = os.path.join(outdir, vidOut, filename)
        nplt += 1
        if os.path.exists(filename):
            continue
        if (np.mod(nplt,5) == 0):
            print("\tvid %d"%(n))

        s5 = f5[raijuInfo.stepStrs[n]]
        time = s5.attrs['time']
        xmin = ru.getVar(s5, 'xmin')[4:,:]
        ymin = ru.getVar(s5, 'ymin')[4:,:]
        active = ru.getVar(s5, 'active')[4:,:]
        topo = ru.getVar(s5, 'topo')[4:,:]
        eta  = ru.getVar(s5, 'eta' )[4:,:]
        bVol = ru.getVar(s5, 'bVol')[4:,:]
        denIn = ru.getVar(s5,'Davg_in')[4:,:,2]
        pressIn = ru.getVar(s5,'Pavg_in')[4:,:,2]
        pressOut = ru.getVar(s5,'Pressure')[4:,:,2]

        bVol = kt.to_center2D(bVol)

        active3d = np.broadcast_to(active[:,:,np.newaxis], eta.shape)
        eta      = np.ma.masked_where(active3d != ru.domain["ACTIVE"], eta     )
        #bVol    = np.ma.masked_where(topo != ru.topo['CLOSED'], bVol   )
        bVol = np.ma.masked_where(active != ru.domain['ACTIVE'], bVol)
        denIn   = np.ma.masked_where(active != ru.domain['ACTIVE'], denIn  )
        pressIn = np.ma.masked_where(active != ru.domain['ACTIVE'], pressIn)
        pressOut = np.ma.masked_where(active != ru.domain["ACTIVE"], pressOut)

        energies = alamc[np.newaxis,np.newaxis,:] * bVol[:,:,np.newaxis]**(-2/3) * 1e-3  # [keV]
        energies = np.abs(energies)

        # energy eval loc
        evalLoc = [-4, 0]
        distSq = (xmin-evalLoc[0])**2 + (ymin-evalLoc[1])**2  # Each slice point's distance from particle
        aMinFlattened = distSq.argmin()  # Index of min distance, if it was a flattened array
        i, j = np.unravel_index(aMinFlattened, xmin.shape)  # Convert back into i,j location
        eEval_H = energies[i, j, k_pMid]
        eEval_E = energies[i, j, k_eMid]

        # MHD den
        ru.plotXYMin(axDen , xmin, ymin, denIn, norm=norm_den, cmap=cm_den, shading='auto', lblsize=None)
        kv.setBndsByAspect(axDen, xBnds)
        kv.addEarth2D(ax=axDen)
        # MHD press
        ru.plotXYMin(axPress , xmin, ymin, pressIn, norm=norm_press, cmap=cm_press, shading='auto', lblsize=None)
        kv.setBndsByAspect(axPress, xBnds)
        kv.addEarth2D(ax=axPress)
        # MHD bVol
        ru.plotXYMin(axBvol , xmin, ymin, bVol, norm=norm_bVol, cmap=cm_bVol, shading='auto', lblsize=None)
        kv.setBndsByAspect(axBvol, xBnds)
        kv.addEarth2D(ax=axBvol)

        textLoc = [0.5, 0.1]
        textSize = 14
        # RAIJU psph
        ru.plotXYMin(axPsph , xmin, ymin, eta[:,:,k_psph ], norm=norm_psph, cmap=cm_eta, shading='auto', lblsize=None)
        t = axPsph.text(textLoc[0], textLoc[1],'Plasmasphere', \
                horizontalalignment='center', verticalalignment='top', \
                transform=axPsph.transAxes,fontsize=textSize)
        t.set_bbox(dict(facecolor='white',edgecolor='black',alpha=0.85))
        kv.setBndsByAspect(axPsph, xBnds)
        kv.addEarth2D(ax=axPsph)
        # RAIJU mid H+
        ru.plotXYMin(axMidH , xmin, ymin, eta[:,:,k_pMid ], norm=norm_midh, cmap=cm_eta, shading='auto', lblsize=None)
        t = axMidH.text(textLoc[0], textLoc[1],'H+, {:0.2f} keV @ {} R$_E$'.format(eEval_H, np.abs(evalLoc[0])), \
                horizontalalignment='center', verticalalignment='top', \
                transform=axMidH.transAxes,fontsize=textSize)
        t.set_bbox(dict(facecolor='white',edgecolor='black',alpha=0.85))
        kv.setBndsByAspect(axMidH, xBnds)
        kv.addEarth2D(ax=axMidH)
        # RAIJU mid e-
        ru.plotXYMin(axMidE , xmin, ymin, eta[:,:,k_eMid ], norm=norm_mide, cmap=cm_eta, shading='auto', lblsize=None)
        t = axMidE.text(textLoc[0], textLoc[1],'e-, {:0.2f} keV @ {} R$_E$'.format(eEval_E, np.abs(evalLoc[0])), \
                horizontalalignment='center', verticalalignment='top', \
                transform=axMidE.transAxes,fontsize=textSize)
        t.set_bbox(dict(facecolor='white',edgecolor='black',alpha=0.85))
        kv.setBndsByAspect(axMidE, xBnds)
        kv.addEarth2D(ax=axMidE)

        #axDen  .xaxis.set_major_formatter(plt.NullFormatter())
        #axPress.xaxis.set_major_formatter(plt.NullFormatter())
        #axBvol .xaxis.set_major_formatter(plt.NullFormatter())
        axDen  .tick_params(axis='x', which='major', pad=-0.9)
        axPress.tick_params(axis='x', which='major', pad=-0.9)
        axBvol .tick_params(axis='x', which='major', pad=-0.9)
        axDen  .yaxis.set_major_formatter(plt.NullFormatter())
        axPress.yaxis.set_major_formatter(plt.NullFormatter())

        axPsph.xaxis.set_major_formatter(plt.NullFormatter())
        axMidH.xaxis.set_major_formatter(plt.NullFormatter())
        axMidE.xaxis.set_major_formatter(plt.NullFormatter())
        axPsph.xaxis.set_ticks_position('top')
        axMidH.xaxis.set_ticks_position('top')
        axMidE.xaxis.set_ticks_position('top')
        axPsph.yaxis.set_major_formatter(plt.NullFormatter())
        axMidH.yaxis.set_major_formatter(plt.NullFormatter())


        fig.suptitle("step# {}; t = {:0.2f} s".format(n, time))
        kv.savePic(filename)

        axBvol.clear()
        axDen.clear()
        axPress .clear()
        axPsph.clear()
        axMidH .clear()
        axMidE .clear()



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

    # checkCorot(raijuInfo, 1)
    
    spc_p = raijuInfo.species[ru.flavs_s['HOTP']]
    spc_k = spc_p.kStart + int(spc_p.N//2)
    # checkGC(raijuInfo, 1, spc_k)


    #checkBVol(raijuInfo)
    # plt.show()

    # plotStep(raijuInfo, 0, 10)

    #makeVid(raijuInfo, dir)
    #makeVid2(raijuInfo, dir, 5)
    makeVid3(raijuInfo, dir, 1)
