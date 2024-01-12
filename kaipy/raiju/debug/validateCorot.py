import os
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import kaipy.kdefs as kd
import kaipy.kaiTools as kt
import kaipy.kaiH5 as kh5
import kaipy.kaiViz as kv
import kaipy.raiju.raijuutils as ru
import kaipy.raiju.debug.testAdvection as rta

def vCorot_rad(Rp_m: float, cPot: float, b_surf: float) -> float:
    """
    Calculate corotation velocity [m/s] in the ionosphere
    Rp: Planetary radius [m]
    cPot: Corotation potential [kV]
    b_surf: surface B-field [nT]
    r: eval equatorial radius [Rp]

    vCorot [rad/s] = cPot * b_surf / Rp^2
    vCorot [m/s] = cPot * b_surf * r [Rp] / Rp^2 = cPos * b_surf * r / Rp
    """
    v = (cPot*1e3) / (b_surf * 1e-9) /  Rp_m**2  # [rad/s]
    return v


def plotSlice(ax: mpl.axes.Axes, raiInfo: ru.RAIJUInfo, n, i=-1, doFollow=False):
    """
        Plot eta along i-slice

    """
    idxPsph = raiInfo.species[ru.flavs_s['PSPH']].kStart
    with h5.File(raiInfo.fname, 'r') as f5:

        y = ru.getVar(f5, 'Y')[0,:]
        ycc = kt.to_center1D(y)
        ycc0 = kt.to_center1D(y)  # [give eta0 its own ycc so we can manipulate it with doFollow

        eta0 = ru.getVar(f5[raiInfo.stepStrs[0]], 'eta')[:,:,idxPsph]
 
        s5 = f5[raiInfo.stepStrs[n]]
        eta = ru.getVar(s5, 'eta')[:,:,idxPsph]
        dt = ru.getVar(s5, 'dtk')[idxPsph]
        nSteps = ru.getVar(s5, 'nStepk')[idxPsph]

        if doFollow:
            # Shift ycc0 so that eta0 follows analytic corotation velocity 
            rp_m  = f5['Planet'].attrs['Rad_surface'   ]  # [m]
            Bs_nt = f5['Planet'].attrs['Mag Moment'    ] * kd.G2nT  # [nT]
            cPot  = f5['Planet'].attrs['Psi Corot']  # [kV]
            v_corot = vCorot_rad(rp_m, cPot, Bs_nt)  # [rad/s]

            t_step = s5.attrs['time']
            phi_shift = v_corot * t_step
            ycc0 += phi_shift
            ycc0[ycc0 > 2*np.pi] -= 2*np.pi
            # Determine where ycc0 wraps back to zero so we can draw as two separate parts and not get weird line along zero
            for j in range(len(ycc0)-1):
                if ycc0[j+1] < ycc0[j]:
                    jRoll = len(ycc0) - j - 1
                    ycc0 = np.roll(ycc0, jRoll)
                    eta0 = np.roll(eta0, jRoll)
                    break



    if i == -1:
        i = int(eta.shape[0]//2)

    ax.plot(ycc0, eta0[i,:], linestyle='--', color='k', label='nStp=0')
    ax.plot(ycc , eta [i,:], alpha=0.7, label='nStp={}'.format(nSteps))
    ax.scatter(ycc, eta[i,:], s=10, alpha=0.7)
    ax.legend()

    ax.set_xlabel('Phi [rad]')
    ax.set_ylabel('eta')
    ax.set_title("timestep = {:1.2f}".format(dt))


def plotPsphEq(raiInfo: ru.RAIJUInfo, xmin: np.ndarray, ymin: np.ndarray, n: int):
    pass


def makeVid(raiInfo: ru.RAIJUInfo, stride = 1):
    pass


if __name__=="__main__":

    # Run this from cirectory with config and scratch
    outDir = "tests"
    kh5.CheckDirOrMake(outDir)
    

    #fname = "raijuSA.raiju.h5"
    fname = "raijuSA_1d.raiju.h5"
    raiInfo = ru.RAIJUInfo.getInfo(fname)
    print(vars(raiInfo))

    rta.checkBVol(raiInfo)
    fname = os.path.join(outDir, "bVolTest.png")
    kv.savePic(fname)

    rta.checkCorot(raiInfo, int(raiInfo.steps[-1]))
    fname = os.path.join(outDir, "vCorot.png")
    kv.savePic(fname)

    # rta.makeVid(raiInfo, outDir, stride=10)

    fig = plt.figure()
    ax = fig.add_subplot()

    plotSlice(ax, raiInfo, raiInfo.steps[1], doFollow=True)
    kv.savePic('testfig_s1.png')
    ax.clear()

    plotSlice(ax, raiInfo, raiInfo.steps[10], doFollow=True)
    kv.savePic('testfig_s10.png')
    ax.clear()

    plotSlice(ax, raiInfo, raiInfo.steps[20], doFollow=True)
    kv.savePic('testfig_s20.png')
    ax.clear()

    plotSlice(ax, raiInfo, raiInfo.steps[-1], doFollow=True)
    kv.savePic('testfig_sEnd.png')
    ax.clear()


