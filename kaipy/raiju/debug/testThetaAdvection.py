"""
Note: this is very coupled to assumptions within raijuIChelpers.F90:initRaijuIC_testThetaAdvection
"""
import os
import h5py as h5
import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import kaipy.kdefs as kd
import kaipy.kaiTools as kt
import kaipy.kaiViz as kv

import kaipy.raiju.raijuutils as ru
import kaipy.raiju.debug.debughelper as dbg


def checkDerivLT(raiI: ru.RAIJUInfo, s5: h5.Group) -> bool:
    """
    dL/dTheta = 2 * sin^-3(Th)*cos(Th)
    """
    xmin = np.abs(ru.getVar(s5, 'xmin')[:,0])
    theta = ru.getVar(s5.file, 'X')[:,0]

    xcc  = kt.to_center1D(xmin )
    thcc = kt.to_center1D(theta)
    dl = xcc[1:] - xcc[:-1]
    dth = thcc[1:] - thcc[:-1]
    dldth =  np.abs(dl / dth)
    rhs = 2*np.sin(theta)**(-3)*np.cos(theta) * raiI.planet['Rad_ionosphere']/raiI.planet['Rad_surface']
    plt.plot(theta[1:-1], dldth, label='dldth')
    plt.plot(theta, rhs, label='rhs')
    plt.legend()
    plt.show()
    
    print(np.abs(dldth - rhs[1:-1])/rhs[1:-1])
    print(np.sum(np.abs(dldth - rhs[1:-1])))

    err = np.abs(dldth - rhs[1:-1])/rhs[1:-1]

    return np.all(err < 0.005)  # 0.5% error


def getEqEField(raiI: ru.RAIJUInfo, s5: h5.Group, doCheck=False) -> np.ndarray:
    """
    Calculates the electric field in the equatorial plane at all faces
    """

    # Corner values
    xEq = ru.getVar(s5, 'xmin')  # (Ni+1, Nj+1)
    yEq = ru.getVar(s5, 'ymin')  # (Ni+1, Nj+1)
    espot = ru.getVar(s5, 'espot')  # (Ni+1, Nj+1) [kV]
    Ni, Nj = xEq.shape
    Ni -= 1; Nj -= 1

    lenFaceEq = np.zeros((Ni+1,Nj+1,2))
    potDrop   = np.zeros((Ni+1,Nj+1,2))
    eFieldEq  = np.zeros((Ni+1,Nj+1,2))

    # Assuming cartesian
    lenFaceEq[:  ,:-1,0] = np.sqrt( (xEq[ :,1:] - xEq[:  ,:-1])**2 + (yEq[: ,1:] - yEq[:  ,:-1])**2 )
    lenFaceEq[:-1,:  ,1] = np.sqrt( (xEq[1:, :] - xEq[:-1,:  ])**2 + (yEq[1:, :] - yEq[:-1,:  ])**2 )

    potDrop[:  ,:-1,0] = espot[ :,1:] - espot[:  ,:-1]
    potDrop[:-1,:  ,1] = espot[1:, :] - espot[:-1,:  ]

    eFieldEq = (potDrop * 1e3) / (lenFaceEq * raiI.planet['Rad_surface'])  # [V/m]

    if doCheck:
        # Do it again in the ionosphere the exact way its done on the Fortran side
        f5 = s5.file
        lenFaceIono = ru.getVar(f5, 'lenFace')
        eFieldIono  = ru.getVar(s5, 'gradPotE')  # V/m
        j = int(Nj//2)
        
        print("eq:")
        print(eFieldEq[:,j,0]*lenFaceEq[:,j,0])
        print("iono:")
        print(eFieldIono[:,j,0]*lenFaceIono[:,j,0])

        plt.figure()
        plt.plot(eFieldEq[:,j,0]*lenFaceEq[:,j,0], label='eq')
        plt.plot(eFieldIono[:,j,0]*lenFaceIono[:,j,0], label='iono')
        plt.legend()
        plt.show()
                   
    return np.ma.masked_invalid(eFieldEq)  # [V/m]

    """
    # Test lenFaceEq
    i = 24; j = int(Nj//2)
    norm = kv.genNorm(0,np.max(lenFaceEq))
    plt.scatter(xEq[i:i+2,j:j+2], yEq[i:i+2,j:j+2])
    linemap = dbg.drawFaces(plt.gca(), xEq[i:i+2,j:j+2], yEq[i:i+2,j:j+2], lenFaceEq[i:i+2,j:j+2,:], norm=norm)
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes('right', size='5%', pad=0.05)
    plt.gcf().colorbar(linemap, cax=cax, orientation='vertical')
    plt.show()

    d_l_th = np.sqrt((xEq[i,j]-xEq[i,j+1])**2 + (yEq[i,j]-yEq[i,j+1])**2)
    d_u_th = np.sqrt((xEq[i+1,j]-xEq[i+1,j+1])**2 + (yEq[i+1,j]-yEq[i+1,j+1])**2)
    d_l_ph = np.sqrt((xEq[i,j]-xEq[i+1,j])**2 + (yEq[i,j]-yEq[i+1,j])**2)
    d_u_ph = np.sqrt((xEq[i,j+1]-xEq[i+1,j+1])**2 + (yEq[i,j+1]-yEq[i+1,j+1])**2)

    print(d_l_th,":",lenFaceEq[i,j,0])
    print(d_u_th,":",lenFaceEq[i+1,j,0])
    print(d_l_ph,":",lenFaceEq[i,j,1])
    print(d_u_ph,":",lenFaceEq[i,j+1,1])
    """

def getEqVelocity(s5: h5.Group, eFieldEq: np.ndarray, doCheck=False) -> np.ndarray:
    """
    Given equatorial e Field, calculate velocity at cell edges/faces
    """

    bmin = ru.getVar(s5, 'bminZ')  # [nT]
    Nic, Njc = bmin.shape  # Number of corners
    # Estimate bmin at face centers
    bminFaces = np.zeros((Nic,Njc,2))

    bminFaces[:  ,:-1,0] = 0.5*(bmin[ :,1:] + bmin[:  ,:-1])
    bminFaces[:-1,:  ,1] = 0.5*(bmin[1:, :] + bmin[:-1,:  ])


    velEq = eFieldEq / (bminFaces * 1e-9)  # V/m / T = V/m * (m^2/V/s) = m/s

    if doCheck:
        print("Velocity check")
        j = int((Njc-1)//2)
        xEq = ru.getVar(s5, 'xmin')[:,j]  # [Rp]
        thIono = ru.getVar(s5.file, 'X')[:,j]  # 1D Thetas
        kPsph = raiI.species[ru.flavs_s['PSPH']].kStart
        iVelIono = ru.getVar(s5,'iVel')[:,:,kPsph,:]  # [m/s]

        # Check vEq
        jStart, jEnd = getGradBnds(s5)
        phi = ru.getVar(s5.file, 'Y')[0,:]
        espot = ru.getVar(s5, 'espot')*1e3  # [V]
        delPhi = phi[jEnd]-phi[jStart]
        delPot = espot[0,jEnd] - espot[0,jStart]
        vEq = delPot/delPhi * xEq**2 / (raiI.planet['Mag Moment']*kd.G2nT*1e-9 * raiI.planet['Rad_surface'])

        plt.figure()
        x = np.arange(dtIono.shape[0])
        plt.plot(x, vEq, label='Partial-analytic')
        plt.scatter(x, velEq[:,j,0], label='Model')
        plt.legend()


        # Check vIono
        jStart, jEnd = getGradBnds(s5)
        phi = ru.getVar(s5.file, 'Y')[0,:]
        espot = ru.getVar(s5, 'espot')*1e3  # [V]
        delPhi = phi[jEnd]-phi[jStart]
        delPot = espot[0,jEnd] - espot[0,jStart]
        Bi = raiI.planet['Mag Moment']*kd.G2nT*1e-9 * (raiI.planet['Rad_ionosphere']/raiI.planet['Rad_surface'])**(-3)
        vI = delPot/delPhi / (2*Bi*np.sin(thIono)*np.cos(thIono)*raiI.planet['Rad_ionosphere'])

        plt.figure()
        plt.plot   (np.arange(0,Nic), vI, label='Partial-Analytic')
        plt.scatter(np.arange(0,Nic), iVelIono[:,j,0], label='model')
        plt.legend()
        plt.title("Velocity comparison")
        plt.xlabel("i")
        plt.ylabel("Velocity [m/s]")


        # Compare dt's from one cell center to the next
        xcc = kt.to_center1D(xEq)
        xDiffEq = raiI.planet['Rad_surface'] * (xcc[1:] - xcc[:-1])

        thIonoCC = kt.to_center1D(thIono)
        xDiffIono = raiI.planet['Rad_ionosphere'] * (thIonoCC[1:] - thIonoCC[:-1])  # [m]

        dtEq   = xDiffEq   / velEq   [1:-1,j,0]
        dtIono = xDiffIono / iVelIono[1:-1,j,0]

        print("dt comp")
        print('eq:')
        print(dtEq)
        print('iono')
        print(dtIono)
        print(dtEq/dtIono)

        plt.figure()
        x = np.arange(dtIono.shape[0])
        plt.plot(x, dtEq, label='eq')
        plt.scatter(x, dtIono, label='iono')
        plt.legend()
        plt.title("Dt in ionosphere vs equator")
        plt.xlabel("i")
        plt.ylabel("dt [s]")
        plt.show()
    

    return np.ma.masked_invalid(velEq)



def getGradBnds(s5: h5.Group):
    """
    Determine the start and end j bounds of the potential drop
    The exact j locations are the bounding values of the potential drop, 
     and all j's in-between has an intermediate potential value
    """

    espot = ru.getVar(s5, 'espot')
    jStart = -1
    jEnd = -1
    for j in range(espot.shape[1]-1):
        if jStart < 0 and espot[0,j+1] != espot[0,j]:
            jStart = j
        if jEnd < 0 and jStart > 0 and espot[0,j+1] == espot[0,j]:
            jEnd = j
            break
    return jStart, jEnd


def showSimpleEq(raiI: ru.RAIJUInfo, s5: h5.Group, faces='Efield', fname=None):
    
    xEq = ru.getVar(s5,'xmin')
    yEq = ru.getVar(s5,'ymin')
    kPsph = raiI.species[ru.flavs_s['PSPH']].kStart
    eta = ru.getVar(s5,'eta')[:,:,kPsph]

    norm_eta = kv.genNorm(0,1)
    cm_eta = 'viridis'

    fig = plt.figure(figsize=(9,6))
    ax = fig.add_subplot()

    map_eta = ax.pcolormesh(xEq, yEq, eta, norm=norm_eta, cmap=cm_eta)
    divider = make_axes_locatable(ax)
    cax_eta = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(map_eta, cax=cax_eta, orientation='vertical', ).set_label('Eta')


    if faces=='Efield':
        eFieldEq = getEqEField(raiI, s5)*1e3  # [mV/m]
        norm_eField = kv.genNorm(np.min(eFieldEq), np.max(eFieldEq))
        #cm_eField = 'RdPu'
        cm_eField = 'PuBuGn'

        jStart, jEnd = getGradBnds(s5)
        map_efd = dbg.drawFaces(ax, xEq[:,jStart:jEnd+1], yEq[:,jStart:jEnd+1], eFieldEq[:,jStart:jEnd+1], direction='theta', norm=norm_eField, cmap=cm_eField, iono=False)
        cax_efd = divider.append_axes('top', size='5%', pad=0.05)
        fig.colorbar(map_efd, cax=cax_efd, orientation='horizontal').set_label('E [mV/m]')
    elif faces=='Vel':
        eFieldEq = getEqEField(raiI, s5)  # [V/m]
        velEq = getEqVelocity(s5, eFieldEq)*1e-3  # [km/s]
        norm_vel = kv.genNorm(np.min(velEq), np.max(velEq))
        #cm_vel = 'RdPu'
        cm_vel = 'PuBuGn'

        jStart, jEnd = getGradBnds(s5)
        map_efd = dbg.drawFaces(ax, xEq[:,jStart:jEnd+1], yEq[:,jStart:jEnd+1], velEq[:,jStart:jEnd+1], direction='theta', norm=norm_vel, cmap=cm_vel, iono=False)
        cax_efd = divider.append_axes('top', size='5%', pad=0.05)
        fig.colorbar(map_efd, cax=cax_efd, orientation='horizontal').set_label('Velocity [km/s]')

    cax_efd.xaxis.set_label_position('top')
    cax_efd.xaxis.set_ticks_position('top')

    # Set bounds to only show interesting area
    espot = ru.getVar(s5,'espot')
    eMax = np.max(espot)
    eMin = np.min(espot)
    isInteresting = np.logical_and(espot > eMin, espot < eMax)
    #xBnds = [np.min(xEq[isInteresting]), np.max(xEq[isInteresting])]
    #yBnds = [-np.max(yEq[isInteresting]), np.max(yEq[isInteresting])]
    jMid = int(xEq.shape[1]//2)
    xBnds = [np.min(xEq), xEq[-1,jMid]]
    yBnds = [yEq[0,jStart], yEq[0,jEnd]]
    ax.set_xlim(xBnds)
    ax.set_ylim(yBnds)
    

def getCellTiming(raiI: ru.RAIJUInfo, s5: h5.Group, iStart=0) -> np.ndarray:
    """
    Main validator tool
    Calculate the time at which stuff should be arriving at each i cell
    Use starting eta position
    using face velocities and cell widths, figure out how long it should take eta to reach that cell
    """

    rp_m = raiI.planet['Rad_surface'] 

    eFieldEq = getEqEField(raiI, s5)  # [V/m]
    velEq = getEqVelocity(s5, eFieldEq)  # [m/s]


    Nic, Njc, _ = velEq.shape
    Ni = Nic-1; Nj = Njc-1
    j = int(Nj//2)

    # Can do 1D from now on
    velEq = velEq[:, j, 0]  # (Ni+1) Velocities along theta faces of our j column
    xmin = ru.getVar(s5, 'xmin')[:,j:j+2] * rp_m # [m]
    xmincc = kt.to_center2D(xmin)[:,0]
    deltaX = np.zeros(Ni)
    deltaX[1:] = xmincc[1:] - xmincc[:-1]

    tArrive = np.zeros(Ni)  # Time at which eta should arrive at each grid cell

    for i in range(iStart, Ni):
        t_this = deltaX[i] / velEq[i] # [s]
        tArrive[i] = tArrive[i-1] + t_this

    """
    print(velEq)
    plt.plot(xmincc/rp_m, tArrive)
    plt.show()
    """

    return tArrive

def getCellTiming_ana(raiI: ru.RAIJUInfo, f5: h5.Group, iStart=0) -> np.ndarray:
    """
    Cell timing, derived analytically
    integral of dx/v(x) from x0 to x1
    """

    s5 = f5[raiI.stepStrs[1]]
    espot = ru.getVar(s5, 'espot')
    phi = ru.getVar(f5, 'Y')

    B0 = raiI.planet['Mag Moment']*kd.G2nT*1e-9  # [T]
    delPot = (espot[0,-1] - espot[0,0])*1e3  # [V]
    jStart, jEnd = getGradBnds(s5)
    delPhi = phi[0,jEnd] - phi[0,jStart]
    print(B0)
    print(delPot)
    print(delPhi)

    Nic, Njc = espot.shape
    Ni = Nic-1; Nj = Njc-1
    j = int(Nj//2)
    rp_m = raiI.planet['Rad_surface'] 
    xmin = ru.getVar(s5, 'xmin')[:,j:j+2] * rp_m # [m]
    xmincc = kt.to_center2D(xmin)[:,0]
    Ni = xmincc.shape[0]
    x0 = xmincc[iStart]

    tArrive = np.zeros(Ni)

    for i in range(Ni):
        #tArrive[i] = -B0*delPhi/delPot* (1/xmincc[i] - 1/x0)  # [s]

        tArrive[i] = -B0*delPhi/delPot* (xmincc[i]**2/((xmincc[i]/rp_m)**3) - x0**2/((x0/rp_m)**3))  # [s]

    """
    plt.plot(tArrive, xmincc/rp_m )
    plt.show()
    quit()
    """

    return tArrive



def compareTiming(raiI: ru.RAIJUInfo, f5: h5.File, tArrive: np.ndarray):

    Ni,Nj = ru.getVar(f5,'BrCC').shape
    j = int(Nj//2)

    rp_m = raiI.planet['Rad_surface'] 
    s1 = f5[raiI.stepStrs[1]]
    xmin = ru.getVar(s1, 'xmin')[:,j:j+2] # [Rp]
    xmincc = kt.to_center2D(xmin)[:,0]
    
    etaBlock = np.zeros((raiI.Nt, Ni))

    kPsph = raiI.species[ru.flavs_s['PSPH']].kStart

    for n in range(raiI.Nt):
        s5 = f5[raiI.stepStrs[n]]
        eta = ru.getVar(s5, 'eta')[:,j,kPsph]
        etaBlock[n,:] = eta

    plt.figure(figsize=(12,8))

    norm_eta = kv.genNorm(0,1)
    plt.pcolormesh(raiI.times, xmincc, etaBlock.T, norm=norm_eta, cmap='viridis')
    plt.colorbar()
    plt.plot(tArrive, xmincc)
    plt.xlabel('Time [s]')
    plt.ylabel('X [R$_p$]')
    plt.show()


if __name__=="__main__":

    inDir = '.'
    oDir = 'test_ta'
    r5fname = 'raijuSA.raiju.h5'
    iStart = 0

    parser = argparse.ArgumentParser(description="Checks pure theta acvection, using the TTA IC in raijuICHelpers.F90")
    parser.add_argument('-indir',type=str,default=inDir,help="(default: %(default)s)")
    parser.add_argument('-odir',type=str,default=oDir ,help="(default: %(default)s)")
    parser.add_argument('-r5',type=str,default=r5fname,help="(default: %(default)s)")
    parser.add_argument('-iStart',type=int,default=iStart,help="(default: %(default)s)")
    args = parser.parse_args()

    indir = args.indir
    odir  = args.odir
    r5fname = args.r5
    iStart = args.iStart

    fullFname = os.path.join(indir, r5fname)
    raiI = ru.RAIJUInfo.getInfo(fullFname)

    # Attr info
    attrs = dbg.getAttrs(fullFname)
    for k in attrs.keys():
        print(k,":",attrs[k])

    f5 = h5.File(fullFname)
    s1 = f5[raiI.stepStrs[1]]

    # Ensure docorot=False (can just check gradPotCorot since it'll always be there)
    if np.sum(np.abs(ru.getVar(s1, 'gradPotCorot'))) != 0:
        print("ERROR: Can't have corotation on for this analysis")
        quit()

    
    showSimpleEq(raiI, f5[raiI.stepStrs[0]], faces='Efield')
    showSimpleEq(raiI, f5[raiI.stepStrs[0]], faces='Vel')
    plt.show()
    

    eFieldEq = getEqEField(raiI, s1, doCheck=False)
    getEqVelocity(s1, eFieldEq, doCheck=False)

    tArrive = getCellTiming(raiI, s1, iStart=iStart)
    tArrive_ana = getCellTiming_ana(raiI, f5, iStart=iStart)
    """
    plt.plot(tArrive,label='discrete')
    plt.plot(tArrive_ana,label='ana')
    plt.show()
    """
    compareTiming(raiI, f5, tArrive_ana)
    #compareTiming(raiI, f5, tArrive)