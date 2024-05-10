import os
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatch
import argparse
from tqdm import tqdm

import kaipy.kaiViz as kv
import kaipy.kaiTools as kt
import kaipy.kaiH5 as kh5
import kaipy.raiju.raijuutils as ru

import fileIO.fileIO as fIO


def makeCompFig(axs, X, Y, var1, var2, normMain, cmapMain, normDiff, cmapDiff):
    axBnd = 15
    rmin = X[-1,0]
    # Var 1
    ax = axs[0][0]
    ax.pcolormesh(X, Y, var1, norm=normMain, cmap=cmapMain)
    ax.set_xlabel('X [R$_p$]')
    ax.set_ylabel('Y [R$_p$]')
    ax.set_title('Raiju solve')
    kv.SetAx([-axBnd, axBnd, -axBnd, axBnd], ax=ax)
    kv.addEarth2D(ax=ax)
    circle_E = mpatch.Circle((0, 0), rmin, color='k', fill=False)
    ax.add_patch(circle_E)
    # Var 2
    ax = axs[1][0]
    ax.pcolormesh(X, Y, var2, norm=normMain, cmap=cmapMain)
    ax.set_xlabel('X [R$_p$]')
    ax.set_title('clawpak solve')
    kv.SetAx([-axBnd, axBnd, -axBnd, axBnd], ax=ax)
    kv.addEarth2D(ax=ax)
    circle_E = mpatch.Circle((0, 0), rmin, color='k', fill=False)
    ax.add_patch(circle_E)
    # Var 3
    diff = (var1 - var2)/2/(var1 + var2)
    ax = axs[2][0]
    #ax.pcolormesh(X, Y, var1 - var2, norm=normDiff, cmap=cmapDiff)
    ax.pcolormesh(X, Y, diff, norm=normDiff, cmap=cmapDiff)
    ax.set_xlabel('X [R$_p$]')
    kv.SetAx([-axBnd, axBnd, -axBnd, axBnd], ax=ax)
    kv.addEarth2D(ax=ax)
    circle_E = mpatch.Circle((0, 0), rmin, color='k', fill=False)
    ax.add_patch(circle_E)


def makeCompVid(raiI, rclI, varAccessStr, outdir, varName=None, normMain=None, cmapMain=None):
    """
    varAccessStr example: ru.getVar({}, 'Density')[:,:,0]
    """
    kh5.CheckDirOrMake(outdir)
    cbarsDone = False

    if varName is None:
        varName = varAccessStr.split("'")[1]
    if cmapMain is None:
        cmapMain = 'viridis'

    rai5 = h5.File(raiI.fname,'r')
    rcl5 = h5.File(rclI.fname,'r')

    Nt = min(raiI.Nt, rclI.Nt)

    # Set up figure
    fig = plt.figure(figsize=(32, 10))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1], wspace=0.1)
    axs = []
    gss = gs[:,0].subgridspec(1,2,width_ratios=[1,0.1], wspace=0.04)
    axs.append([fig.add_subplot(gss[0,0]), fig.add_subplot(gss[0,1])])
    gss = gs[:,1].subgridspec(1,2,width_ratios=[1,0.1], wspace=0.04)
    axs.append([fig.add_subplot(gss[0,0]), fig.add_subplot(gss[0,1])])
    gss = gs[:,2].subgridspec(1,2,width_ratios=[1,0.1], wspace=0.04)
    axs.append([fig.add_subplot(gss[0,0]), fig.add_subplot(gss[0,1])])
    cmapDiff = 'RdBu_r'

    # How many 0's do we need for filenames?
    n_pad = int(np.log10(Nt)) + 1
    n_pad = max(n_pad, 4)

    for i in tqdm(range(Nt)):
        s_rai = rai5[raiI.stepStrs[i]]
        s_rcl = rcl5[raiI.stepStrs[i]]
        xmin = ru.getVar(s_rai, 'xmin')
        ymin = ru.getVar(s_rai, 'ymin')
        time = s_rai.attrs['time']

        var1 = eval(varAccessStr.format('s_rai'))
        var2 = eval(varAccessStr.format('s_rcl'))

        if normMain is None:
            nMax = np.max((var1, var2))
            nMin = 1e-3*nMax
            #normMain = kv.genNorm(np.min((var1, var2)), np.max((var1, var2)))
            normMain = kv.genNorm(nMin, nMax, doLog=True)
        if not cbarsDone:
            #normDiff = kv.genNorm(1e-2*normMain.vmin, normMain.vmax, doLog=True)
            normDiff = kv.genNorm(-0.5, 0.5, doLog=False)
            kv.genCB(axs[0][1], normMain, cbT=varName, cM=cmapMain, doVert=True)
            kv.genCB(axs[1][1], normMain, cbT=varName, cM=cmapMain, doVert=True)
            kv.genCB(axs[2][1], normDiff, cbT='Diff', cM=cmapDiff, doVert=True)
            cbarsDone = True

        axs[0][0].clear()
        axs[1][0].clear()
        axs[2][0].clear()
        makeCompFig(axs, xmin, ymin, var1, var2, normMain, cmapMain, normDiff, cmapDiff)
        fig.suptitle(time)

        filename = "{}.{:0>{n}d}.png".format("vid", i, n=n_pad)
        outPath = os.path.join(outdir, filename)
        kv.savePic(outPath)


if __name__=="__main__":

    rundir = '/glade/derecho/scratch/sciola/raijudev/RaiRcmComp/SA/bell_1deg'
    outdir = 'comp'
    raiFname = 'raijuSA.raiju.h5'
    rclFname = 'raijuSA_claw.raiju.h5'

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d',default=rundir,
        help="Directory contianing two output files (default: %(default)s)"
    )
    parser.add_argument(
        '-out',default=outdir,
        help="Directory to store output plots (default: rundir+%(default)s)"
    )
    parser.add_argument(
        '-rai',default=raiFname,
        help="raiju solver filename (default: %(default)s)"
    )
    parser.add_argument(
        '-rcl',default=rclFname,
        help="raiju w/ clawpak solver filename (default: %(default)s)"
    )
    args = parser.parse_args()
    rundir = args.d
    outdir = os.path.join(rundir, args.out)
    raiFname = args.rai
    rclFname = args.rcl

    print("Getting raiI")
    raiI = kh5.H5Info.getInfo(raiFname)
    print("Getting rclI")
    rclI = kh5.H5Info.getInfo(rclFname)

    # Density first
    print("Making Density video")
    varAccessStr = "ru.getVar({}, 'Density')[:,:,0]"
    vidOut = os.path.join(outdir, "vid_den0")
    #norm = kv.genNorm(0.01, 5, doLog=True)
    norm = kv.genNorm(1e-9, 1e-6, doLog=True)
    makeCompVid(raiI, rclI, varAccessStr, vidOut, 'Total Density [#/cc]', normMain=norm, cmapMain='viridis')

    #print("Making eta video")
    #varAccessStr = "ru.getVar({}, 'eta')[:,:,27]"
    #vidOut = os.path.join(outdir, "vid_eta27")
    #norm = kv.genNorm(1e-15, 1e6, doLog=True)
    #makeCompVid(raiI, rclI, varAccessStr, vidOut, 'eta [#/cc*(Rp/nT)]', normMain=norm, cmapMain='viridis')

