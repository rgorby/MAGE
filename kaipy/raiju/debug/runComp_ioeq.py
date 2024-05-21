import os
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from tqdm import tqdm
import cmasher as cmr
import argparse

import kaipy.kaiViz as kv
import kaipy.kaiH5 as kh5
import kaipy.kaiTools as kt

import kaipy.raiju.raijuutils as ru



def getInfo(fname, model):
    if model=='raiju':
        return ru.RAIJUInfo.getInfo(fname)
    if model=='mhdrcm':
        return kh5.H5Info.getInfo(fname)
    

def getXYiono(fname, model):
    f5 = h5.File(fname, 'r')
    if model == 'raiju':    
        xiono = ru.getVar(f5, 'X')*180/np.pi
        yiono = ru.getVar(f5, 'Y')*180/np.pi
        
    elif model == 'mhdrcm':
        # !! Swap x and y so they are the same th, ph as raiju
        yiono = f5['X'][:]
        xiono = 90 - f5['Y'][:]
    
    xicc = kt.to_center2D(xiono)
    yicc = kt.to_center2D(yiono)
    f5.close()

    return {'X': xiono, 'Y': yiono, 'Xcc': xicc, 'Ycc': yicc}


def getStepData(s5, model, iono):
    if model == 'raiju':
        xmin = ru.getVar(s5, 'xmin')
        ymin = ru.getVar(s5, 'ymin')
        #Ni = ymin.shape[0] - 1
        #Nj = ymin.shape[1] - 1
        topo = ru.getVar(s5, 'topo')
        active = ru.getVar(s5, 'active')
        mask_corner = topo != ru.topo['CLOSED']
        mask_cc = active != ru.domain['ACTIVE']
        
        press = ru.getVar(s5, 'Pressure', mask=mask_cc,broadcast_dims=(2,))[:,:,0]
        
    elif model == 'mhdrcm':
        xmin = s5['xMin'][:]
        ymin = s5['yMin'][:]
        #Nj = ymin_rcm.shape[0]
        iopen = s5['IOpen'][:]
        mask_cc = iopen > -0.5

        xiono = iono['X']
        mask_corner = np.full(xiono.shape, False)
        for i in range(1,xiono.shape[0]):
            for j in range(1,xiono.shape[1]):
                mask_corner[i,j] = np.any(mask_cc[i-1:i+1,j-1:j+1])
        
        press = s5['P'][:]
        press = np.ma.masked_where(mask_cc, press)

        
    xcc = kt.to_center2D(xmin)
    ycc = kt.to_center2D(ymin)
    rmin = np.sqrt(xmin**2 + ymin**2)
    rcc  = np.sqrt(xcc**2 + ycc**2)

    return {
        'model': model,
        'xmin': xmin,
        'ymin': ymin,
        'xmincc': xcc,
        'ymincc': ycc,
        'rmin': rmin,
        'rcc': rcc,
        'mask_corner': mask_corner,
        'mask_cc': mask_cc,
        'press': press
    }


def plotTP(ax, ionoData, stepData, norm_press, cmap_press, norm_rmin, cmap_rmin, lvls, xlim, ylim):
    Xiono = ionoData['X']
    Yiono = ionoData['Y']
    Xcc = ionoData['Xcc']
    Ycc = ionoData['Ycc']

    ax.pcolormesh(  Xiono, Yiono, stepData['press'], norm=norm_press, cmap=cmap_press)
    if stepData['model']=='raiju':
        cs = ax.contour(Xiono, Yiono, stepData['rmin'], alpha=0.5, levels=lvls,cmap=cmap_rmin, norm=norm_rmin)
    elif stepData['model']=='mhdrcm':
        cs = ax.contour(Xcc, Ycc, stepData['rmin'] ,alpha=0.5, levels=lvls,cmap=cmap_rmin, norm=norm_rmin)
    ax.clabel(cs)
    ru.drawGrid(ax, Xiono, Yiono, mask=stepData['mask_corner'],color='black')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(alpha=0.25)

    
def plotEq(ax, stepData, norm, cmap, xlim, fname, usrMask=None):
    xmin = stepData['xmin']
    ymin = stepData['ymin']
    press = stepData['press']

    if usrMask is not None:
        press = np.ma.masked_where(usrMask, press)

    ax.pcolormesh(xmin, ymin, press, norm=norm, cmap=cmap)
    if stepData['model']=='raiju':
        ru.drawGrid(ax, xmin, ymin, mask=stepData['mask_corner'], alpha=0.4, color='black')  
    elif stepData['model']=='mhdrcm':
        ru.drawGrid(ax, xmin, ymin, mask=stepData['mask_cc'], alpha=0.4, color='black')
    kv.addEarth2D(ax=ax)
    kv.setBndsByAspect(ax, xlim)
    ax.grid(alpha=0.25)

    textstr = fname

    # place a text box in upper left in axes coords
    props = dict(facecolor='white', alpha=0.5)
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=14,
            verticalalignment='top',horizontalalignment='right', bbox=props)


def plotModelComp(r1_info, r1_model, r2_info, r2_model, stepList, outdir):

    kh5.CheckDirOrMake(outdir)
    fname = os.path.join(outdir, "vid.{:0>4d}.png")

    fig = plt.figure(figsize=(18,12))
    nCol = 4
    iCol = nCol-1
    gs = gridspec.GridSpec(2, nCol, wspace=0.15, hspace=0.07, width_ratios=[0.1,0.4,1,0.1])

    norm_press_hot  = kv.genNorm(1e-2,100,doLog=True)
    norm_rmin = kv.genNorm(2, 20)
    lvls_rmin = np.arange(2, 20, 1)

    cmap_hot = 'viridis'
    cmap_rmin=cmr.torch


    xlim_TP = [15, 45]
    ylim_TP = [360,0]
    bndsXeq = [-15, 10]

    axCB = fig.add_subplot(gs[:,0])
    kv.genCB(axCB, norm_rmin, cbT='L [R$_E$]', cM=cmap_rmin, doVert=True)
    axCB.yaxis.set_ticks_position('left')
    axCB.yaxis.set_label_position('left')
    axCB = fig.add_subplot(gs[:,iCol])
    kv.genCB(axCB, norm_press_hot, cbT='Pressure [nPa]', cM=cmap_hot, doVert=True)

    ax_r1 = [fig.add_subplot(gs[0,iCol-2]), \
             fig.add_subplot(gs[0,iCol-1])]
    ax_r2 = [fig.add_subplot(gs[1,iCol-2]), \
             fig.add_subplot(gs[1,iCol-1])]
    

    # Get base data
    r1_io = getXYiono(r1_info.fname, r1_model)
    r2_io = getXYiono(r2_info.fname, r2_model)
    
    r1_f5 = h5.File(r1_info.fname, 'r')
    r2_f5 = h5.File(r2_info.fname, 'r')


    #for n in tqdm(raiI.Nt):
    for n in tqdm(stepList):
        fname_n = fname.format(n)
        if os.path.exists(fname_n):
            continue

        ax_r1[0].clear()
        ax_r1[1].clear()
        ax_r2[0].clear()
        ax_r2[1].clear()

        s5_r1 = r1_f5[r1_info.stepStrs[n]]
        t_r2 = np.abs(r2_info.times - r2_info.times[n]).argmin()
        s5_r2 = r2_f5[r2_info.stepStrs[t_r2]]

        r1_step = getStepData(s5_r1, r1_model, r1_io)
        r2_step = getStepData(s5_r2, r2_model, r2_io)

        # r1 TP
        ax = ax_r1[0]
        plotTP(ax, r1_io, r1_step, norm_press_hot, cmap_hot, norm_rmin, cmap_rmin, lvls_rmin, xlim_TP, ylim_TP)
        
        # r1 Eq
        ax = ax_r1[1]
        plotEq(ax, r1_step, norm_press_hot, cmap_hot, bndsXeq, r1_info.fname)

        # r2 TP
        ax = ax_r2[0]
        plotTP(ax,r2_io, r2_step, norm_press_hot, cmap_hot, norm_rmin, cmap_rmin, lvls_rmin, xlim_TP, ylim_TP)

        # r2 Eq
        ax = ax_r2[1]
        plotEq(ax, r2_step, norm_press_hot, cmap_hot, bndsXeq, r2_info.fname)
        

        time = s5_r1.attrs['time']
        UT = kt.MJD2UT(s5_r1.attrs['MJD'])
        titleStr = "T = {:1.0f} s, UT={}".format(time, UT)
        ax_r1[1].set_title(titleStr)

        kv.savePic(fname_n,dpiQ=350)


if __name__=="__main__":
    indir = '.'
    outdir = 'runComp'
    r1_fname = 'raijuOWD.raiju.h5'
    r1_model = 'raiju'
    r2_fname = 'msphere.mhdrcm.h5'
    r2_model = 'mhdrcm'
    stepStride = 10

    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-d',default=indir,
        help="Directory contianing two output files (default: %(default)s)"
    )
    parser.add_argument(
        '-out',default=outdir,
        help="Directory to store output plots (default: rundir+'%(default)s')"
    )
    parser.add_argument(
        '-r1f',default=r1_fname,
        help="Run 1 fname (default: %(default)s)"
    )
    parser.add_argument(
        '-r1m',default=r1_model,
        help="Run 1 model name ('raiju' or 'mhdrcm', default: %(default)s)"
    )
    parser.add_argument(
        '-r2f',default=r2_fname,
        help="Run 2 fname (default: %(default)s)"
    )
    parser.add_argument(
        '-r2m',default=r2_model,
        help="Run 2 model name ('raiju' or 'mhdrcm', default: %(default)s)"
    )
    parser.add_argument(
        '-s',default=stepStride,
        help="step stride (default: %(default)s)"
    )
    args = parser.parse_args()
    indir = args.d
    outdir = args.out
    r1_fname = args.r1f
    r1_model = args.r1m
    r2_fname = args.r2f
    r2_model = args.r2m
    stepStride = int(args.s)


    r1_fname = os.path.join(indir, r1_fname)
    f2_fname = os.path.join(indir, r2_fname)

    idx_psph = ru.flavs_s['PSPH']+1
    idx_e    = ru.flavs_s['HOTE']+1
    idx_p    = ru.flavs_s['HOTP']+1

    print("Getting model 1 info")
    r1_info = getInfo(r1_fname, r1_model)
    print("Getting model 2 info")
    r2_info = getInfo(r2_fname, r2_model)

    print("Steps:",r1_info.Nt)
    stepList = np.arange(0,r1_info.Nt, stepStride)
    plotModelComp(r1_info, r1_model, r2_info, r2_model, stepList, outdir)