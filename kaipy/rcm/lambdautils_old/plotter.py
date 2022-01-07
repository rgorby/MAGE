import matplotlib
defaultBackend = matplotlib.get_backend()  # Save whatever the default loaded backend is for future GUI use
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import ticker, colors
import numpy as np

from kaipy.rcm.lambdautils.helperdefs import L_to_bVol

titlesize = 16
ylabelsize = 22
xlabelsize = 20
xlabelpad = 8
ylabelpad = 8
ticklabelsize = 16
xticklabelsize = 15
legendsize = 20

def makeShowable():
    matplotlib.use(defaultBackend)
    import matplotlib.pyplot as plt
    from matplotlib import ticker, colors

def plotChannels(alams, doShow=False):
    titlesize = 22
    ylabelsize = 22
    xlabelsize = 20
    xlabelpad = 8
    ylabelpad = 8
    ticklabelsize = 16
    xticklabelsize = 15
    legendsize = 20

    if(doShow): makeShowable()
    fig, ax = plt.subplots(1)

    # Divide alams by e and i species
    iproton = len(alams)-1
    while alams[iproton-1] > 0:
        iproton -= 1
    alams_e = alams[:iproton]
    alams_p = alams[iproton:]
    k_arr_e = np.array([k for k in range(len(alams_e))])
    k_arr_p = np.array([k for k in range(len(alams_p))])
    aSpace_e = alams_e[1:]-alams_e[:-1]
    aSpace_p = alams_p[1:]-alams_p[:-1]

    # Plot alam values first
    ax.plot(k_arr_p, alams_p, label="ions")
    ax.plot(k_arr_e, np.abs(alams_e), label="electrons")
    #ax.plot(k_arr_p[:-1], aSpace_p, '--')
    #ax.plot(k_arr_e[:-1], np.abs(aSpace_e), '--')
    ax.set_xlabel("Channel #", fontsize=xlabelsize)
    ax.set_ylabel("$\\lambda$", fontsize = ylabelsize)
    #ax.set_yscale('log')
    ax.legend(prop={'size': legendsize})
    ax.grid(True)
    ax.tick_params(axis="x", labelsize=ticklabelsize)
    ax.tick_params(axis="y", labelsize=ticklabelsize)
    plt.suptitle("Energy Channel Distribution", fontsize=titlesize)
    

    if(doShow): plt.show()


def plotEtas(alamData, etaData, vm, title, doShow=False):

    if(doShow): makeShowable()
    fig, ax = plt.subplots()
        
    s = "spec2"  # species label for protons
    legendText = "Protons"

    widths = (alamData.amaxs[s]-alamData.amins[s])
    mids = (alamData.amaxs[s]+alamData.amins[s])/2.0
    ax.bar(mids, etaData[s]/widths, widths*0.9, label=legendText)
    #for i in range(len(alamData.alams[s])):
    #    ax.plot([alamData.alams[s][i], alamData.alams[s][i]], [0, etaData[s][i]/widths[i]], 'r-')
    
    ax.set_xlabel("Energy Invariant $\lambda$", fontsize=xlabelsize)
    ax.set_ylabel("$\\frac{\\eta}{\\Delta \\lambda}$", fontsize=28)
    ax.set_ylim((1E15, 5E17))
    #ax.set_yscale('log')
    ax.set_title(title, fontsize=titlesize)

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(ticklabelsize)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(ticklabelsize)
    
    if(doShow): plt.show()

    return ax

def plotLambdaVsLShell(alamData, Llow=1.5, Lhigh=10, doShow=False):
    nSamples = 100
    lshell_arr = np.linspace(Llow, Lhigh, nSamples, endpoint=True)
    vm_arr = np.array([L_to_bVol(L) for L in lshell_arr])**(-2/3)
    iLams_e = np.array([0, len(alamData.alams['spec1'])-1])
    if np.abs(alamData.alams['spec1'][0]) < 1E-5:
        iLams_e[0] += 1
    iLams_p = np.array([0, len(alamData.alams['spec2'])-1])

    lams_e = np.array([np.abs(alamData.alams['spec1'][i]) for i in iLams_e])
    lams_p = np.array([alamData.alams['spec2'][i] for i in iLams_p])
    print(lams_e)
    energies_e = np.zeros((len(iLams_e), nSamples))
    energies_p = np.zeros((len(iLams_p), nSamples))

    print(energies_e[0])
    for i in range(nSamples):
        #print(energies_e[:][:])
        energies_e[:,i] = lams_e*vm_arr[i]
        energies_p[:,i] = lams_p*vm_arr[i]

    print(energies_e)
    if(doShow): makeShowable()
    fig, ax = plt.subplots()

    for e_e in energies_e:
        line_e, = ax.plot(lshell_arr, e_e, color="red")
    for e_p in energies_p:
        line_p, = ax.plot(lshell_arr, e_p, color="blue")

    ax.set_yscale('log')
    ax.set_xlabel("L shell [$R_E$]", fontsize=xlabelsize)
    ax.set_ylabel("Channel Energy [eV]", fontsize=ylabelsize)
    ax.tick_params(axis="x", labelsize=ticklabelsize)
    ax.tick_params(axis="y", labelsize=ticklabelsize)
    ax.grid(True)
    ax.legend([line_e, line_p], ["Electron min/max", "Ion min/max"])
    plt.suptitle("RCM Channel energy range vs. L shell")
    

    if(doShow): plt.show()

    return ax

def plotRatios_tempRange1D(kt_arr, dRatios, pRatios, tag, fracChannelUse=None, doShow=False):
    
    #Calc plot bounds based on tolerance of pressure
    itolMin = 0
    itolMax = len(dRatios)-1
    
    if(doShow): makeShowable()
    fig, ax = plt.subplots()

    ax.plot(kt_arr[itolMin:itolMax], dRatios[itolMin:itolMax], label="D'/D")
    ax.plot(kt_arr[itolMin:itolMax], pRatios[itolMin:itolMax], label="P'/P")
    ax.plot([kt_arr[itolMin], kt_arr[itolMax]], [1,1], 'k--')
    ax.set_xlabel("kT [keV]", fontsize=xlabelsize)
    ax.set_ylabel("X'/X", fontsize=ylabelsize)
    ax.set_xscale("log")
    ax.set_xlim((kt_arr[itolMin], kt_arr[itolMax]))
    ax.set_ylim((0.90, 1.1))
    ax.legend(loc='upper left', prop={'size': legendsize})

    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(ticklabelsize)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(ticklabelsize)

    if fracChannelUse is not None:
        ax2=ax.twinx()
        ax2.plot(kt_arr, fracChannelUse, '-.', color="red", label="frac. channel use")
        ax2.set_ylabel("Fraction of channel use", color="red", fontsize=ylabelsize)
        #ax2.set_ylim((0,1))
        ax2.tick_params(axis="y", colors="red", labelsize=ticklabelsize)
        #for tick in ax2.yaxis.get_major_ticks():
        #    tick.label.set_fontsize(ticklabelsize)
    
    plt.suptitle(tag, fontsize=titlesize)
    

    if(doShow): plt.show()

    return ax, ax2

def plotRatios_dpRange2D(dens_arr, press_arr, dRatios, pRatios, tag, doShow=False):
    
    #dd, pp = np.meshgrid(dens_arr, press_arr, indexing='ij')
    dd, pp = np.meshgrid(dens_arr, press_arr)
    
    if(doShow): makeShowable()
    fig, (ax1, ax2) = plt.subplots(1,2)

    lvl_min = .85
    lvl_max = 1.15

    levels = np.linspace(lvl_min, lvl_max, 10)
    level_ticks = np.linspace(lvl_min, lvl_max, 5, endpoint=True)
    level_tick_str = ["{:1.2f}".format(x) for x in level_ticks]
    
    CS_d = ax1.contourf(dd*1E-6, pp*1E9, dRatios, levels, cmap=plt.cm.plasma)
    cbar_d = fig.colorbar(CS_d)
    cbar_d.set_ticks(level_ticks)
    cbar_d.set_ticklabels(level_tick_str)
    ax1.set_xlabel("Density [1/cc]")
    ax1.set_ylabel("Pressure [nPa]")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_title("D'/D")

    CS_p = ax2.contourf(dd*1E-6, pp*1E9, pRatios, levels, cmap=plt.cm.plasma)
    #cbar_p = fig.colorbar(CS_p)
    #cbar_p.set_ticks(level_ticks)
    #cbar_p.set_ticklabels(level_tick_str)
    ax2.set_xlabel("Density [1/cc]")
    ax2.set_ylabel("Pressure [nPa]")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_title("P'/P")

    plt.suptitle(tag)
    

    if(doShow): plt.show()

    return ax1, ax2

#def superPlot(filename=None):
#    fig, axs = plt.subplots()

#Graveyard

    """ log levels
    #lev_exp = np.linspace(np.floor(np.log10(z.min())-1), np.ceil(np.log10(z.max())+1), 100)
        lev_exp = np.linspace(np.floor(np.log10(lvl_min)), np.ceil(np.log10(lvl_max)), 100)
        levels = np.power(10, lev_exp)

        #lev_tick_exp = np.linspace(np.floor(np.log10(z.min())-1), np.ceil(np.log10(z.max())+1), 5)
        lev_tick_exp = np.linspace(np.floor(np.log10(lvl_min)), np.ceil(np.log10(lvl_max)), 5)
        level_ticks = np.power(10, lev_tick_exp)
        level_tick_str = ["{:1.2e}".format(x) for x in level_ticks]

        CS_d = ax.contourf(dd*1E-6, pp*1E9, z, levels, norm=colors.LogNorm(), cmap=plt.cm.plasma)
    """

#Adjust x range for ratio_1D plot
    """
    tol_min = 0.75
    tol_max = 1.25
    for i in range(len(kt_arr)-1):
        #Try to capture transitions from in to out of tolerance, and visa-versa
        if (pRatios[i] > tol_min and pRatios[i+1] < tol_min) \
            or (pRatios[i] < tol_min and pRatios[i+1] > tol_min) \
            or (pRatios[i] > tol_max and pRatios[i+1] < tol_max) \
            or (pRatios[i] < tol_max and pRatios[i+1] > tol_max):
                if itolMin == 0:
                    itolMin = i
                    print("set itolMin to" + str(itolMin))
                    break
    for i in range(len(kt_arr)-2, 0, -1):
        if (pRatios[i] > tol_min and pRatios[i+1] < tol_min) \
            or (pRatios[i] < tol_min and pRatios[i+1] > tol_min) \
            or (pRatios[i] > tol_max and pRatios[i+1] < tol_max) \
            or (pRatios[i] < tol_max and pRatios[i+1] > tol_max):
                itolMax = i
                print("set itolMax to" + str(itolMax))
                break
    """
    #itolMax = len(dRatios)-1
    #print("{} (dR = {}, pR = {})".format(itolMin, dRatios[itolMin], pRatios[itolMin]))
    #print("{} (dR = {}, pR = {})".format(itolMax, dRatios[itolMax], pRatios[itolMax]))
