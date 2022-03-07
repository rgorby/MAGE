import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

def L_to_bVol(L):  # L shell [Re] to V [Re/nT]
    bsurf_nT = 3.11E4
    colat = np.arcsin(np.sqrt(1.0/L))

    cSum = 35*np.cos(colat) - 7*np.cos(3*colat) +(7./5.)*np.cos(5*colat) - (1./7.)*np.cos(7*colat)
    cSum /= 64.
    s8 = np.sin(colat)**8
    V = 2*cSum/s8/bsurf_nT
    return V


def plotLambdas_Val_Spac(specDataList, yscale='log',L=None):
	doEnergy = True if L is not None else False
	if doEnergy: 
		bVol = L_to_bVol(L)
		vm = bVol**(-2/3)
	if not isinstance(specDataList, list):
		specDataList = [specDataList]

	fig = plt.figure(figsize=(10,5))
	gs = gridspec.GridSpec(1,2)
	AxAlams = fig.add_subplot(gs[:,0])
	AxDiffs = fig.add_subplot(gs[:,1])

	for specData in specDataList:
		chNum = np.array([i for i in range(specData.n)])
		if doEnergy:
			energies = np.abs(specData.alams)*vm*1E-3  # [ev -> keV]
			energyDiff = np.diff(energies)
			AxAlams.step(chNum, energies, label=specData.name)
			AxDiffs.step(chNum[:-1], energyDiff, label=specData.name)
		else:
			alams = np.abs(specData.alams)
			alamDiff = np.diff(alams)
			AxAlams.step(chNum, alams, label=specData.name)
			AxDiffs.step(chNum[:-1], alamDiff, label=specData.name)
	AxAlams.set_yscale(yscale)
	AxAlams.legend()
	AxAlams.set_xlabel('Channel Number')
	AxAlams.title.set_text("Values")

	AxDiffs.set_yscale(yscale)
	AxDiffs.legend()
	AxDiffs.set_xlabel('Channel Number')
	AxDiffs.title.set_text("Spacing")

	if doEnergy:
		AxAlams.set_ylabel("E [keV]")
		AxDiffs.set_ylabel(r"$\Delta$E [keV]")
		plt.suptitle("L={} vm={:2.2f}".format(L, vm))
	else:
		AxAlams.set_ylabel(r"$\lambda$")
		AxDiffs.set_ylabel(r"$\Delta\lambda$")
	plt.show()

def plotLambdasBySpec(specDataList, yscale='log',L=None):
	doEnergy = True if L is not None else False
	if doEnergy:
		bVol = L_to_bVol(L)
		vm = bVol**(-2/3)
	if not isinstance(specDataList, list):
		specDataList = [specDataList]
	specPlotList = []
	for i in range(len(specDataList)):
		if specDataList[i].name != "Plasmasphere":
			specPlotList.append(specDataList[i])
			
	nSpecs = len(specPlotList)

	fig = plt.figure(figsize=(10,5))
	gs = gridspec.GridSpec(1,nSpecs)
	#AxAlams = fig.add_subplot(gs[:,0])
	#AxDiffs = fig.add_subplot(gs[:,1])

	for i in range(nSpecs):
		specData = specPlotList[i]
		chNum = np.array([i for i in range(specData.n)])
		Ax = fig.add_subplot(gs[:,i])
		if doEnergy:
			energies = np.abs(specData.alams)*vm*1E-3  # [ev -> keV]
			energyDiff = np.diff(energies)
			Ax.step(chNum, energies, label='Values')
			Ax.step(chNum[:-1], energyDiff, label="Spacing")
			Ax.set_ylabel("E [keV]")
		else:
			alams = np.abs(specData.alams)
			alamDiff = np.diff(alams)
			Ax.step(chNum, alams, label="Values")
			Ax.step(chNum[:-1], alamDiff, label="Spacing")
			Ax.set_ylabel(r"$\lambda$")
		Ax.set_yscale(yscale)
		Ax.legend()
		Ax.grid(True)
		Ax.set_xlabel('Channel Number')
		Ax.title.set_text(specData.name)
	if doEnergy:
		plt.suptitle("L={} vm={:2.2f}".format(L, vm))
	plt.show()