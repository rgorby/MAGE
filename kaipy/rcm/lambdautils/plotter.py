import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

import kaipy.kaiTools as kT


def plotLambdas_Val_Spac(specDataList, yscale='log',L=None):
	doEnergy = True if L is not None else False
	if doEnergy: 
		bVol = kT.L_to_bVol(L)
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
		bVol = kT.L_to_bVol(L)
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
		plt.suptitle("Energy grid values and spacing\n@ L={}, bVol={:2.2f}, vm={:2.2f}".format(L, bVol, vm))
	plt.show()

def plotEnergyRange(specDataList, rInner=1.5, rOuter=15, rRes=100):
	"""Plot energy range of each species as a function of L shell
	"""
	rVals = np.linspace(rInner, rOuter, rRes)
	bVols = np.array([kT.L_to_bVol(r) for r in rVals])
	vms = bVols**(-2/3)

	allMin = 1e8
	allMax = -1
	specPlotList = []
	for i in range(len(specDataList)):
		specData = specDataList[i]
		if specData.name != "Plasmasphere":
			specPlotList.append(specData)
			allMin = np.min((allMin, np.abs(specData.alams[0])*vms[-1]*1E-3))
			allMax = np.max((allMax, np.abs(specData.alams[-1])*vms[0]*1E-3))
	nSpecs = len(specPlotList)

	allMax


	fig = plt.figure(figsize=(10,5))
	gs = gridspec.GridSpec(1,nSpecs)

	for i in range(nSpecs):	
		specData = specPlotList[i]

		lamMin = np.abs(specData.alams[0])
		lamMax = np.abs(specData.alams[-1])
		eMins  = lamMin*vms*1E-3  # [eV -> keV]
		eMaxs  = lamMax*vms*1E-3  # [eV -> keV]

		Ax = fig.add_subplot(gs[:,i])
		Ax.fill_between(rVals, eMins, eMaxs)
		Ax.set_yscale('log')
		Ax.set_xlabel('L Shell [R$_p$]')
		Ax.set_ylabel('Energy [keV]')
		Ax.set_ylim([0.8*allMin, 1.2*allMax])
		Ax.grid(True)
		Ax.set_title(specData.name)
	plt.suptitle("Covered enery range w.r.t. dipole L-shell")
	plt.show()


