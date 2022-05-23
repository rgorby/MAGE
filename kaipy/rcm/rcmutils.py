import numpy as np
import h5py as h5

import kaipy.kdefs as kd


#------
# Factors
#------
massi = kd.Mp_cgs*1e-3 # mass of ions in kg
masse = kd.Me_cgs*1e-3 # mass of lectrons in kg
nt    = 1.0e-9 # nt
ev    = kd.eCharge # 1ev in J
rx    = kd.Re_cgs*1e-2 # radius of earth in m

pressure_factor   = 2./3.*ev/rx*nt
specFlux_factor_i = 1/np.pi/np.sqrt(8)*np.sqrt(ev/massi)*nt/rx/1.0e1  # [units/cm^2/keV/str]
specFlux_factor_e = 1/np.pi/np.sqrt(8)*np.sqrt(ev/masse)*nt/rx/1.0e1  # [units/cm^2/keV/str]
TINY = 1.0e-8

def updateFactors(rxNew):
	global rx
	global pressure_factor
	global specFlux_factor_i
	global specFlux_factor_e

	rx = rxNew  # [m]

	pressure_factor   = 2./3.*ev/rx*nt
	specFlux_factor_i = 1/np.pi/np.sqrt(8)*np.sqrt(ev/massi)*nt/rx/1.0e1  # [units/cm^2/keV/str]
	specFlux_factor_e = 1/np.pi/np.sqrt(8)*np.sqrt(ev/masse)*nt/rx/1.0e1  # [units/cm^2/keV/str]

#------
# Lambda Channels
#------

#electrons: kStart = kStart, kEnd = kIon
#ions: kStart = kIon, kEnd = len(rcmS0['alamc'])
def getSpecieslambdata(rcmS0, species='ions'):
	#Determine what our bounds are
	kStart = 1 if rcmS0['alamc'][0] == 0 else 0  # Check channel 0 for plasmasphere
	kIon = (rcmS0['alamc'][kStart:]>0).argmax()+kStart
	if species == 'electrons':
		kEnd = kIon
	elif species == 'ions':
		kStart = kIon
		kEnd = len(rcmS0['alamc'])

	ilamc = rcmS0['alamc'][kStart:kEnd]  # Cell centers
	Nk = len(ilamc)
	ilami = np.zeros(Nk+1)  # Cell interfaces
	for n in range(0, Nk-1):
		ilami[n+1] = 0.5*(ilamc[n]+ilamc[n+1])
	ilami[Nk] = ilamc[-1] + 0.5*(ilamc[-1]-ilamc[-2])
	
	ilamc = np.abs(ilamc)
	ilami = np.abs(ilami)
	lamscl = np.diff(ilami)*np.sqrt(ilamc)

	result = {	'kStart': kStart,
				'kEnd' : kEnd,
				'ilamc' : ilamc,
				'ilami' : ilami,
				'lamscl' : lamscl}
	return result

#------
# Domain
#------

def getClosedRegionMask(s5):

	rmin = np.sqrt(s5['rcmxmin'][:]**2 + s5['rcmymin'][:]**2)

	mask = np.full(s5['rcmxmin'].shape, False)
	mask = (s5['rcmvm'][:] <= 0) | (rmin < 1.0)

	return mask

#------
# Cumulative pressure
#------
def getCumulPress(ilamc, vm, eetas, doFraction=False):
	""" Returns 3D array of cumulative presures
		ilamc: [Nk]       lambda values (Nk NOT the full number of energy channels, can be subset)
		vm   : [Nj,Ni]    vm value
		eetas: [Nk,Nj,Ni] eetas corresponding to ilamc values
		doFraction: return cumulative pressure fraction instead of cumulative pressure
	"""

	Nk = len(ilamc)
	Nj,Ni = vm.shape
	pCumul = np.zeros((Nk,Nj,Ni))

	ilam_kji = ilamc[:,np.newaxis,np.newaxis]
	vm_kji = vm[np.newaxis, :, :]

	pPar = pressure_factor*ilam_kji*eetas*vm_kji**2.5 * 1E9  # [Pa -> nPa], partial pressures for each channel

	pCumul[0] = pPar[0] # First bin's partial press = cumulative pressure
	for k in range(1,Nk):
		pCumul[k] = pCumul[k-1] + pPar[k]  # Get current cumulative pressure by adding this bin's partial onto last bin's cumulative
	
	if not doFraction:
		return pCumul
	else:
		return pCumul/pCumul[-1]

def getValAtLoc_linterp(val, xData, yData, getAxis='y'):
	""" Use linear interpolation to find desired x/y values in data
		val: x/y value to find y/x location of
		xData/yData: 1D arrays of equal length
		getAxis: desired axis. If 'y', assumes targets are x values, and visa versa
	"""
	
	idx = 0
	if getAxis == 'y': # target is x-axis value, find its location in xData
		while xData[idx+1] < val: idx += 1
	elif getAxis == 'x': # target is y-axis value, find its location in yData
		while yData[idx+1] < val: idx += 1
	m = (yData[idx+1] - yData[idx])/(xData[idx+1] - xData[idx])
	b = yData[idx] - m*xData[idx]
	if getAxis == 'y':
		loc = m*val + b
	elif getAxis == 'x':
		loc = (val - b)/m

	return loc





