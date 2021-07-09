import numpy as np
import h5py as h5

from kaipy.rcm.lambdautils.AlamData import AlamData
from kaipy.rcm.lambdautils.helperdefs import L_to_bVol 
import kaipy.rcm.lambdautils.dpRangeEval as dpRangeEval
#Constants
radius_earth_m = 6.371E6
ev = 1.6e-19  # electron volt

# Enter testing through one of these
# They all just use the input data to make an AlamData object, then send to testAlams
def testAlam_arr(alam_arr):
	alamdict = {}
	iproton = len(alam_arr)-1
	while alam_arr[iproton-1] > 0:
		iproton -= 1
	alamdict['spec1'] = alam_arr[:iproton]  # Electrons
	alamdict['spec2'] = alam_arr[iproton:]  # Assuming sincle ion species
	testAlam_dict(alamdict)

def testAlam_dict(alamdict):
	alamdata = AlamData(alamdict=alamdict)
	testAlams(alamdata)

def testAlam_config(fname):
	alamdata = AlamData(filename=fname)
	testAlams(alamdata)

def testAlam_AD(alamdata):
	testAlams(alamdata)

# Actual alam testing start
def testAlams(alamdata, LSmearTest=2.5, LRangeTest=10, ratioTolerance=0.02):
	kap = 5
	smearResult = smearTest(alamdata, L=LSmearTest)
	maxwellRangeMin, maxwellRangeMax, maxwell_vsq_d, maxwell_vsq_p = ktRangeTest(alamdata, ratioTolerance, L=LRangeTest)
	print(ktRangeTest(alamdata, ratioTolerance, L=LRangeTest, kap=kap))
	kapRangeMin, kapRangeMax, kap_vsq_d, kap_vsq_p = ktRangeTest(alamdata, ratioTolerance, L=LRangeTest, kap=kap)

	if smearResult < 1:
		resultstr = "Passed"
		
	else:
		resultstr = "Failed"
	print("Smear test at L = {}: {}".format(LSmearTest, resultstr))
	print("  Worst smear/cellWidth: {:2.2e}".format(smearResult))
	print("kT min/max range within {:3.1f}% tolerance:".format(ratioTolerance*100))
	print("  Maxwellian : {:1.2e}/{:1.2e} [keV]".format(maxwellRangeMin, maxwellRangeMax))
	print("    Variance: D = {:1.2e} P = {:1.2e}".format(maxwell_vsq_d, maxwell_vsq_p))
	print("  kappa = {:3.0f}: {:1.2e}/{:1.2e} [keV]".format(kap, kapRangeMin, kapRangeMax))
	print("    Variance: D = {:1.2e} P = {:1.2e}".format(kap_vsq_d, kap_vsq_p))


"""
	Quantify "smearing" of energy channel by calculating
	the difference in drift velocities between lamba+ and 
	lambda- values over some time dt.

	If smear > individual cell width, lambda channel width
	is too large.

	Inputs: 
		alamdata: AlamData object
		L       : L shell in R_e
		dt      : time-delta in seconds
		deg     : longitudinal RCM resolution in

	Outputs:
		Pass/Fail value
		  0 = Pass, non-zero = maxSmear/cellWidth > 1
"""
def smearTest(alamdata, L=2.5, dt=10, deg=1, doVerbose=False):
	# Get width of a grid cell in meters
	cellWidth = (L*radius_earth_m)*(np.pi/180*deg)
	if doVerbose: print("Cell width: {:1.2e}".format(cellWidth))
	bVol = L_to_bVol(L)
	vm = bVol**(-2./3.)
	if doVerbose: print("bVol = {}, vm = {}".format(bVol, vm))

	smeardict = {}
	maxSmear = 0
	for s in alamdata.alams.keys():
		alams = alamdata.alams[s]
		amins = alamdata.amins[s]
		amaxs = alamdata.amaxs[s]

		smear_arr = np.array([])

		#  Drift velocity ~ W_k*delB/q/B^2
		bMag = 3.03E-5*(L)**(-3)  # [T]
		delB = 9.09E-5*(L)**(-4)/radius_earth_m  # [T/m]
		const = delB/ev/bMag**2

		driftMin = const*alams[1]*vm*ev
		driftMax = const*alams[-1]*vm*ev
		if doVerbose: print("Drift min/max ({}): {:1.2e}/{:1.2e}".format(s, driftMin, driftMax))

		for i in range(len(alams)):

			wMinus = amins[i]*vm*ev  # [J]
			wPlus = amaxs[i]*vm*ev  # [J]

			smear = const*(wPlus-wMinus) * dt  # [m]

			smear_arr = np.append(smear_arr, smear)
			if smear > maxSmear:
				maxSmear = smear

			energy = alams[i]*vm*1E-3
			drift = const*alams[i]*vm*ev
			if(doVerbose):
				print("Energy = {:1.2E} [keV] Drift: {:1.2e} [m/s] Smear: {:1.2e} [m] Smear/CellWidth: {:1.2e}".format(energy, drift, smear, smear/cellWidth))

	#if maxSmear/cellWidth < 1:
	#	return 0  # Pass
	#else:
	return maxSmear/cellWidth  # Return ratio of worst offender

# Similar setup in main of dpRangeEval.py
# Calculate X'/X for range ot kT's, then determine min/max kT within tolerance range
def ktRangeTest(alamdata, tol, L=10, kap=-1):

	# Defaults for now
	nSamples = 100
	ktmin = 1E-3 # [keV]
	ktmax = 1000 # [keV]
	dconst = 10 # [1/cc]
	kt_arr = np.logspace(np.log10(ktmin), np.log10(ktmax), nSamples)
	dens_arr = np.ones((nSamples))*dconst*1E6 # [1/m^3]
	press_arr = np.array([kt*dconst/6.25 * 1E-9 for kt in kt_arr]) # [Pa]
	bVol = L_to_bVol(L)
	vm = bVol**(-2./3.)

	dRatios, pRatios = dpRangeEval.getRatios1D(alamdata, dens_arr, press_arr, vm, kap)

	# This next part is not efficient
	
	ktTolMin = kt_arr[0]
	ktTolMax = kt_arr[-1]
	# Find point closest to 1keV
	iStart = np.abs(kt_arr - 1).argmin()
	if (np.abs(dRatios[iStart] - 1) > tol) or (np.abs(pRatios[iStart] - 1) > tol):
		return -1, -1, -1, -1
	
	# Scan to find kT values where X'/X falls out of tolerance range,
	iTolMin = 0 
	iTolMax = 0
	iScan = iStart
	while iScan < len(kt_arr):
		if (np.abs(dRatios[iScan] - 1) > tol) or (np.abs(pRatios[iScan] - 1) > tol): 
			iTolMax = iScan - 1
			ktTolMax = kt_arr[iTolMax]
			iScan = 2*len(kt_arr)
		else:
			iScan += 1
	iScan = iStart
	while iScan > -1:
		if (np.abs(dRatios[iScan] - 1) > tol) or (np.abs(pRatios[iScan] - 1) > tol): 
			iTolMin = iScan + 1
			ktTolMin = kt_arr[iTolMin]
			iScan = -1
		else:
			iScan -= 1
	# Make a metric for how much the ratios deviate from 1 within the tolerance range
	var_sq_d = 0
	var_sq_p = 0
	for k in range(iTolMin, iTolMax):
		var_sq_d += np.abs(dRatios[k]-1)**2
		var_sq_p += np.abs(pRatios[k]-1)**2
	var_sq_d /= (iTolMax-iTolMin)
	var_sq_p /= (iTolMax-iTolMin)

	return ktTolMin, ktTolMax, var_sq_d, var_sq_p




