import h5py as h5
import numpy as np
import os
import progressbar
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.time import Time
import datetime

import kaipy.kaiH5 as kh5
import kaipy.kaiViz as kv
import kaipy.kaijson as kj
import kaipy.kaiTools as kT

import kaipy.satcomp.scutils as scutils


#Optionals
doProgressBar = True

#Constants
massi = 1.67e-27 # mass of ions in kg
masse = 9.11e-31 # mass of lectrons in kg
ev = 1.607e-19 # 1ev in J
nt = 1.0e-9 # nt
re = 6.380e6 # radius of earth in m

pressure_factor   = 2./3.*ev/re*nt
specFlux_factor_i = 1/np.pi/np.sqrt(8)*np.sqrt(ev/massi)*nt/re/1.0e1  # [units/cm^2/keV/str]
specFlux_factor_e = 1/np.pi/np.sqrt(8)*np.sqrt(ev/masse)*nt/re/1.0e1  # [units/cm^2/keV/str]
TINY = 1.0e-8

#Hard-coded filenames
RCM_TIME_JFNAME = 'rcm_times.json'
RCM_TKL_JFNAME = 'rcm_tkl.json'
RCM_EQLATLON_JFNAME = 'rcm_eqlatlon.json'
MHDRCM_TIME_JFNAME = 'mhdrcm_times.json'

#Spacecraft strings for cdaweb retrieval
SC_str = {
	'RBSP': {
	 'E_RBSPICE' : {
	  'DsetName': "RBSP-%s-RBSPICE_LEV-3_ESRHELT",
	  'DsetVar': ["FEDU", "MLT", 'L_Eq', 'L_Star'],
	  'dfStr': "FEDU",
	  'nrgStr': "FEDU_Energy_Fixed",
	  'epochStr': 'Epoch'}, 
	 "E-PAP_RBSPICE" : {
	  'DsetName': "RBSP-%s-RBSPICE_LEV-3-PAP_ESRHELT",
	  'DsetVar': ["FEDU", "FEDU_OmniFlux", "MLT", 'L_Eq', 'L_Star'],
	  'dfStr': "FEDU_OmniFlux",
	  'nrgStr': "FEDU_Energy",
	  'epochStr': 'Epoch'},
	 "H_RBSPICE": {
	  'DsetName': "RBSP-%s-RBSPICE_LEV-3_TOFXEH",
	  'DsetVar': ["bad"],  # Leave this here so things fail correctly
	  'dfStr': "FPDU",
	  'nrgStr': "FPDU_Energy",
	  'epochStr': 'Epoch'},
	 "H-PAP_RBSPICE": {
	  'DsetName': "RBSP-%s-RBSPICE_LEV-3-PAP_TOFXEH",
	  'DsetVar': ["bad"],  # Leave this here so things fail correctly
	  'dfStr': "FPDU_OmniFlux",
	  'nrgStr': "FPDU_Energy",
	  'epochStr': 'Epoch'},
	 'E_HOPE' : {
	  'DsetName': 'RBSP%s_REL04_ECT-HOPE-PA-L3',
	  'DsetVar': ["FEDU, MLT_Ele"],
	  'dfStr': 'FEDU',
	  'nrgStr': 'HOPE_ENERGY_Ele',
	  'epochStr': 'Epoch_Ele'},
	 'H_HOPE': {
	  'DsetName': 'RBSP%s_REL04_ECT-HOPE-PA-L3',
	  'DsetVar': ["FPDU, MLT_ION"],
	  'dfStr': 'FPDU',
	  'nrgStr': 'HOPE_ENERGY_Ion',
	  'epochStr': 'Epoch_Ion'},
	 'EPH': {
	 'DsetName': "RBSP-%s_MAGNETOMETER_1SEC-GSM_EMFISIS-L3",
	 'DsetVar': "coordinates"}
	}
}

#======
#Helpers
#======

#Generate json filename for given spacecraft dataset
def genSCD_jfname(jdir, scName, dSetName):
	jfname = scName + "_" + dSetName + ".json"
	return os.path.join(jdir, jfname)
def genRCMTrack_jfname(jdir, scName):
	jfname = "rcmTrack_" + scName + ".json"
	return os.path.join(jdir, jfname)

#Return a 3D cube of a given variable to use in interpolation
def getVarCube(rcm5, varName, stepLow, stepHigh, ilon, ilat, k=None):
	if k is None:
		v1 = np.expand_dims(rcm5[stepLow  ][varName][ilon:ilon+2, ilat:ilat+2], axis=2)
		v2 = np.expand_dims(rcm5[stepHigh][varName][ilon:ilon+2, ilat:ilat+2], axis=2)
	else:
		v1 = np.expand_dims(rcm5[stepLow  ][varName][k, ilon:ilon+2, ilat:ilat+2], axis=2)
		v2 = np.expand_dims(rcm5[stepHigh][varName][k, ilon:ilon+2, ilat:ilat+2], axis=2)
	return np.append(v1, v2, axis=2)

#electrons: kStart = kStart, kEnd = kIon
#ions: kStart = kIon, kEnd = len(rcmS0['alamc'])
def getSpecieslambdata(rcmS0, kStart, kEnd):
	ilamc = rcmS0['alamc'][kStart:kEnd]  # Cell centers
	Nk = len(ilamc)
	ilami = np.zeros(Nk+1)  # Cell interfaces
	for n in range(0, Nk-1):
		ilami[n+1] = 0.5*(ilamc[n]+ilamc[n+1])
	ilami[Nk] = ilamc[-1] + 0.5*(ilamc[-1]-ilamc[-2])
	
	ilamc = np.abs(ilamc)
	ilami = np.abs(ilami)
	lamscl = np.diff(ilami)*np.sqrt(ilamc)

	result = {'ilamc' : ilamc,
				'ilami' : ilami,
				'lamscl' : lamscl}
	return result

def get_aspect(ax):
		from operator import sub
		# Total figure size
		figW, figH = ax.get_figure().get_size_inches()
		# Axis size on figure
		_, _, w, h = ax.get_position().bounds
		# Ratio of display units
		disp_ratio = (figH * h) / (figW * w)
		# Ratio of data units
		# Negative over negative because of the order of subtraction
		data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())

		return disp_ratio #/ data_ratio

#======
#Main work
#======

#Given sc and dataset name (according to above dict), grab specifically omnidirecitonal differential flux
def getSCOmniDiffFlux(scName, dSetName, t0, t1, jdir=None):
	if jdir is not None:
		dojson = True
		jfname = genSCD_jfname(jdir, scName, dSetName)
		if os.path.exists(jfname):
			print("Grabbing spacecraft data from " + jfname)
			data = kj.load(jfname)
			ephdata = data['ephdata']
			dataset = data['dataset']
			return ephdata, dataset
	else:
		dojson = False

	print("Pulling spacecraft data from cdaweb")
	#TODO: Add all desired datasets
	if 'RBSP' in scName:
		#Extract RBSP identifier (A or B)
		scTag = scName.split('RBSP')[1]
		if '-' in scTag:
			scTag = scName.split('-')[1]

		#Pull data
		ephStrs = SC_str['RBSP']['EPH']
		dsStrs = SC_str['RBSP'][dSetName]
		ephdata = scutils.getCdasData(ephStrs['DsetName']%(scTag), ephStrs['DsetVar'], t0, t1)
		data = scutils.getCdasData(dsStrs['DsetName']%(scTag), dsStrs['DsetVar'], t0, t1)
		#data = getCdasData(SC_str['RBSP']['E-PAP_RBSPICE'])
		
		dfStr = dsStrs['dfStr']
		nrgStr = dsStrs['nrgStr']
		epochStr = dsStrs['epochStr']
		
		dataset = {}
		dataset['name'] = scName
		species = 'electrons' if 'E' == dSetName[0] else 'ions'
		dataset['species'] = species
		dataset['epoch'] = data[epochStr]
		#Turn each dataset's data into omniflux
		if dSetName == 'E-PAP_RBSPICE' or dSetName == 'H-PAP_RBSPICE':
			#Already got omni flux, no problem
			dataset['OmniDiffFlux'] = data[dfStr]*1E-3  # Diferential flux [1/(MeV-cm^2-s-sr]*[1/MeV -> 1/keV]
			dataset['energies'] = data[nrgStr]*1E3  # [MeV] -< [keV]

	#Pause to save to json
	if dojson:
		print("Saving to file")
		jdata = {'ephdata' : ephdata, 'dataset' : dataset}
		kj.dump(jfname, jdata)

	return ephdata, dataset

def getRCMtimes(rcmf5,mhdrcmf5,jdir=None):
	"""Grab RCM times, sIDs, and MJDs
		If jdir given, will try to find the files there
		If not found there, will pull the data from the hdf5's
		If those aren't provided, we can't do much
	"""
	if jdir is not None:
		dojson = True
		rcmjfname = os.path.join(jdir, RCM_TIME_JFNAME)
		mhdrcmjfname = os.path.join(jdir, MHDRCM_TIME_JFNAME)
	else:
		dojson = False

	if dojson:
		if os.path.exists(rcmjfname):
			print("Grabbing RCM time data from " + rcmjfname)
			#That's all that's needed, done
			return kj.load(rcmjfname)
		else:
			#Must create RCM MJD's based on MHDRCM file, so get MHDRCM data
			if os.path.exists(mhdrcmjfname):
				print("Grabbing MHDRCM time data from " + mhdrcmjfname)
				mhdrcmTimes = kj.load(mhdrcmjfname)
				mhdrcmT = mhdrcmTimes['T']
				mhdrcmMJDs = mhdrcmTimes['MJD']
			else:
				print("No usable files in " + jdir + ", grabbing all from hdf5's")

	if 'mhdrcmTimes' not in locals():
		print("Grabbing MHDRCM time data from " + mhdrcmf5)
		Nt, sIDs = kh5.cntSteps(mhdrcmf5)
		sIDs = np.sort(sIDs)
		mhdrcmT = kh5.getTs(mhdrcmf5, sIDs, aID='time')
		mhdrcmMJDs = kh5.getTs(mhdrcmf5, sIDs, aID='MJD')

		mhdrcmTimes = {'Nt': Nt,
						'sIDs' : sIDs,
						'T' : mhdrcmT,
						'MJD' : mhdrcmMJDs}
		#Pause to save to file
		if dojson: kj.dump(mhdrcmjfname, mhdrcmTimes)

	print("Grabbing RCM time data from " + rcmf5)
	Nt,sIDs = kh5.cntSteps(rcmf5)
	sIDs = np.sort(sIDs)
	sIDstrs = np.array(['Step#'+str(s) for s in sIDs])
	rcmT = kh5.getTs(rcmf5,sIDs,aID="time")

	if (mhdrcmT[:] == rcmT[:]).all():
		rcmMJDs = mhdrcmMJDs
	else:
		# !! This only works if all rcm steps are also in mhdrcm
		#    In the future: As long as a single timestep is in both, that MJD can be used to get all MJDs for rcm steps
		rcmMJDs = np.zeros((len(rcmT)))
		for i in range(len(rcmT)):
			idx = np.where(mhdrcmT == rcmT[i])[0]
			rcmMJDs[i] = mhdrcmMJDs[idx]

	rcmTimes = {'Nt' : Nt,
				'sIDs' : sIDs,
				'sIDstrs' : sIDstrs,
				'T' : rcmT,
				'MJD' : rcmMJDs}
	#Pause to save rcm jfile
	if dojson: kj.dump(rcmjfname, rcmTimes)

	return rcmTimes

#TODO: Add scName to trackfile attrs so we can pull directly from there
def getRCM_scTrack(trackf5, rcmf5, rcmTimes, jdir=None, scName=""):
	"""Pull RCM data along a given spacecraft track
		trackfile: big spacecraft trajectory hdf5, generated from sctrack.x
		rcmf5: <tag>.rcm.h5 file
		jdir: If included, will try to do json saving and loading to save time

		returns: dictionary containing along track: time, mjd, mlat, mlon, vm, e and i energy and eetas
	"""
	if jdir is not None:
		dojson = True
		jfname = genRCMTrack_jfname(jdir, scName)
	else:
		dojson = False

	if dojson:
		if os.path.exists(jfname):
			print("Grabbing RCM track data from " + jfname)
			return kj.load(jfname)

	print("Extracting RCM track data from " + rcmf5)
	kh5.CheckOrDie(trackf5)
	scMLATs = kh5.PullVar(trackf5, 'MLAT')
	scMLONs = kh5.PullVar(trackf5, 'MLON')
	scTs = kh5.PullVar(trackf5, 'T')
	scMJDs = kh5.PullVar(trackf5, 'MJDs')
	Nsc = len(scTs)

	#Get information for mirror ratio
	Bx = kh5.PullVar(trackf5,"Bx")
	By = kh5.PullVar(trackf5,"By")
	Bz = kh5.PullVar(trackf5,"Bz")
	Bmag = np.sqrt(Bx**2.0 + By**2.0 + Bz**2.0)
	Beq = kh5.PullVar(trackf5,"Beq")

	J0 = scutils.getJScl(Bmag,Beq)

	#Unpack what we need from rcmTimes
	sIDstrs = rcmTimes['sIDstrs']
	rcmMJDs = rcmTimes['MJD']

	#Init rcm h5 info
	rcm5 = h5.File(rcmf5,'r')
	rcmS0 = rcm5[sIDstrs[0]]
	Ni, Nj = rcmS0['aloct'].shape
	rcmMLAT = 90.0-rcmS0['colat'][0,:]*180/np.pi
	rcmMLON = rcmS0['aloct'][:,0]*180/np.pi
	rcmMLAT_min = np.min(rcmMLAT)
	rcmMLAT_max = np.max(rcmMLAT)

	kStart = 1 if rcmS0['alamc'][0] == 0 else 0  # Check channel 0 for plasmasphere
	kIon = (rcm5[sIDstrs[0]]['alamc'][kStart:]>0).argmax()+kStart

	#Init electron and ion dicts
	sdata = {}
	sdata['electrons'] = getSpecieslambdata(rcmS0, kStart, kIon               )
	sdata['ions'     ] = getSpecieslambdata(rcmS0, kIon  , len(rcmS0['alamc']))
	Nk_e = kIon - kStart
	Nk_i = len(rcmS0['alamc']) - kIon
	
	#Collect data along spacecraft track
	vms = np.zeros((Nsc))
	xmin = np.zeros((Nsc))
	ymin = np.zeros((Nsc))
	zmin = np.zeros((Nsc))
	energies_e = np.zeros((Nsc, Nk_e))
	energies_i = np.zeros((Nsc, Nk_i))
	eeta_e = np.zeros((Nsc, Nk_e))
	eeta_i = np.zeros((Nsc, Nk_i))
	diffFlux_e = np.zeros((Nsc, Nk_e))
	diffFlux_i = np.zeros((Nsc, Nk_i))
	
	if doProgressBar: bar = progressbar.ProgressBar(max_value=Nsc)

	for n in range(Nsc):
		if doProgressBar: bar.update(n)

		mjd_sc = scMJDs[n]
		mlat_sc = scMLATs[n]
		mlon_sc = scMLONs[n]
		#Make sure track and rcm domain overlap
		if mjd_sc < rcmMJDs[0] or mjd_sc > rcmMJDs[-1] or \
			mlat_sc < rcmMLAT_min or mlat_sc > rcmMLAT_max:
			continue

		# Get bounds in rcm space
		ilat = len(rcmMLAT)-1  # mlat_rcm goes from high to low
		while rcmMLAT[ilat] < mlat_sc: ilat -= 1
		ilon = 2
		while rcmMLON[ilon+1] < mlon_sc: ilon += 1
		imjd = 0
		while rcmMJDs[imjd+1] < mjd_sc: imjd += 1

		latbnd = [rcmMLAT[ilat], rcmMLAT[ilat+1]]
		lonbnd = [rcmMLON[ilon], rcmMLON[ilon+1]]
		mjdbnd = [rcmMJDs[imjd], rcmMJDs[imjd+1]]		
		stepLow = sIDstrs[imjd]
		stepHigh = sIDstrs[imjd+1]

		vmcube = getVarCube(rcm5, 'rcmvm', stepLow, stepHigh, ilon, ilat)
		vms[n] = scutils.trilinterp(lonbnd, latbnd, mjdbnd, vmcube, mlon_sc, mlat_sc, mjd_sc)
		#Do the same for xeq, yeq, zeq
		xmincube = getVarCube(rcm5, 'rcmxmin', stepLow, stepHigh, ilon, ilat)
		ymincube = getVarCube(rcm5, 'rcmymin', stepLow, stepHigh, ilon, ilat)
		zmincube = getVarCube(rcm5, 'rcmzmin', stepLow, stepHigh, ilon, ilat)
		xmin[n] = scutils.trilinterp(lonbnd, latbnd, mjdbnd, xmincube, mlon_sc, mlat_sc, mjd_sc)
		ymin[n] = scutils.trilinterp(lonbnd, latbnd, mjdbnd, ymincube, mlon_sc, mlat_sc, mjd_sc)
		zmin[n] = scutils.trilinterp(lonbnd, latbnd, mjdbnd, zmincube, mlon_sc, mlat_sc, mjd_sc)

		def getSpecEetas(kOffset, Nk):
			eetas = np.zeros((Nk))
			for k in range(Nk):
				kr = k + kOffset
				eetacube = getVarCube(rcm5, 'rcmeeta', stepLow, stepHigh, ilon, ilat, kr)
				eetas[k] = scutils.trilinterp(lonbnd, latbnd, mjdbnd, eetacube, mlon_sc, mlat_sc, mjd_sc)
			return eetas

		eeta_e[n,:] = getSpecEetas(kStart, Nk_e)
		eeta_i[n,:] = getSpecEetas(kIon, Nk_i)
		energies_e[n,:] = vms[n]*sdata['electrons']['ilamc']
		energies_i[n,:] = vms[n]*sdata['ions']['ilamc']

		diffFlux_e[n,:] = J0[n]*specFlux_factor_e*energies_e[n,:]*eeta_e[n,:]/sdata['electrons']['lamscl']
		diffFlux_i[n,:] = J0[n]*specFlux_factor_i*energies_i[n,:]*eeta_i[n,:]/sdata['ions'     ]['lamscl']

	# Package everything together
	sdata['electrons']['energies'] = energies_e*1E-3  # [eV -> keV]
	sdata['electrons']['eetas'   ] = eeta_e
	sdata['electrons']['diffFlux'] = diffFlux_e
	sdata['ions'     ]['energies'] = energies_i*1E-3  # [eV -> keV]
	sdata['ions'     ]['eetas'   ] = eeta_i
	sdata['ions'     ]['diffFlux'] = diffFlux_i
	
	result = {}
	result['T'        ] = scTs
	result['MJD'      ] = scMJDs
	result['MLAT'     ] = scMLATs
	result['MLON'     ] = scMLONs
	result['vm'       ] = vms
	result['xmin'     ] = xmin
	result['ymin'     ] = ymin
	result['zmin'     ] = zmin
	#TODO: Map to full equator using zmin I think
	result['eqmin'    ] = np.sqrt(xmin**2+ymin**2)
	result['electrons'] = sdata['electrons']
	result['ions'     ] = sdata['ions']

	if dojson: kj.dump(jfname, result)

	return result

#TODO: Energy grid mapping in a nice, jsonizable way
#      Right now, just need to call this whenever you want it
def consolidateODFs(scData, rcmTrackData, eGrid=None, jdir=None):
	"""Prepare the spacecraft and rcm track data for comparison
		Match up energy grids, save all the needed info in one place
	"""
	#scData determines which species we're using
	species = scData['species']
	rcmSpec = rcmTrackData[species]


	scEGrid = scData['energies']  # Might be 1D (fixed bins) or 2D (different energies over time)
	rcmEGrid = rcmSpec['energies']

	#Manipulate sc EGrid, if 2D, so that time dimension is first index
	sce_shape = scEGrid.shape
	if len(sce_shape) > 1:
		timeAxis = -1
		
		for i in range(len(shape)):
			if sce_shape[i] == Nt_sc:
				timeAxis = i
		if timeAxis == -1:
			print("Finding right time axis in sc data didn't work")
			return
		scEGrid = np.rollaxis(scEGrid, timeAxis, 0)
		
	#If given grid to use, do that
	#Else, make fixed bins spanning full energy range
	if eGrid is None:
		eMax = np.max([scEGrid.max(), rcmEGrid.max()])
		eMin = np.min([scEGrid[scEGrid>0].min(), rcmEGrid[rcmEGrid>0].min()])
		numPoints = 150
		eGrid = np.logspace(np.log10(eMin), np.log10(eMax), numPoints, endpoint=True)

	Nt_sc = len(scData['epoch'])
	sc_odf = np.zeros((Nt_sc, len(eGrid)))
	for n in range(Nt_sc):
		#Might not need this check
		if len(sce_shape) > 1:
			sc_odf[n,:] = scutils.varMap_1D(scEGrid[n,:], eGrid, scData['OmniDiffFlux'][n,:])
		else:
			sc_odf[n,:] = scutils.varMap_1D(scEGrid, eGrid, scData['OmniDiffFlux'][n,:])

	Nt_rcm = len(rcmTrackData['T'])
	rcm_odf = np.zeros((Nt_rcm, len(eGrid)))
	for n in range(Nt_rcm):
		rcm_odf[n,:] = scutils.varMap_1D(rcmEGrid[n,:], eGrid, rcmSpec['diffFlux'][n,:])

	result = {}
	result['energyGrid'] = eGrid
	result['sc'] = {
		'name' : scData['name'],
		'time' : scData['epoch'],
		'diffFlux' : sc_odf}
	result['rcm'] = {
		'time' : rcmTrackData['MJD'],
		'origEGrid' : rcmEGrid,
		'origODF' : rcmSpec['diffFlux'],
		'diffFlux' : rcm_odf}

	return result

def getIntensitiesVsL(rcmf5, mhdrcmf5, sStart, sEnd, sStride, species='ions', eGrid=None, AxLvT=None, jdir=None, joverwrite=True):
	"""Calculate rcm intensities (summed diff flux)
		rcmf5: rcm hdf5 filename to pull xmin, ymin, zmin from
		mhdrcmf5: mhdrcm hdf5 filename to pull IOpen from
		AxLvT: If given, will plot the resulting L shell vs. time intensity
		jdir: Give json directory to enable json usage (read/write results to file)
		joverwrite: If dataset already found in file, re-calc anyways and overwrite it
	"""

	if jdir is not None:
		dojson = True
		rcmjfname = os.path.join(jdir, RCM_TKL_JFNAME)
	else:
		dojson = False

	if dojson:
		if os.path.exists(rcmjfname):
			print("Grabbing RCM tkl data from " + rcmjfname)
			return kj.load(rcmjfname)

	#Setup

	rcmTimes = getRCMtimes(rcmf5,mhdrcmf5,jdir=jdir)
	iTStart = np.abs(rcmTimes['sIDs']-sStart).argmin()
	iTEnd = np.abs(rcmTimes['sIDs']-sEnd).argmin()

	sIDstrs = rcmTimes['sIDstrs'][iTStart:iTEnd+1:sStride]
	nSteps = len(sIDstrs)
	rcm5 = h5.File(rcmf5,'r')
	rcmS0 = rcm5[sIDstrs[0]]
	mhdrcm5 = h5.File(mhdrcmf5,'r')


	kStart = 1 if rcmS0['alamc'][0] == 0 else 0  # Check channel 0 for plasmasphere
	kIon = (rcm5[sIDstrs[0]]['alamc'][kStart:]>0).argmax()+kStart
	if species == 'electrons':
		kEnd = kIon
		sf_factor = specFlux_factor_e
	elif species == 'ions':
		kStart = kIon
		kEnd = len(rcm5[sIDstrs[0]]['alamc'])-1
		sf_factor = specFlux_factor_i
	alamData = getSpecieslambdata(rcmS0, kStart, kEnd)

	alams_kxx = alamData['ilamc'][:, np.newaxis, np.newaxis]
	Nk = kEnd - kStart

	nlbins = 50
	lbins = np.linspace(2, 15, nlbins)

	if eGrid is None:
		eMin = 1E2  # [eV]
		eMax = 1E6  # [eV]
		Ne = 50
		eGrid = np.logspace(np.log10(eMin), np.log10(eMax), Ne, endpoint=True)
	Ne = len(eGrid)

	#Add endpoints to catch values outside of desired ranges
	lbins = np.concatenate(([lbins[0]-TINY], lbins, [lbins[-1]+TINY]))
	eGrid = np.concatenate(([eGrid[0]-TINY], eGrid, [eGrid[-1]+TINY]))

	rcmodf_tkl = np.zeros((nSteps, Ne+2, nlbins+2))
	rcmpress_tkl = np.zeros((nSteps, Ne+2, nlbins+2))

	if doProgressBar: bar = progressbar.ProgressBar(max_value=nSteps)
	for n in range(nSteps):
		if doProgressBar: bar.update(n)

		rcmS = rcm5[sIDstrs[n]]
		mhdrcmS = mhdrcm5[sIDstrs[n]]

		IOpen = mhdrcmS['IOpen']
		Ni, Nj = IOpen.shape  # Default should be (179, 90)
		#Shorten rcm data to match Ni, Nj of mhdrcm5
		vms = rcmS['rcmvm'][2:, :]
		xmins = rcmS['rcmxmin'][2:,:]
		ymins = rcmS['rcmymin'][2:,:]
		zmins = rcmS['rcmzmin'][2:,:]
		eetas = rcmS['rcmeeta'][kStart:kEnd,2:,:]

		vms_xij = vms[np.newaxis,:,:]

		#Calculate L shell for whole plane
		L_arr = scutils.xyz_to_L(xmins, ymins, zmins)  # [i,j]

		le_counts = np.zeros((Ne+2, nlbins+2))  # Keep track of how many entries into each bin so we can divide by it to get the average later

		energies = vms_xij * alams_kxx  # Should be [Nk, Ni, Nj]

		iL_arr = np.array([np.abs(lbins-i).argmin() for i in L_arr.flatten()]).reshape((Ni, Nj))
		iE_arr = np.array([np.abs(eGrid-e).argmin() for e in energies.flatten()]).reshape((Nk, Ni, Nj))

		#diffFlux_Nk = sf_factor*energies*eetas/alamData['lamscl'][:, np.newaxis, np.newaxis]  # [k,i,j]

		pressure_kij = pressure_factor*alams_kxx*eetas*vms_xij**2.5 * 1E9  # [Pa -> nPa]

		for i in range(Ni):
			for j in range(Nj):
				for k in range(Nk):
					#rcmodf_tkl[n, iE_arr[k,i,j], iL_arr[i,j]] += diffFlux_Nk[k,i,j]
					rcmpress_tkl[n, iE_arr[k,i,j], iL_arr[i,j]] += pressure_kij[k,i,j]
					le_counts[iE_arr[k,i,j], iL_arr[i,j]] += 1

		#Normalize to get avg. pressure per count
		rcmpress_tkl[n,:,:] /= le_counts

	#Trim off the extra values
	lbins = lbins[1:-1]
	eGrid = eGrid[1:-1]
	#rcmodf_tkl = rcmodf_tkl[:,1:-1,1:-1]
	rcmpress_tkl = rcmpress_tkl[:,1:-1,1:-1]
	
	# Sum across energy to get total avg. pressures per L shell
	rcmpress_tl = np.ma.sum(np.ma.masked_invalid(rcmpress_tkl), axis=1)


	result = {}
	result['T']	         = rcmTimes['T'  ][iTStart:iTEnd+1:sStride]
	result['MJD']        = rcmTimes['MJD'][iTStart:iTEnd+1:sStride]
	result['lambda']     = alamData['ilamc']
	result['L_bins']     = lbins
	result['energyGrid'] = eGrid
	#result['nrg_tkl'] = rcmnrg_tkl
	#result['odf_tkl']    = rcmodf_tkl
	result['press_tkl']  = rcmpress_tkl
	result['press_tl']  = rcmpress_tl

	if dojson: kj.dump(rcmjfname, result)

	return result

#TODO: Take list of variable strings to pull from rcm.h5 file
def getRCM_eqlatlon(mhdrcmf5, rcmTimes, jdir=None):
	"""Grab certain variables along with equatorial and lat-lon grid
		Can use json but there's not much point if you already have the (mhd)rcm file(s)
	"""
	if jdir is not None:
		dojson = True
		rcmjfname = os.path.join(jdir, RCM_EQLATLON_JFNAME)
	else:
		dojson = False

	if dojson:
		if os.path.exists(rcmjfname):
			print("Grabbing RCM eq_lat-lon data from " + rcmjfname)
			return kj.load(rcmjfname)

	Nt = rcmTimes['Nt']
	sIDs = rcmTimes['sIDs']
	sIDstrs = rcmTimes['sIDstrs']

	mhdrcm5 = h5.File(mhdrcmf5,'r')
	
	mhdrcmS0 = mhdrcm5[sIDstrs[0]]
	Ni,Nj = mhdrcmS0['aloct'].shape

	mlat_rcm = 90.0-mhdrcmS0['colat'][0,:]*180/np.pi
	mlon_rcm = mhdrcmS0['aloct'][:,0]*180/np.pi
	mlatrcm_min = np.min(mlat_rcm)
	mlatrcm_max = np.max(mlat_rcm)

	#Get desired variables for all timesteps
	xmin_arr = np.ma.zeros((Nt, Ni, Nj))
	ymin_arr = np.ma.zeros((Nt, Ni, Nj))
	press_arr = np.ma.zeros((Nt, Ni, Nj))
	dens_arr = np.zeros((Nt, Ni, Nj))
	
	scLoc_eq = np.zeros((Nt, 2))
	scLoc_latlon = np.zeros((Nt, 2))

	rMin = 1.25
	rMax = 35.0
	ioCut = -0.5
	pCut = 1E-8
	
	print("Grabbing data...")
	
	for t in range(len(sIDstrs)):
		xm = mhdrcm5[sIDstrs[t]]['xMin'][:]
		ym = mhdrcm5[sIDstrs[t]]['yMin'][:]
		bmR = np.sqrt(xm*xm + ym*ym)
		pm = mhdrcm5[sIDstrs[t]]['P'][:]
		#Turn IOpen into a big true/false map
		iopen_t = mhdrcm5[sIDstrs[t]]['IOpen'][:] 

		Ir = (bmR < rMin) | (bmR > rMax)
		I_m = Ir | (iopen_t > ioCut) | (pm < pCut)
		
	import kaipy.gamera.gampp as gampp
	import kaipy.gamera.rcmpp as rcmpp
	#import kaipy.gamera.msphViz as mviz
	rcmdata = gampp.GameraPipe('',mhdrcmf5.split('.h5')[0])
	for t in range(1, len(sIDs)):
		bmX, bmY = rcmpp.RCMEq(rcmdata, sIDs[t], doMask=True)
		I = rcmpp.GetMask(rcmdata, sIDs[t])
		pm = rcmpp.GetVarMask(rcmdata, sIDs[t], 'P', I)

		xmin_arr[t,:,:] = np.transpose(bmX)
		ymin_arr[t,:,:] = np.transpose(bmY)
		press_arr[t,:,:] = np.transpose(pm)

	"""
	for t in range(tStart, len(sIDstrs)):
		#Linterp sc location based on time
		rcmTime = rcm_ut[t]
		if ut[0] > rcmTime or ut[-1] < rcmTime:
			scLoc_eq[t,:] = [0,0]
			scLoc_latlon[t,:] = [0,0]
		else:
			itime = 0
			while ut[itime] < rcmTime: itime += 1

			scLoc_eq[t,:] = [xeq_sc[itime], yeq_sc[itime]]
			scLoc_latlon[t,:] = [mlat_sc[itime], mlon_sc[itime]]
	"""
	result = {}
	result['T']     = rcmTimes['T']
	result['MJD']   = rcmTimes['MJD']
	result['MLAT']  = mlat_rcm
	result['MLON']  = mlon_rcm
	result['xmin']  = xmin_arr
	result['ymin']  = ymin_arr
	result['press'] = press_arr

	if dojson: kj.dump(rcmjfname, result)

	return result

#======
#Plotting
#======

def plt_ODF_Comp(AxSC, AxRCM, AxCB, odfData, mjd=None, cmapName='CMRmap', norm=None):
	axIsPopulated = not AxSC.get_ylabel() == ''

	eGrid   = odfData['energyGrid']
	scTime  = odfData['sc']['time']
	scODF   = odfData['sc']['diffFlux']
	rcmTime = odfData['rcm']['time']
	rcmODF  = odfData['rcm']['diffFlux']

	#ut = scutils.mjd_to_ut(rcmTime)
	ut = kT.MJD2UT(rcmTime)
	
	if norm is None:
		vMax = np.max([scODF.max(), rcmODF.max()])
		vMin = np.max([scODF[scODF>0].min(), rcmODF[rcmODF>0].min()])
		norm = kv.genNorm(vMin,vMax,doLog=True)


	if not axIsPopulated:
		kv.genCB(AxCB,norm,r'Differential Flux [$cm^{-2} sr^{-1} s^{-1} keV^{-1}$]',cM=cmapName,doVert=True)

		AxSC.pcolormesh(scTime, eGrid, np.transpose(scODF), norm=norm, shading='nearest', cmap=cmapName)
		AxSC.set_xlim([ut[0], ut[-1]])
		AxSC.set_ylabel("%s Energy [keV]"%odfData['sc']['name'])
		AxSC.set_yscale('log')
		AxSC.xaxis.set_major_formatter(plt.NullFormatter())

		AxRCM.pcolormesh(ut, eGrid, np.transpose(rcmODF), norm=norm, shading='nearest', cmap=cmapName)
		#for n in range(len(rcmTime)):
		#	AxRCM.plot(odfData['rcm']['origEGrid'][n,:], odfData['rcm']['origODF'][n,:])
		AxRCM.set_ylabel("RCM Energy [keV]")
		AxRCM.set_yscale('log')


	if mjd is not None:
		if mjd < rcmTime[0] or mjd > rcmTime[-1]:
			print(str(mjd) + "not in rcm data, exiting")
			return
		iMJD = np.abs(rcmTime - mjd).argmin()
		lineUT = ut[iMJD]

		if len(AxSC.lines) != 0:
			AxSC.lines.pop(0)  #Remove previous mjd line
			AxRCM.lines.pop(0)
		yMin, yMax = AxSC.get_ylim()
		AxSC.plot([lineUT, lineUT], [yMin, yMax], '-k')
		yMin, yMax = AxRCM.get_ylim()
		AxRCM.plot([lineUT, lineUT], [yMin, yMax], '-k')
		
def plt_tkl(AxTL, AxTKL, AxCB, tkldata, mjd=None, cmapName='CMRmap', norm=None):
	"""If 'mjd' is not provided, make all plots that are vs. time
	   If 'mjd' is provided:
	     If we were also given a populated AxTL and AxCB, update with an mjd scroll line
	     Also generate AxTKL for this mjd step
	"""
	tlIsPopulated = not AxTL.get_ylabel() == ''

	k_arr = tkldata['energyGrid']
	L_arr = tkldata['L_bins']

	press_tl = np.array(tkldata['press_tl'][:], dtype=float)  # Need to do this to handle masked stuff
	press_tl = np.ma.masked_invalid(press_tl)
	
	#ut = scutils.mjd_to_ut(tkldata['MJD'])
	ut = kT.MJD2UT(tkldata['MJD'])
	#utstr = [t.strftime('%m-%d\n%H') for t in ut]

	if norm is None:
		vMin = np.min(press_tl[press_tl>0])
		vMax = np.max(press_tl)
		print(vMin)
		print(vMax)
		norm = kv.genNorm(vMin, vMax, doLog=True)


	#Initialize static plots if hasn't been done yet
	if not tlIsPopulated:
		#L vs. Time
		AxTL.pcolormesh(ut, L_arr, np.transpose(press_tl), norm=norm, shading='nearest', cmap=cmapName)
		AxTL.set_xlim([ut[0], ut[-1]])
		AxTL.set_ylabel('L shell')
		AxTL.set_xlabel('UT')
		kv.genCB(AxCB, norm, r'Total pressure [$nPa$]', cM=cmapName, doVert=False)

	#L vs. k for a specific mjd
	if mjd is not None:
		if mjd < tkldata['MJD'][0] or mjd > tkldata['MJD'][-1]:
			print(str(mjd) + "not in tkl data, exiting")
			return
		iMJD = np.abs(tkldata['MJD'] - mjd).argmin()
		lineUT = ut[iMJD]

		if len(AxTL.lines) != 0:
			AxTL.lines.pop(0)  #Remove previous mjd line
		yMin, yMax = AxTL.get_ylim()
		AxTL.plot([lineUT, lineUT], [yMin, yMax], '-k')


		#energy_arr = tkldata['nrg_tkl'][]
		klslice = tkldata['press_tkl'][iMJD,:,:]
		AxTKL.pcolormesh(k_arr*1E-3, L_arr, np.transpose(klslice), norm=norm, shading='nearest', cmap=cmapName)
		#AxTKL.pcolormesh(k_arr, L_arr, np.transpose(klslice), shading='nearest', cmap=cmapName)
		AxTKL.set_xlabel('Energy [keV]')
		AxTKL.set_ylabel('L shell')

def plt_rcm_eqlatlon(AxLatlon, AxEq, rcmData, satTrackData, AxCB=None, mjd=None, norm=None, cmapName='viridis'):
	
	mjd_arr   = rcmData['MJD']
	xmin_arr  = rcmData['xmin']
	ymin_arr  = rcmData['ymin']
	mlat_arr  = rcmData['MLAT']
	mlon_arr  = rcmData['MLON']
	press_arr = rcmData['press']

	x_sc  = satTrackData['xmin']
	y_sc  = satTrackData['ymin']
	z_sc  = satTrackData['zmin']
	eq_sc = satTrackData['eqmin']

	#ut = scutils.mjd_to_ut(rcmData['MJD'])
	ut = kT.MJD2UT(rcmData['MJD'])
	Nt = len(ut)

	if norm is None:
		vMin = np.min(press_arr[press_arr>0])
		vMax = np.max(np.ma.masked_invalid(press_arr))
		print(vMin)
		print(vMax)
		norm = kv.genNorm(vMin, vMax, doLog=True)


	#Initialize static plots if hasn't been done yet
	if AxCB is not None:
		AxCB = kv.genCB(AxCB, norm, r'Pressure [$nPa$]', cM=cmapName, doVert=True)

	if mjd is not None:
		if mjd < mjd_arr[0] or mjd > mjd_arr[-1]:
			print(str(mjd) + "not in rcm data, exiting")
			return
		iMJD = np.abs(mjd_arr - mjd).argmin()
		lineUT = ut[iMJD]

		AxLatlon.clear()
		AxEq.clear()

		#Prep rcm lat/lons for polar plotting
		riono = np.cos(mlat_arr*np.pi/180.)
		tiono = np.concatenate((mlon_arr, [mlon_arr[0]]))*np.pi/180.
		AxLatlon.pcolor(tiono, riono, np.transpose(press_arr[iMJD]),norm=norm, shading='auto', cmap=cmapName)
		AxLatlon.axis([0, 2*np.pi, 0, 0.7])

		AxEq.pcolor(xmin_arr[iMJD], ymin_arr[iMJD], press_arr[iMJD], norm=norm, shading='auto', cmap=cmapName)

		#Draw satellite location
		if eq_sc[iMJD] > 1E-8:
			leadMax = iMJD
			while leadMax < min(iMJD+80, Nt) and eq_sc[leadMax] > 1E-8: leadMax += 1 #????
			AxEq.plot(x_sc[iMJD:leadMax], y_sc[iMJD:leadMax], 'k-')
			
			satCircle = plt.Circle((x_sc[iMJD], y_sc[iMJD]), 0.15, color='black')
			AxEq.add_patch(satCircle)
		kv.addEarth2D(ax=AxEq)

		#Set bounds
		boxRatio = get_aspect(AxEq)
		xbMin = -10
		xbMax = 10
		dx = xbMax - xbMin
		dy = boxRatio*dx
		ybMin = -dy/2
		ybMax = dy/2
		AxEq.set_xlim([xbMin, xbMax])
		AxEq.set_ylim([ybMin, ybMax])
