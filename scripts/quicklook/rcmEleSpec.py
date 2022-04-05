#This is mostly a copy/paste from miscwork/quicklook/satTrackSlice.py

import numpy as np
from numpy import ma
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import h5py as h5
import argparse
import os
import progressbar
import datetime

import kaipy.kaiViz as kv
import kaipy.kaiTools as kT
import kaipy.kaiH5 as kh5
import kaipy.satcomp.scRCM as scRCM

#internal stuff
FLAG = 1E-8
CF_PRESS_I = "cf_Pi"
CF_PRESS_E = "cf_Pe"
C_PRESS_I = "c_Pi"
C_PRESS_E = "c_Pe"
ENERGY_I = "energy_i"
ENERGY_E = "energy_e"

#Physics stuff
EarthMag = 0.31*1.0e-4  # [nT]

to_center = lambda A: 0.25*(A[:-1,:-1]+A[1:,:-1]+A[:-1,1:]+A[1:,1:])

def plt_rcmeq(AxEq, rcmData, xPnts, yPnts, AxCB=None, mjd=None, norm=None, cmapName='viridis'):
	
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

	mjd_arr   = rcmData['MJD']
	xmin_arr  = rcmData['xmin']
	ymin_arr  = rcmData['ymin']
	press_arr = rcmData['press']
	#press_arr = rcmData['pressE']
	
	#ut = scutils.mjd_to_ut(rcmData['MJD'])
	ut = kT.MJD2UT(rcmData['MJD'])
	Nt = len(ut)

	if norm is None: 
		norm = genVarNorm(press_arr, doLog=True)

	#Initialize static plots if hasn't been done yet
	if AxCB is not None:
		AxCB = kv.genCB(AxCB, norm, r'Pressure [$nPa$]', cM=cmapName, doVert=False)

	if mjd is not None:
		if mjd < mjd_arr[0] or mjd > mjd_arr[-1]:
			print(str(mjd) + "not in rcm data, exiting")
			return
		iMJD = np.abs(mjd_arr - mjd).argmin()
		lineUT = ut[iMJD]

		AxEq.clear()
		AxEq.pcolor(xmin_arr[iMJD], ymin_arr[iMJD], to_center(press_arr[iMJD]), norm=norm, shading='auto', cmap=cmapName)

		kv.addEarth2D(ax=AxEq)

		#Draw probe points
		for xPnt,yPnt in zip(xPnts,yPnts):
			satCircleB = plt.Circle((xPnt, yPnt), 0.12, color='black',zorder=10)
			satCircleW = plt.Circle((xPnt, yPnt), 0.18, color='white',zorder=10)
			AxEq.add_patch(satCircleW)
			AxEq.add_patch(satCircleB)

			r = np.sqrt(xPnt**2+yPnt**2)
			t = Ax.text(xPnt,yPnt+0.05,"{}".format(r), horizontalalignment='center', \
				verticalalignment='center',fontsize=18)
			t.set_bbox(dict(facecolor='white',edgecolor='black',alpha=0.85))

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



def plotTSCResult(x2D,y2D,pntX,pntY,vWeights,vGrid,pf):

	plt.scatter(x2D,y2D)
	plt.scatter(pntX,pntY)
	plt.xlabel('X [R]')
	plt.ylabel('Y [R]')
	plt.text(pntX,pntY+0.005, '({:4.2f},{:4.2f})'.format(pntX,pntY), fontsize=9)
	plt.text(pntX,pntY-0.01,'Pf={:4.2f}'.format(pf), fontsize=9)
	for i in range(3):
		for j in range(3):
			plt.text(x2D[i,j],y2D[i,j]+0.005,'w: {:4.3f}'.format(vWeights[i,j]), fontsize=8)
			plt.text(x2D[i,j],y2D[i,j]-0.01,'P: {:4.2f}'.format(vGrid[i,j]), fontsize=8)
	plt.show()

def interpTSCWeights(gridX, gridY, x, y):
	""" gridX/gridY: 3-element x & y grid vals (center value is closest point to desired lat/lon)
	 	var: 3x3 values of desired variable
	 	x: dim1 of point of interest
	 	y: dim2 of point of interest
	"""

	dx = np.abs(gridX[0]-gridX[1])
	dy = np.abs(gridY[0]-gridY[1])

	eta  = (x - gridX[1])/dx
	zeta = (y - gridY[1])/dy
	#print('eta : ' + str(eta))
	#print('zeta: ' + str(zeta))

	def weight1D(eta):
		return np.array([
				0.5*(0.5-eta)**2,
				0.75 - eta**2,
				0.5*(0.5+eta)**2
				
			])

	wX = weight1D(eta)
	wY = weight1D(zeta)

	w2D = np.zeros((3,3))
	for i in range(3):
		for j in range(3):
			w2D[i,j] = wX[i]*wY[j]

	return w2D

def getVarDims(s5, mxstr, vName):
	
	#Easy part, figure out size of special vars
	if vName in [C_PRESS_I, CF_PRESS_I, ENERGY_I]: #Assume s5 is an rcm.h5 step
		lamData = scRCM.getSpecieslambdata(s5, species='ions')
		Nk = lamData['kEnd']-lamData['kStart']
		return Nk
	elif vName in [C_PRESS_E, CF_PRESS_E, ENERGY_E]: #Assume s5 is an rcm.h5 step
		lamData = scRCM.getSpecieslambdata(s5, species='electrons')
		Nk = lamData['kEnd']-lamData['kStart']
		return Nk
	elif vName == ENERGY_I:
		lamData = scRCM.getSpecieslambdata(s5, species='ions')
		Nk = lamData['kEnd']-lamData['kStart']
		return Nk
	elif vName == ENERGY_E: #Assume s5 is an rcm.h5 step
		lamData = scRCM.getSpecieslambdata(s5, species='electrons')
		Nk = lamData['kEnd']-lamData['kStart']
		return Nk

	if vName not in s5.keys():
		print("Error, '{}' not in s5 keys: {}".format(vName, s5.keys()))
		return 0

	#If here, we are simply returning a variable in the file
	#Its either size 1 if 2d, or var.shape[0] if 2d (rcm file)
	vShape = s5[vName].shape
	if len(vShape) == 2:
		return 1
	elif len(vShape == 3):
		return vShape[0]
	else:
		print("Idk how to get dimension for '{}'".format(vName))
		return 0

def getCumPress(ilamc, eetas, vm):
	Nk = len(ilamc)
	pCum = np.zeros((Nk,3,3))

	ilam_kji = ilamc[:,np.newaxis,np.newaxis]
	vm_kji = vm[np.newaxis, :, :]

	#TODO: Make sure dims work out (3x3 ilamc and vm, kx3x3 eetas)
	pPar = scRCM.pressure_factor*ilam_kji*eetas*vm_kji**2.5 * 1E9  # [Pa -> nPa], partial pressures for each channel

	pCum[0] = pPar[0] # First bin's partial press = cumulative pressure
	for k in range(1,Nk):
		pCum[k] = pCum[k-1] + pPar[k]  # Get current cumulative pressure by adding this bin's partial onto last bin's cumulative
	return pCum

def getVarOnStencil(s5, ciM, ciP, cjM, cjP, vName):
	
	#Simple choice: just return data straight from the file
	if vName in s5.keys():
		vShape = s5[vName].shape
		if len(vShape) == 2:
			return s5[vName][ciM:ciP, cjM:cjP:-1]
		elif len(vShape) == 3:
			return s5[vName][:, ciM:ciP, cjM:cjP:-1]
		else:
			print("Idk how to get stencil for '{}'".format(vName))
			return 0

	#Harder option: we want some derived quantity
	if vName in [C_PRESS_I, C_PRESS_E, CF_PRESS_I, CF_PRESS_E]:  # Can assume s5 = rcm file
		species = 'ions' if vName in [C_PRESS_I, CF_PRESS_I] else 'electrons'
		lamData = scRCM.getSpecieslambdata(s5, species=species)
		kStart = lamData['kStart']
		kEnd = lamData['kEnd']
		ciM += 2; ciP += 2; # To account for rcmGrid offset
		if species == 'ions':
			eetas = np.array(s5['rcmeeta'])[kStart:kEnd+1, ciM:ciP, cjM:cjP:-1]
		else:
			eetas = np.array(s5['rcmeeta'])[kStart:kEnd, ciM:ciP, cjM:cjP:-1]
		vms = np.array(s5['rcmvm'])[ciM:ciP, cjM:cjP:-1]
		cPress = getCumPress(lamData['ilamc'], eetas, vms)

		return cPress if vName in [C_PRESS_I, C_PRESS_E] else cPress/cPress[-1,:,:]

	if vName == ENERGY_I or vName == ENERGY_E:
		species = 'ions' if vName == ENERGY_I else 'electrons'
		lamData = scRCM.getSpecieslambdata(s5, species=species)
		vms = np.array(s5['rcmvm'])[ciM:ciP, cjM:cjP:-1]

		ilamc_kji = lamData['ilamc'][:, np.newaxis, np.newaxis]
		vm_kji = vms[np.newaxis, :, :]
		return ilamc_kji*vm_kji

def getValsAlongSlice(s5, IOpen, lGrid, theta_rad, varList):
	""" 
		s5: hdf5 step object
		IOpen: iopen var from (must have same dims as other values in s5)
		rGrid: 1D slice grid (Assuming its radius in equator)
		theta: theta val of slice
		varList: list of strings indicating values to grab/calc and return
	"""
	if 'xMin' in s5.keys():
		mxstr = 'xMin'
		mystr = 'yMin'
		modelX = np.array(s5[mxstr])
		modelY = np.array(s5[mystr])
	elif 'rcmxmin' in s5.keys():
		mxstr = 'rcmxmin'
		mystr = 'rcmymin'
		modelX = np.array(s5[mxstr][2:,:])
		modelY = np.array(s5[mystr][2:,:])
	Nmi = modelX.shape[0]
	Nmj = modelX.shape[1]
	Ng = len(lGrid)  # number of grid points

	varDict = {}
	for v in varList: #Populate arrays for each variable depending on number of dims needed
		varDict[v] = np.full((Ng, getVarDims(s5,mxstr,v)), FLAG)

	for i_R in range(Ng):
		pntR = lGrid[i_R]
		pntT = (theta_rad*180/np.pi)%360  # [deg], accounting for wrap-around
		pntX = pntR*np.cos(theta_rad)
		pntY = pntR*np.sin(theta_rad)

		#Get closest point
		distSq = (modelX-pntX)**2 + (modelY-pntY)**2
		_i = distSq.argmin()
		ci = int(_i/modelX.shape[1])
		cj = _i%modelX.shape[1]

		#Nvm, very lazily ignore end points
		if ci+1 >= Nmi or ci-1 <= 0:
			continue
		if cj+1 >= Nmj or cj-1 <= 0:
			continue

		if np.any(IOpen[ci-1:ci+2, cj+1:cj-2:-1] > -0.5):
			break

		#Take note of j indexing (i dim in RCM). Must correct so that r values are increasing with increasing j 
		rGrid = np.sqrt(modelX[ci, cj+1:cj-2:-1]**2 + modelY[ci, cj+1:cj-2:-1]**2)
		tGrid = (np.arctan2(modelY[ci-1:ci+2, cj], modelX[ci-1:ci+2, cj])*180/np.pi)%360  # modulo to 'add 360' to negative values
		#vGrid = modelV[ci-1:ci+2, cj+1:cj-2:-1]

		#Now do calculations
		vWeights = interpTSCWeights(tGrid,rGrid,pntT,pntR)
		for v in varList:
			vGrid = getVarOnStencil(s5, ci-1, ci+2, cj+1, cj-2, v)

			val = 0
			for i in range(3):
				for j in range(3):
					#TODO: Handle 2d or 3d vGrid
					if len(vGrid.shape) == 2:
						val += vGrid[i,j]*vWeights[i,j]
					elif len(vGrid.shape) == 3:
						val += vGrid[:,i,j]*vWeights[i,j]
			varDict[v][i_R] = val

	return varDict

def findLoc_cumulFracVal(cfVals, cfData, energies):
	"""Find the energy at the locations of desired values within given cumul. frac. data
		cfVals: 1D array of special values
		cfData: nxk array of cumul. frac. data (n = number of samples, k = number of energies)
		energies: nxk energies
	"""
	shape = (len(cfData), len(cfVals))
	locs = np.full(shape, FLAG)  # Energies at location of desired values
	for s in range(len(cfData)):
		for i in range(len(cfVals)):
			cfv = cfVals[i]
			i_e = 0
			while cfData[s,i_e+1] < cfv: i_e += 1
			#Linear interp between points
			m = (cfData[s,i_e+1] - cfData[s,i_e])/(energies[s,i_e+1]-energies[s,i_e])
			b = cfData[s,i_e] - m*energies[s,i_e]
			loc = (cfv - b)/m
			locs[s,i] = loc

	return locs

def linterp_findLoc(targets, xData, yData, getAxis='y'):
	""" Use linear interpolation to find desired x/y values in data
		targets: 1D array of x/y values to find y/x locations of
		xData/yData: sxk array of cumul. frac. data (s = number of samples, k = number of energies)
		getAxis: desired axis. If 'y', assumes targets are x values, and visa versa
	"""
	shape = (len(xData), len(targets))  # sxt , t = number of target points
	locs = np.full(shape, FLAG)
	for s in range(shape[0]): 
		for i in range(shape[1]):
			#If y data == flag then this sample point should be ignored
			if yData[s][0] == FLAG: continue
			val = targets[i]
			idx = 0
			if getAxis == 'y': # target is x-axis value, find its location in xData
				while xData[s][idx+1] < val: idx += 1
			elif getAxis == 'x': # target is y-axis value, find its location in yData
				while yData[s][idx+1] < val: idx += 1
			m = (yData[s][idx+1] - yData[s][idx])/(xData[s][idx+1] - xData[s][idx])
			b = yData[s][idx] - m*xData[s][idx]
			if getAxis == 'y':
				locs[s][i] = m*val + b
			elif getAxis == 'x':
				locs[s][i] = (val - b)/m
	return locs


if __name__=="__main__":

	fdir  = os.getcwd()
	ftag  = "msphere"
	sStart = -1
	sEnd = -1
	sStride = 10
	jdir = 'jstore'
	pIDs = ['times', 'track', 'vcalc']
	npIDs = ['main', 'tseries']  # NoPlot IDs
	vidOut = "vid_cumulFracKeo"

	vCalcNames = [CF_PRESS_I, CF_PRESS_E]

	parser = argparse.ArgumentParser(description="Plot var from (mhd)rcm over L cut following satellite trajectory")
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of model data (default: %(default)s)")
	parser.add_argument('-jdir',type=str,metavar="directory",default=jdir,help="Directory to store and find json files (default: %(default)s)")
	parser.add_argument('-sStart',type=int, default=sStart,help="Starting time step for L vs. E calculation (default: First step in RCM data)")
	parser.add_argument('-sEnd',type=int, default=sEnd,help="Ending time step for L vs. E calculation (default: Last step in RCM data)")
	parser.add_argument('-sStride',type=int, default=sStride,help="Time step stride for L vs. E calculation (default: %(default)s)")
	parser.add_argument('-forceCalc',type=str,metavar=pIDs,default="",help="Comma-separated process IDs to force recalculation for given process")
	parser.add_argument('-noPlot',type=str,metavar=npIDs,default="",help="Comma-separated process IDs to turn off plotting")
	parser.add_argument('-vidOut',type=str,default=vidOut,help="Output directory (relative to -d) for video images (default: %(default)s)")
	parser.add_argument('-doE',action='store_true')

	args = parser.parse_args()
	fdir  = args.d
	ftag  = args.id
	sStart = args.sStart
	sEnd = args.sEnd
	sStride = args.sStride
	jdir = args.jdir
	fcStr = args.forceCalc
	npStr = args.noPlot
	vidOut = args.vidOut
	doElectrons = args.doE

	fcIDs = [] if (fcStr == "") else fcStr.split(',')
	npIDs = [] if (npStr == "") else npStr.split(',')

	#Init file/directory stuff
	outdir = os.path.join(fdir, vidOut)
	kh5.CheckDirOrMake(outdir)
	mhdrcm_fname = os.path.join(fdir, ftag+'.mhdrcm.h5')
	rcm_fname    = os.path.join(fdir, ftag+'.rcm.h5'   )
	kh5.CheckOrDie(mhdrcm_fname)
	kh5.CheckOrDie(rcm_fname)
	kh5.CheckDirOrMake(jdir)

	#Determine start, stop, stride times and steps
	print('Getting (MHD)RCM times')
	rcmTimes = scRCM.getRCMtimes(rcm_fname,mhdrcm_fname,jdir=jdir,forceCalc=('times' in fcIDs))
	rcmSIDs = rcmTimes['sIDs']
	if sStart == -1:
		sStart = rcmSIDs[0]
	elif sStart < rcmSIDs[0]:
		print("Step '{}' not in RCM times, starting from {}".format(sStart, rcmSIDs[0]))
		sStart = rcmSIDs[0]
	if sEnd == -1:
		sEnd = rcmSIDs[-1]
	elif sEnd > rcmSIDs[-1]:
		print("Step '{}' not in RCM times, ending at {}".format(sEnd, rcmSIDs[-1]))
		sEnd = rcmSIDs[-1]

	iTStart = np.abs(rcmSIDs-sStart).argmin()
	iTEnd = np.abs(rcmSIDs-sEnd).argmin()
	sIDs = rcmSIDs[iTStart:iTEnd+1:sStride]
	sIDstrs = rcmTimes['sIDstrs'][iTStart:iTEnd+1:sStride]
	Nt = len(sIDs)

	print('Getting eq plot info')
	rcm_eqlatlon = scRCM.getRCM_eqlatlon(mhdrcm_fname, rcmTimes, sStart, sEnd, sStride)

	#Some default grids
	#if vName == CF_PRESS_I or vName == CF_PRESS_E:
	rGrid = np.array([3.8,4.8,5.8])
	Ng = len(rGrid)
	theta_rad = 90*np.pi/180
	xPnts = rGrid*np.cos(theta_rad)
	yPnts = rGrid*np.sin(theta_rad)

	#Hard coding varList for now
	if not doElectrons:
		varList = [C_PRESS_I, CF_PRESS_I, ENERGY_I]
		plotVarList = [C_PRESS_I, CF_PRESS_I]
		cfpStr = CF_PRESS_I
		cpStr = C_PRESS_I
		eStr = ENERGY_I
		egBndLow = 1; egBndHigh = 6E2
		pressBndLow = 1E-2; pressBndHigh = 75
	else:
		varList = [C_PRESS_E, CF_PRESS_E, ENERGY_E]
		plotVarList = [C_PRESS_E, CF_PRESS_E]
		cfpStr = CF_PRESS_E
		cpStr = C_PRESS_E
		eStr = ENERGY_E
		egBndLow = .01; egBndHigh=1000 # in keV
		pressBndLow = 0E-8; pressBndHigh = 10
	#eHlVals = np.array([80])*1E3  # energy highlight [eV]
	

	r5 = h5.File(rcm_fname,'r')
	mr5 = h5.File(mhdrcm_fname,'r')

	#Grab all the data we need
	critE_arr = np.full((Nt,len(rGrid)), FLAG)
	critPct_arr = np.full((Nt,len(rGrid)), FLAG)
	vdN = []  # Will hold all the returned data (will be Nt x valDist size)
	ut = kT.MJD2UT(rcm_eqlatlon['MJD'])
	bar = progressbar.ProgressBar(max_value=Nt)
	for n in range(Nt):
		bar.update(n)
		rS5 = r5[sIDstrs[n]]
		mrS5 = mr5[sIDstrs[n]]

		IOpen = mrS5['IOpen'][:]

		valDict = getValsAlongSlice(rS5, IOpen, rGrid, theta_rad, varList)
		vdN.append(valDict)
		#Determine which points were/weren't populated so we can color code dots nicely later
		pntsUsed = np.full(len(rGrid), False)
		for i in range(len(rGrid)):
			pnt = valDict[varList[0]][i]
			if pnt[0] != FLAG:
				pntsUsed[i] = True
	bar.finish()


	#Calc energy at which cumul. frac. is x % of total pressure
	cfHlVals = np.array([0.75, 0.5, 0.25])  # cumul. frac. highlight val(s). Align backwards so plot looks nicer
	NcfHl = len(cfHlVals)
	cfHlPlt  = ['-', '--', '-.']
	cfHlEs   = np.ma.zeros((Ng, NcfHl,  Nt))
	pressureMax = 0
	for p in range(Ng):
		cfData = np.array([vdN[n][cfpStr][p] for n in range(Nt)])
		mask = [cfData[n][-1] < 0.1 for n in range(Nt)]

		energies = [vdN[n][eStr][p] for n in range(Nt)]
		cHlEnergies = linterp_findLoc(cfHlVals, energies, cfData, getAxis='x')
		cfHlEs[p,:,:] = cHlEnergies.T
		for c in range(NcfHl):
			cfHlEs[p,c,:] = np.ma.masked_where(mask, cfHlEs[p,c,:])

		cData = [vdN[n][cpStr][p] for n in range(Nt)]
		pressureMax = np.max((pressureMax, np.max(cData)))

	#Map onto regular grid for the sake of plotting
	eGrid = np.logspace(np.log10(egBndLow), np.log10(egBndHigh), 200, endpoint=True)  # [kev]
	Ne = len(eGrid)
	cpGrid = np.full((Ng,Nt,Ne),np.nan)
	cfpGrid = np.zeros((Ng,Nt,Ne))
	for p in range(Ng):
		for n in range(len(vdN)):
			c = vdN[n][cpStr][p]
			cf = vdN[n][cfpStr][p]
			e = vdN[n][eStr][p]*1E-3
			if cf[-1] < 0.1: continue  # If this value isn't properly populated, we don't want to include it
			for en in range(Ne):
				iE = np.abs(e-eGrid[en]).argmin()
				cpGrid[p,n,en] = c[iE]
				cfpGrid[p,n,en] = cf[iE]
	cpGrid_m = np.ma.masked_where(cfpGrid>0.999, cpGrid)


	isotfmt = '%Y-%m-%d %H:%M:%S'
	hrMarks = ['06:00', '09:15', '12:00', '15:00']
	utMarks = [datetime.datetime.strptime('2013-03-17 {}:00'.format(hr), isotfmt) for hr in hrMarks]

	#------
	#Plotting
	#------

	figMain = plt.figure(figsize=(15,9))
	gs = gridspec.GridSpec(Ng,20, wspace=0.8, hspace=0.2)
	#AxVar1 = figMain.add_subplot(gsMain[:,:])
	AxKContainer = figMain.add_subplot(gs[:,:19])
	AxKeyoList = [figMain.add_subplot(gs[i,:19]) for i in range(Ng)]
	AxCB = figMain.add_subplot(gs[:,19:])

	#pressnorm = kv.genNorm(1E-2, 75, doLog=False)
	#pressnorm = kv.genNorm(0, 20, doLog=False)
	pressnorm = kv.genNorm(pressBndLow, pressBndHigh, doLog=False)
	cmap_press = 'gnuplot'

	for p in range(Ng):
		"""
		for n in range(len(vdN)):
			Ax = AxKeyoList[p]
			c = vdN[n][C_PRESS_I][p]
			cf = vdN[n][CF_PRESS_I][p]
			e = vdN[n][ENERGY_I][p]

			c_m = np.ma.masked_where(cf > 0.999, c)

			#AxVar1.plot(e, c)
			#AxVar1.set_xscale('log')
			#Ax.scatter(np.full(len(e),ut[n]), e*1E-3, c=c_m, cmap='gnuplot')
			Ax.scatter(np.full(len(e),ut[n]), e*1E-3, c=c_m, norm=pressnorm, marker='s',s=7,cmap=cmap_press)
		"""
		Ax = AxKeyoList[p]
		#Ax.pcolor(ut, eGrid, to_center(cpGrid_m[p].T),norm=pressnorm,shading='nearest')
		Ax.pcolor(ut, eGrid,cpGrid_m[p].T,norm=pressnorm,shading='nearest',cmap=cmap_press)
		Ax.set_yscale('log')
		#Ax.set_ylim([5E0,3E2])
		Ax.set_ylim([egBndLow,egBndHigh])

		t = Ax.text(0.99,0.95,"{} $R_E$".format(rGrid[p]), horizontalalignment='right', \
				verticalalignment='top',transform=Ax.transAxes,fontsize=18)
		t.set_bbox(dict(facecolor='white',edgecolor='black',alpha=0.85))

		#Ax.set_ylabel('Energy [keV]')
		Ax.grid(True)

		#Add quartile cf press lines
		for iHlv in range(len(cfHlEs[0,:])):
			lblStr = '{:2.0f}%'.format(cfHlVals[iHlv]*100)
			Ax.plot(ut,cfHlEs[p,iHlv,:]*1E-3,cfHlPlt[iHlv],color='white',label=lblStr)
		if p == 0:
			leg = Ax.legend(loc='upper left',fontsize=14,framealpha=0.95)
			for i in range(len(cfHlVals)):
				leg.legendHandles[i].set_color('black')
		if p < Ng-1:
				Ax.xaxis.set_major_formatter(plt.NullFormatter())


		#Draw UT lines
		for uM in utMarks:
			ymin,ymax = Ax.get_ylim()
			Ax.plot([uM,uM],[ymin,ymax],'--',color='cyan')
			Ax.set_ylim([ymin,ymax])

	kv.genCB(AxCB,pressnorm,r'Cumulative Pressure [nPa]',cM=cmap_press,doVert=True)


	AxKContainer.yaxis.set_major_formatter(plt.NullFormatter())
	AxKContainer.xaxis.set_major_formatter(plt.NullFormatter())
	AxKContainer.tick_params(axis='x',bottom=False,top=False)
	AxKContainer.tick_params(axis='y',left=False,right=False)
	AxKContainer.spines['top'].set_visible(False)
	AxKContainer.spines['right'].set_visible(False)
	AxKContainer.spines['bottom'].set_visible(False)
	AxKContainer.spines['left'].set_visible(False)
	AxKContainer.yaxis.labelpad = 20
	AxKContainer.xaxis.labelpad = 15
	AxKContainer.set_ylabel(r'Energy [keV]',fontsize=18)
	AxKContainer.set_xlabel('UT',fontsize=18)


	filename = 'cumulFrac.png'
	ofname = os.path.join(outdir, filename)
	kv.savePic(ofname)

	#Now make RCM eq plots to show along-side

	figRCM = plt.figure(figsize=(8,8))
	gs = gridspec.GridSpec(7,1, wspace=1.4, hspace=0.6)
	AxRCMEq = figRCM.add_subplot(gs[:6,:])
	AxCB_rcmpress = figRCM.add_subplot(gs[6:,:])

	pressnorm = kv.genNorm(1E-2, 75, doLog=False)
	cmap_press='viridis'

	for iU in range(len(utMarks)):
		uM = utMarks[iU]
		iUT = np.array([np.abs((ut[n]-uM).total_seconds()) for n in range(len(ut))]).argmin()
		pltmjd = rcm_eqlatlon['MJD'][iUT]
		AxRCMEq.clear()
		plt_rcmeq(AxRCMEq, rcm_eqlatlon, xPnts, yPnts, AxCB=AxCB_rcmpress, mjd=pltmjd, norm=pressnorm, cmapName=cmap_press)
		AxRCMEq.title.set_text(ut[iUT].strftime(isotfmt))
		AxRCMEq.title.set_size(18)
		filename = 'eq_{}.png'.format(hrMarks[iU])
		ofname = os.path.join(outdir, filename)
		kv.savePic(ofname)
