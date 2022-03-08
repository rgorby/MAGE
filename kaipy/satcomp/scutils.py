import numpy as np
import datetime
import os
import glob
import sys
import subprocess
from xml.dom import minidom

from astropy.time import Time


from cdasws import CdasWs
from cdasws import TimeInterval

import kaipy.kaijson as kj
import kaipy.kaiTools as kaiTools

import spacepy
from spacepy.coordinates import Coords
from spacepy.time import Ticktock
import spacepy.datamodel as dm

#import hdf5
import h5py

TINY = 1.0e-8

package_directory = os.path.dirname(os.path.abspath(__file__))

scstrs_fname = os.path.join(package_directory, 'sc_cdasws_strs.json')

#======
#General
#======

def trilinterp(xbnd, ybnd, zbnd, valbnd, x, y, z):
	"""3D linear interpolation
		xbnd,ybnd,zbnd: 2-element arrays each of the bounding dimensions
		valbnd: 2x2x2 array of variable to be interpolated
		x, y, z: point inside bounds
	"""

	xd = (x - xbnd[0])/(xbnd[1]-xbnd[0])
	yd = (y - ybnd[0])/(ybnd[1]-ybnd[0])
	zd = (z - zbnd[0])/(zbnd[1]-zbnd[0])
	v00 = valbnd[0,0,0]*(1-xd) + valbnd[1,0,0]*xd
	v01 = valbnd[0,0,1]*(1-xd) + valbnd[1,0,1]*xd
	v10 = valbnd[0,1,0]*(1-xd) + valbnd[1,1,0]*xd
	v11 = valbnd[0,1,1]*(1-xd) + valbnd[1,1,1]*xd
	v0 = v00*(1-yd) + v10*yd
	v1 = v01*(1-yd) + v11*yd
	v = v0*(1-zd) + v1*zd
	
	return v

def varMap_1D(og, ng, var):
	"""Map variable from one grid to another
	og: old grid
	ng: new grid
	var: variable to re-map
	"""
	varnew = np.zeros((len(ng)))

	for e in range(len(ng)):
		if ng[e] < og[0] or ng[e] > og[-1]:
			continue

		idx = 0
		while og[idx+1] < ng[e]: idx += 1

		glow = og[idx]
		ghigh = og[idx+1]
		d = (ng[e] - glow)/(ghigh-glow)
		

		varnew[e] = var[idx]*(1-d) + var[idx+1]*d
	return varnew

def getWeights_ConsArea(og, og_lower, og_upper, ng, ng_lower, ng_upper):
	"""Calculate overlap (weights) to map values on one grid to another,
		where total width in grid dimension are conserved
		(i.e. properly map RCM eetas to uniform grid)
		og, og_lower, og_upper: old grid and lower/upper bounds of each grid point
		ng, ng_lower, ng_upper: new grid and lower/upper bounds of each grid point
	"""
	"""Example:
		og: || |  |  |   |   |     |     |
		ng: |  |  |  |  |  |  |  |  |  |  |
		For each cell center on ng, calc which og cells overlap and the fraction of overlap
	"""

	Nog = len(og)
	Nng = len(ng)
	weightMap = [[] for e in range(Nng) ]  # Ne x (nx2)
	for iNG in range(Nng):
		ng_l = ng_lower[iNG]
		ng_u = ng_upper[iNG]
		ng_w = ng_u - ng_l  # cell width
		frac_arr = []
		for k in range(Nog):
			#Do these two cells overlap
			if ng_l <= og_upper[k] and ng_u >= og_lower[k]:
				#Get overlap bounds and width
				ovl_lower = og_lower[k] if og_lower[k] > ng_l else ng_l
				ovl_upper = og_upper[k] if og_upper[k] < ng_u else ng_u
				ovl_width = ovl_upper - ovl_lower
				#og_width = og_upper[k] - og_lower[k]
				frac_arr.append([k, ovl_width/ng_w])
		weightMap[iNG] = frac_arr
	return weightMap

def computeErrors(obs,pred):
	MAE = 1./len(obs) *np.sum(np.abs(obs-pred))
	MSE = 1./len(obs) * np.sum((obs-pred)**2)
	RMSE = np.sqrt(MSE)
	MAPE = 1./len(obs) * np.sum(np.abs(obs-pred)/
		np.where(abs(obs) < TINY,TINY,abs(obs)))
	RSE = (np.sum((obs-pred)**2)/
		np.where(np.sum((obs-np.mean(obs))**2)<TINY,TINY,np.sum((obs-np.mean(obs))**2)))
	PE = 1-RSE
	return MAE,MSE,RMSE,MAPE,RSE,PE

#======
#Cdaweb-related
#======

def getScIds(doPrint=False):
	"""Load info from stored file containing strings needed to get certain spacefract datasets from cdaweb
	"""
	scdict = kj.load(scstrs_fname)

	if doPrint:
		print("Retrievable spacecraft data:")
		for sc in scdict.keys():
			print('	 ' + sc)
			for v in scdict[sc].keys():
				print('	   ' + v)
	return scdict

def getCdasDsetInterval(dsName):
	cdas = CdasWs()

	data = cdas.get_datasets(idPattern=dsName)
	if len(data) == 0:
		return None, None
	tInt = data[0]['TimeInterval']
	return tInt['Start'], tInt['End']

def pullVar(cdaObsId,cdaDataId,t0,t1,deltaT=60,epochStr="Epoch",doVerbose=False):
	"""Pulls info from cdaweb
		cdaObsId  : [str] Dataset name
		cdaDataId : [str or list of strs] variables from dataset
		t0	  : [str] start time, formatted as '%Y-%m-%dT%H:%M:%S.%f'
		t1	  : [str] end time, formatted as '%Y-%m-%dT%H:%M:%S.%f'
		deltaT	  : [float] time cadence [sec], used when interping through time with no data
		epochStr  : [str] name of Epoch var in dataset. Used when needing to build day-by-day
		doVerbose : [bool] Helpful for debugging/diagnostics
	"""

	binData={'interval' : deltaT, 
			 'interpolateMissingValues' : True,
			 'sigmaMultipler' : 4}

	cdas = CdasWs()
	status,data =  cdas.get_data(cdaObsId,cdaDataId,t0,t1,binData=binData)

	if status['http']['status_code'] != 200 or data is None:
		# Handle the case where CdasWs just doesn't work if you give it variables in arg 2
		# If given empty var list instead, it'll return the full day on day in t0, and that's it
		# So, call for as many days as we need data for and build one big data object
		if doVerbose: print("Bad pull, trying to build day-by-day")

		if '.' in t0:
			t0dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%S.%fZ")
			t1dt = datetime.datetime.strptime(t1, "%Y-%m-%dT%H:%M:%S.%fZ")
		else:
			t0dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%SZ")
			t1dt = datetime.datetime.strptime(t1, "%Y-%m-%dT%H:%M:%SZ")
		numDays = t1dt.day-t0dt.day + 1 #Number of days we want data from
		if doVerbose: print("numDays: " + str(numDays))

		tstamp_arr = []
		tstamp_deltas = []
		for i in range(numDays):
			tstamp_arr.append((t0dt + datetime.timedelta(days=i)).strftime("%Y-%m-%dT%H:%M:%SZ"))
			tstamp_deltas.append((t0dt + datetime.timedelta(days=i+1)).strftime("%Y-%m-%dT%H:%M:%SZ"))
		if doVerbose: print("Tstamp_arr: " + str(tstamp_arr))
		#Get first day
		status, data = cdas.get_data(cdaObsId, [], tstamp_arr[0], tstamp_deltas[0], binData=binData)
		if doVerbose: print("Pulling " + t0)
		
		if status['http']['status_code'] != 200:
			# If it still fails, its some other problem and we'll die
			if doVerbose: print("Still bad pull. Dying.")
			return status,data
		if data is None:
			if doVerbose: print("Cdas responded with 200 but returned no data")
			return status,data
		if epochStr not in data.keys():
			if doVerbose: print(epochStr + " not in dataset, can't build day-by-day")
			data = None
			return status,data
		
		#Figure out which axes are the epoch axis in each dataset so we can concatenate along it
		dk = list(data.keys())
		nTime = len(data[epochStr])
		cataxis = np.array([-1 for i in range(len(dk))])
		for k in range(len(dk)):
			shape = np.array(data[dk[k]].shape)
			for i in range(len(shape)):
				if shape[i] == nTime:
					cataxis[k] = i
					continue

		#Then append rest of data accordingly
		for i in range(1,numDays):
			if doVerbose: print("Pulling " + str(tstamp_arr[i]))
			status, newdata = cdas.get_data(cdaObsId, [], tstamp_arr[i], tstamp_deltas[i], binData=binData)
			for k in range(len(dk)):
				if cataxis[k] == -1:
					continue
				else:
					key = dk[k]
					data[key] = np.concatenate((data[key], newdata[key]), axis=cataxis[k])
	else:
		if doVerbose: print("Got data in one pull")


	return status,data

def addVar(mydata,scDic,varname,t0,t1,deltaT,epochStr='Epoch'):
	#print(scDic,varname,idname,dataname,scDic[idname])
	if scDic[varname]['Id'] is not None:
		status,data = pullVar(scDic[varname]['Id'],scDic[varname]['Data'],t0,t1,deltaT,epochStr=epochStr)
		#print(status)
		if status['http']['status_code'] == 200 and data is not None:
			mydata[varname] = dm.dmarray(data[scDic[varname]['Data']],
										 attrs=data[scDic[varname]['Data']].attrs)
			#mydata.tree(attrs=True)
	else:
		#Mimic the cdasws return code for case when id isn't provided
		status = {'http': {'status_code': 404}}
	return status

def getSatData(scDic,t0,t1,deltaT):
	#First get the empheris data if it doesn't exist return the failed status code and
	#go no further
	status,data = pullVar(scDic['Ephem']['Id'],scDic['Ephem']['Data'],
						  t0,t1,deltaT)
	if status['http']['status_code'] != 200 or data is None:
		print('Unable to get data for ', scDic['Ephem']['Id'])
		return status,data
	else:
		#data.tree(attrs=True)
		mydata = dm.SpaceData(attrs={'Satellite':data.attrs['Source_name']})
		if 'Epoch_bin' in data.keys():
			#print('Using Epoch_bin')
			mytime = data['Epoch_bin']
			epochStr = 'Epoch_bin'
		elif 'Epoch' in data.keys():
			#print('Using Epoch')
			mytime = data['Epoch']
			epochStr = 'Epoch'
		elif ([key for key in data.keys() if key.endswith('_state_epoch')]):
			epochStr = [key for key in data.keys() if key.endswith('_state_epoch')][0]
			#mytime = data[[key for key in data.keys()
			#if key.endswith('_state_epoch')][0]]
			mytime = data[epochStr]
		else:
			print('Unable to determine time type')
			status = {'http': {'status_code': 404}}
			return status,data
		mydata['Epoch_bin'] = dm.dmarray(mytime,
										 attrs=mytime.attrs)
		mydata['Ephemeris'] = dm.dmarray(data[scDic['Ephem']['Data']],
										 attrs= data[scDic['Ephem']['Data']].attrs)
		keys = ['MagneticField','Velocity','Density','Pressure']
		for key in keys:
			if key in scDic:
				status1 = addVar(mydata,scDic,key,t0,t1,deltaT,epochStr=epochStr)

		#Add any metavar since they might be needed for unit/label determination
		search_key = 'metavar'
		res = [key for key,val in data.items() if search_key in key]
		for name in res:
			try:
				len(mydata[name])
			except:
				mydata[name] = dm.dmarray([data[name]],attrs=data[name].attrs)
			else:
				mydata[name] = dm.dmarray(data[name],attrs=data[name].attrs)

	return status,mydata


#======
#Shared data derivations
#======

def xyz_to_L(x, y, z):
	r = np.sqrt(x**2 + y**2 + z**2)
	lat = np.arctan2(z, np.sqrt(x**2 + y**2))
	#Convert sc location to L shell, assuming perfect dipole
	lat = lat*np.pi/180.0  # deg to rad
	return r/np.cos(lat)**2

def getJScl(Bmag,Beq,en=2.0):
	#Given sin^n(alpha) dep. on intensity calculate fraction based on accessible Alpha

	Na = 360
	A = np.linspace(0,0.5*np.pi,Na)
	da = A[1]-A[0]
	Ia = np.sin(A)**en
	Ic = np.zeros(Ia.shape)
	Nt = len(Bmag)
	I0 = Ia.sum()

	It = np.zeros(Nt)
	for n in range(Nt):
		if (Bmag[n]<TINY):
			It[n] = 0.0
		else:
			Ac = np.arcsin(np.sqrt(Beq[n]/Bmag[n]))
			Ic[:] = Ia[:]
			Icut = (A>Ac)
			Ic[Icut] = 0.0
			It[n] = Ic.sum()/I0
	return It

def genSCXML(fdir,ftag,
	scid="sctrack_A",h5traj="sctrack_A.h5",numSegments=1):

	(fname,isMPI,Ri,Rj,Rk) = kaiTools.getRunInfo(fdir,ftag)
	root = minidom.Document()
	xml = root.createElement('Kaiju')
	root.appendChild(xml)
	chimpChild = root.createElement('Chimp')
	scChild = root.createElement("sim")
	scChild.setAttribute("runid",scid)
	chimpChild.appendChild(scChild)
	fieldsChild = root.createElement("fields")
	fieldsChild.setAttribute("doMHD","T")
	fieldsChild.setAttribute("grType","LFM")
	fieldsChild.setAttribute("ebfile",ftag)
	if isMPI:
		fieldsChild.setAttribute("isMPI","T")
	chimpChild.appendChild(fieldsChild)
	if isMPI:
		parallelChild = root.createElement("parallel")
		parallelChild.setAttribute("Ri","%d"%Ri)
		parallelChild.setAttribute("Rj","%d"%Rj)
		parallelChild.setAttribute("Rk","%d"%Rk)
		chimpChild.appendChild(parallelChild)
	unitsChild = root.createElement("units")
	unitsChild.setAttribute("uid","EARTH")
	chimpChild.appendChild(unitsChild)
	trajChild = root.createElement("trajectory")
	trajChild.setAttribute("H5Traj",h5traj)
	trajChild.setAttribute("doSmooth","F")
	chimpChild.appendChild(trajChild)
	if numSegments > 1:
		parInTimeChild = root.createElement("parintime")
		parInTimeChild.setAttribute("NumB","%d"%numSegments)
		chimpChild.appendChild(parInTimeChild)
	xml.appendChild(chimpChild)
	return root


#======
#SCTrack
#======

def convertGameraVec(x,y,z,ut,fromSys,fromType,toSys,toType):
	invec = Coords(np.column_stack((x,y,z)),fromSys,fromType)
	invec.ticks = Ticktock(ut)
	outvec = invec.convert(toSys,toType)
	return outvec

def createInputFiles(data,scDic,scId,mjd0,sec0,fdir,ftag,numSegments):
	Re = 6380.0
	toRe = 1.0
	if 'UNITS' in data['Ephemeris'].attrs:
		if "km" in data['Ephemeris'].attrs['UNITS']:
			toRe = 1.0/Re
	elif 'UNIT_PTR' in data['Ephemeris'].attrs:
		if data[data['Ephemeris'].attrs['UNIT_PTR']][0]:
			toRe = 1.0/Re
	if 'SM' == scDic['Ephem']['CoordSys']:
		smpos = Coords(data['Ephemeris'][:,0:3]*toRe,'SM','car')
		smpos.ticks = Ticktock(data['Epoch_bin'])
	elif 'GSM' == scDic['Ephem']['CoordSys'] :
		scpos = Coords(data['Ephemeris'][:,0:3]*toRe,'GSM','car')
		scpos.ticks = Ticktock(data['Epoch_bin'])
		smpos = scpos.convert('GSE','car')
		scpos = Coords(data['Ephemeris'][:,0:3]*toRe,'GSM','car')
		scpos.ticks = Ticktock(data['Epoch_bin'])
		smpos = scpos.convert('SM','car')
	elif 'GSE'== scDic['Ephem']['CoordSys']:
		scpos = Coords(data['Ephemeris'][:,0:3]*toRe,'GSE','car')
		scpos.ticks = Ticktock(data['Epoch_bin'])
		smpos = scpos.convert('SM','car')
	else:
		print('Coordinate system transformation failed')
		return
	elapsedSecs = (smpos.ticks.getMJD()-mjd0)*86400.0+sec0
	scTrackName = os.path.join(fdir,scId+".sc.h5")
	with h5py.File(scTrackName,'w') as hf:
		hf.create_dataset("T" ,data=elapsedSecs)
		hf.create_dataset("X" ,data=smpos.x)
		hf.create_dataset("Y" ,data=smpos.y)
		hf.create_dataset("Z" ,data=smpos.z)
	chimpxml = genSCXML(fdir,ftag,
		scid=scId,h5traj=os.path.basename(scTrackName),numSegments=numSegments)
	xmlFileName = os.path.join(fdir,scId+'.xml')
	with open(xmlFileName,"w") as f:
		f.write(chimpxml.toprettyxml())

	return (scTrackName,xmlFileName)

def addGAMERA(data,scDic,h5name):
	h5file = h5py.File(h5name, 'r')
	ut = kaiTools.MJD2UT(h5file['MJDs'][:])

	bx = h5file['Bx']
	by = h5file['By']
	bz = h5file['Bz']

	if not 'MagneticField' in scDic:
		toCoordSys = 'GSM'
	else:
		toCoordSys = scDic['MagneticField']['CoordSys']
	lfmb_out = convertGameraVec(bx[:],by[:],bz[:],ut,
		'SM','car',toCoordSys,'car')
	data['GAMERA_MagneticField'] = dm.dmarray(lfmb_out.data,
		attrs={'UNITS':bx.attrs['Units'],
		'CATDESC':'Magnetic Field, cartesian'+toCoordSys,
		'FIELDNAM':"Magnetic field",'AXISLABEL':'B'})
	vx = h5file['Vx']
	vy = h5file['Vy']
	vz = h5file['Vz']
	if not 'Velocity' in scDic:
		toCoordSys = 'GSM'
	else:
		toCoordSys = scDic['Velocity']['CoordSys']
	lfmv_out = convertGameraVec(vx[:],vy[:],vz[:],ut,
		'SM','car',toCoordSys,'car')
	data['GAMERA_Velocity'] = dm.dmarray(lfmv_out.data,
		attrs={'UNITS':vx.attrs['Units'],
		'CATDESC':'Velocity, cartesian'+toCoordSys,
			   'FIELDNAM':"Velocity",'AXISLABEL':'V'})
	den = h5file['D']
	data['GAMERA_Density'] = dm.dmarray(den[:],
		attrs={'UNITS':den.attrs['Units'],
		'CATDESC':'Density','FIELDNAM':"Density",'AXISLABEL':'n'})
	pres = h5file['P']
	data['GAMERA_Pressure'] = dm.dmarray(pres[:],
		attrs={'UNITS':pres.attrs['Units'],
		'CATDESC':'Pressure','FIELDNAM':"Pressure",'AXISLABEL':'P'})
	temp = h5file['T']
	data['GAMERA_Temperature'] = dm.dmarray(temp[:],
		attrs={'UNITS':pres.attrs['Units'],
		'CATDESC':'Temperature','FIELDNAM':"Temperature",
		'AXISLABEL':'T'})
	inDom = h5file['inDom']
	data['GAMERA_inDom'] = dm.dmarray(inDom[:],
		attrs={'UNITS':inDom.attrs['Units'],
		'CATDESC':'In GAMERA Domain','FIELDNAM':"InDom",
		'AXISLABEL':'In Domain'})
	return

def matchUnits(data):
	vars = ['Density','Pressure','Temperature','Velocity','MagneticField']
	for var in vars:
		try:
			data[var]
		except:
			print(var,'not in data')
		else:
			if (data[var].attrs['UNITS'] == data['GAMERA_'+var].attrs['UNITS'].decode()):
				print(var,'units match')
			else:
				if 'Density' == var:
					if (data[var].attrs['UNITS'] == 'cm^-3' or data[var].attrs['UNITS'] == '/cc'):
						data[var].attrs['UNITS'] = data['GAMERA_'+var].attrs['UNITS']
						print(var,'units match')
					else:
						print('WARNING ',var,'units do not match')
				if 'Velocity' == var:
					if (data[var].attrs['UNITS'] == 'km/sec'):
						data[var].attrs['UNITS'] = data['GAMERA_'+var].attrs['UNITS']
						print(var,'units match')
					else:
						print('WARNING ',var,'units do not match')
				if 'MagneticField' == var:
					if (data[var].attrs['UNITS'] == '0.1nT'):
						print('Magnetic Field converted from 0.1nT to nT')
						data[var]=data[var]/10.0
						data[var].attrs['UNITS'] = 'nT'
					else:
						print('WARNING ',var,'units do not match')
				if 'Pressure' == var:
					print('WARNING ',var,'units do not match')
				if 'Temperature' == var:
					print('WARNING ',var,'units do not match')

	return

def extractGAMERA(data,scDic,scId,mjd0,sec0,fdir,ftag,cmd,numSegments,keep):

	(scTrackName,xmlFileName) = createInputFiles(data,scDic,scId,
		mjd0,sec0,fdir,ftag,numSegments)

	if 1 == numSegments:
		sctrack = subprocess.run([cmd, xmlFileName], cwd=fdir,
							stdout=subprocess.PIPE, stderr=subprocess.PIPE,
							text=True)
		#print(sctrack)
		h5name = os.path.join(fdir, scId + '.sc.h5')

	else:
		process = []
		for seg in range(1,numSegments+1):
			process.append(subprocess.Popen([cmd, xmlFileName,str(seg)],
							cwd=fdir,
							stdout=subprocess.PIPE, stderr=subprocess.PIPE,
							text=True))
		for proc in process:
			proc.communicate()
		h5name = mergeFiles(scId,fdir,numSegments)


	addGAMERA(data,scDic,h5name)

	if not keep:
		subprocess.run(['rm',h5name])
		subprocess.run(['rm',xmlFileName])
		subprocess.run(['rm',scTrackName])
		if numSegments > 1:
			h5parts = os.path.join(fdir,scId+'.*.sc.h5')
			subprocess.run(['rm',h5parts])
	return

def copy_attributes(in_object, out_object):
	'''Copy attributes between 2 HDF5 objects.'''
	for key, value in list(in_object.attrs.items()):
		out_object.attrs[key] = value


def createMergeFile(fIn,fOut):
	iH5 = h5py.File(fIn,'r')
	oH5 = h5py.File(fOut,'w')
	copy_attributes(iH5,oH5)
	for Q in iH5.keys():
		oH5.create_dataset(Q,data=iH5[Q],maxshape=(None,))
		copy_attributes(iH5[Q],oH5[Q])
	iH5.close()
	return oH5

def addFileToMerge(mergeH5,nextH5):
	nS = nextH5.attrs['nS']
	nE = nextH5.attrs['nE']
	for varname in mergeH5.keys():
		dset = mergeH5[varname]
		dset.resize(dset.shape[0]+nextH5[varname].shape[0],axis=0)
		dset[nS-1:nE]=nextH5[varname][:]
	return

def mergeFiles(scId,fdir,numSegments):
	seg = 1
	inH5Name = os.path.join(fdir,scId+'.%04d'%seg+'.sc.h5')
	mergeH5Name = os.path.join(fdir,scId+'.sc.h5')
	mergeH5 = createMergeFile(inH5Name,mergeH5Name)
	#print(inH5Name,mergeH5Name)
	for seg in range(2,numSegments+1):
		nextH5Name = os.path.join(fdir,scId+'.%04d'%seg+'.sc.h5')
		nextH5 = h5py.File(nextH5Name,'r')
		addFileToMerge(mergeH5,nextH5)

	return mergeH5Name

def genSatCompPbsScript(scId,fdir,cmd,account='P28100045'):
	headerString = """#!/bin/tcsh
#PBS -A %s
#PBS -N %s
#PBS -j oe
#PBS -q casper
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1
"""
	moduleString = """module purge
module load git/2.22.0 intel/18.0.5 hdf5/1.10.5 impi/2018.4.274
module load ncarenv/1.3 ncarcompilers/0.5.0 python/3.7.9 cmake/3.14.4
module load ffmpeg/4.1.3 paraview/5.8.1 mkl/2018.0.5
ncar_pylib /glade/p/hao/msphere/gamshare/casper_satcomp_pylib
module list
"""
	commandString = """cd %s
setenv JNUM ${PBS_ARRAY_INDEX}
date
echo 'Running analysis'
%s %s $JNUM
date
"""
	xmlFileName = os.path.join(fdir,scId+'.xml')
	pbsFileName = os.path.join(fdir,scId+'.pbs')
	pbsFile = open(pbsFileName,'w')
	pbsFile.write(headerString%(account,scId))
	pbsFile.write(moduleString)
	pbsFile.write(commandString%(fdir,cmd,xmlFileName))
	pbsFile.close()

	return pbsFileName

def genSatCompLockScript(scId,fdir,account='P28100045'):
	headerString = """#!/bin/tcsh
#PBS -A %s
#PBS -N %s
#PBS -j oe
#PBS -q casper
#PBS -l walltime=0:15:00
#PBS -l select=1:ncpus=1
"""
	commandString = """cd %s
touch %s
"""
	pbsFileName = os.path.join(fdir,scId+'.done.pbs')
	pbsFile = open(pbsFileName,'w')
	pbsFile.write(headerString%(account,scId))
	pbsFile.write(commandString%(fdir,scId+'.lock'))
	pbsFile.close()

	return pbsFileName

def errorReport(errorName,scId,data):

	keysToCompute = []
	keys = data.keys()

	print('Writing Error to ',errorName)
	f = open(errorName,'w')
	if 'Density' in keys:
		keysToCompute.append('Density')
	if 'Pressue' in keys:
		keysToCompute.append('Pressue')
	if 'Temperature' in keys:
		keysToCompute.append('Temperature')
	if 'MagneticField' in keys:
		keysToCompute.append('MagneticField')
	if 'Velocity' in keys:
		keysToCompute.append('Velocity')

	for key in keysToCompute:
		#print('Plotting',key)
		if 'MagneticField' == key or 'Velocity' == key:
			for vecComp in range(3):
				maskedData = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
					data[key][:,vecComp])
				maskedGamera = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
					data['GAMERA_'+key][:,vecComp])
				MAE,MSE,RMSE,MAPE,RSE,PE = computeErrors(maskedData,maskedGamera)
				f.write(f'Errors for: {key},{vecComp}\n')
				f.write(f'MAE: {MAE}\n')
				f.write(f'MSE: {MSE}\n')
				f.write(f'RMSE: {RMSE}\n')
				f.write(f'MAPE: {MAPE}\n')
				f.write(f'RSE: {RSE}\n')
				f.write(f'PE: {PE}\n')
		else:
				maskedData = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
					data[key][:])
				maskedGamera = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
					data['GAMERA_'+key][:])
				MAE,MSE,RMSE,MAPE,RSE,PE = computeErrors(maskedData,maskedGamera)
				f.write(f'Errors for: {key}\n')
				f.write(f'MAE: {MAE}\n')
				f.write(f'MSE: {MSE}\n')
				f.write(f'RMSE: {RMSE}\n')
				f.write(f'MAPE: {MAPE}\n')
				f.write(f'RSE: {RSE}\n')
				f.write(f'PE: {PE}\n')
	f.close()
	return









