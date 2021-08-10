import os, sys
import datetime
import numpy as np
from astropy.time import Time
import datetime

from cdasws import CdasWs
from cdasws import TimeInterval

import kaipy.kaijson as kj

TINY = 1.0e-8

package_directory = os.path.dirname(os.path.abspath(__file__))

scstrs_fname = os.path.join(package_directory, 'sc_cdasws_strs.json')

#======
#General
#======

def mjd_to_ut(mjd_arr):
	isotfmt = '%Y-%m-%dT%H:%M:%S.%f'
	UT = Time(mjd_arr,format='mjd').isot
	return [datetime.datetime.strptime(UT[n],isotfmt) for n in range(len(UT))]

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
			print('  ' + sc)
			for v in scdict[sc].keys():
				print('    ' + v)
	return scdict

def getCdasData(dsName, dsVars, t0, t1, epochStr="Epoch", doVerbose=False):
	"""Pull dataset from CdasWs
		dsName: dataset name
		dsVars: list of variable (strings) from dataset
		t0: start time
		t1: end time
		epochStr: name of Epoch variable in dataset (usually "Epoch", but not always)
	"""

	cdas = CdasWs()

	status,data = cdas.get_data(dsName, dsVars, t0, t1)
	if status['http']['status_code'] != 200:
		# Handle the case where CdasWs just doesn't work if you give it variables in arg 2
		# If given empty var list instead, it'll return the full day on day in t0, and that's it
		# So, call for as many days as we need data for and build one big data object
		if doVerbose: print("Bad pull, trying to build day-by-day")

		t0dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%SZ")
		t1dt = datetime.datetime.strptime(t1, "%Y-%m-%dT%H:%M:%SZ")
		numDays = t1dt.day-t0dt.day + 1 #Number of days we want data from
		if doVerbose: print("numDays: " + str(numDays))

		tstamp_arr = []
		for i in range(numDays):
			tstamp_arr.append((t0dt + datetime.timedelta(days=i)).strftime("%Y-%m-%dT%H:%M:%SZ"))
		
		#Get first day
		status, data = cdas.get_data(dsName, [], tstamp_arr[0], tstamp_arr[1])
		if doVerbose: print("Pulling " + t0)
		
		if status['http']['status_code'] != 200:
			# If it still fails, its some other problem and we'll die
			if doVerbose: print("Still bad pull. Dying.")
			return {}
		if data is None:
			if doVerbose: print("Cdas responded with 200 but returned no data")
			return {}
		
		#Figure out which axes are the epoch axis in each dataset so we can concatenate along it
		nTime = len(data[epochStr])
		dk = list(data.keys())
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
			status, newdata = cdas.get_data(dsName, [], tstamp_arr[i], tstamp_arr[i])
			for k in range(len(dk)):
				if cataxis[k] == -1:
					continue
				else:
					key = dk[k]
					data[key] = np.concatenate((data[key], newdata[key]), axis=cataxis[k])
	else:
		if doVerbose: print("Got data in one pull")

	return data


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


