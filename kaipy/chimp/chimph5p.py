#Various routines to deal with h5p data
import numpy as np
import h5py

#Return time at step n
def tStep(fname,nStp=0):
	with h5py.File(fname,'r') as hf:
		gID = "Step#%d"%(nStp)
		grp = hf.get(gID)
		t = grp.attrs.get("time")
	return t
	
#Count number of timesteps in an h5p file
def cntSteps(fname):
	with h5py.File(fname,'r') as hf:
		grps = hf.values()
		grpNames = [str(grp.name) for grp in grps]
		#Steps = [stp if "/Step#" in stp for stp in grpNames]
		Steps = [stp for stp in grpNames if "/Step#" in stp]
		nSteps = len(Steps)

		return nSteps
#Count number of particles in an h5p file
def cntTPs(fname):
	with h5py.File(fname,'r') as hf:
		grp = hf.get("Step#0")
		ids = (grp.get("id")[()])
		Np = ids.shape[0]
	return Np

def bndTPs(fname):
	with h5py.File(fname,'r') as hf:
		grp = hf.get("Step#0")
		ids = (grp.get("id")[()])
		Np = ids.shape[0]
		nS = ids.min()
		nE = ids.max()
	return Np,nS,nE

#Find array index for a given particle ID (ie if block doesn't start at 1)
def locPID(fname,pid):
	with h5py.File(fname,'r') as hf:
		grp = hf.get("Step#0")
		ids = grp.get("id")[()]
		isP = (ids == pid)
		loc = isP.argmax()
		if (ids[loc] != pid):
			print("Didn't find particle %d ..."%(pid))
			loc = None
			quit()
	return loc
		
#Given an h5part file, create a time series for a single particle w/ ID = pid
def getH5pid(fname,vId,pid):
	p0 = locPID(fname,pid) #Find particle in array
	Nt = cntSteps(fname) #Find number of slices

	V = np.zeros(Nt)
	t = np.zeros(Nt)
	with h5py.File(fname,'r') as hf:
		for n in range(Nt):
			#Create gId
			gId = "Step#%d"%(n)
			grp = hf.get(gId)
			t[n] = grp.attrs.get("time")
			V[n] = (grp.get(vId)[()])[p0]
	return t,V

#Given an h5part file, create a time series from an input string
def getH5p(fname,vId,Mask=None):
	Nt = cntSteps(fname)
	Np = cntTPs(fname)
	t = np.zeros(Nt)
	V = np.zeros((Nt,Np))	
	with h5py.File(fname,'r') as hf:
		for n in range(Nt):
			#Create gId
			gId = "Step#%d"%(n)
			grp = hf.get(gId)
			t[n] = grp.attrs.get("time")
			V[n,:] = grp.get(vId)[()]
	if (Mask is None):
		return t,V

	else:
		return t,V[:,Mask]
#Given an h5p file and step, get one slice of data
def getH5pT(fname,vID="isIn",nStp=0,cutIn=False):
	with h5py.File(fname,'r') as hf:
		gID = "Step#%d"%(nStp)
		V = hf.get(gID).get(vID)[()]

	return V
	