import h5py
import numpy as np

def cntSteps(fname):
		with h5py.File(fname,'r') as hf:
				grps = hf.values()
				grpNames = [str(grp.name) for grp in grps]
				#Steps = [stp if "/Step#" in stp for stp in grpNames]
				Steps = [stp for stp in grpNames if "/Step#" in stp]
				nSteps = len(Steps)
				
				sIds = np.array([str.split(s,"#")[-1] for s in Steps],dtype=np.int)
				return nSteps,sIds

def getTs(fname,sIds):
	Nt = len(sIds)
	T = np.zeros(Nt)
	i0 = sIds.min()
	i1 = sIds.max()
	with h5py.File(fname,'r') as hf:
		for n in range(i0,i1+1):
			gId = "/Step#%d"%(n)
			if ("time" in hf[gId].attrs.keys()):
				T[n-i0] = hf[gId].attrs["time"]
			else:
				T[n-i0] = -np.inf
	return T

#Get shape/dimension of grid
def getDims(fname,doFlip=True):
	with h5py.File(fname,'r') as hf:
		Dims = hf["/"]["X"].shape
	Ds = np.array(Dims,dtype=np.int)
	if (doFlip):
		Ds = np.flip(Ds,axis=0)
	return Ds

#Get root variables
def getRootVars(fname):
	with h5py.File(fname,'r') as hf:
		vIds = []
		for k in hf.keys():
			#Don't include stuff that starts with step
			if "Step" not in str(k):
				vIds.append(str(k))
	#Remove coordinates from list of root variables
	xyzS = ["X","Y","Z"]
	for s in xyzS:
		if s in vIds:
			vIds.remove(s)
	return vIds

#Get variables in initial Step
def getVars(fname,s0):
	with h5py.File(fname,'r') as hf:
		gId = "/Step#%d"%(s0)
		stp0 = hf[gId]
		vIds = []
		for k in stp0.keys():
			L = stp0[k].shape[0]
			vIds.append(str(k))
	return vIds

#Get variable data
def PullVar(fname,vID,s0=None):
	with h5py.File(fname,'r') as hf:
		if (s0 is None):
			V = hf[vID].value.T
		else:
			gId = "/Step#%d"%(s0)
			V = hf[gId][vID].value.T
			#V = hf[gId][vID][:].T
	return V
