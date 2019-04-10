#Various routines to make 3D particle trajectories

import numpy as np
import h5py

#Nt = # of time slices
#Np = # of particles
#Nv = # of variables
#xyz = |Nt,Np,3|
#V = |Nt,Np,Nv|
#vIDs = |Nv|
def AddTPs(xyz,V=None,vIDs=None,fOut="tps.h5",nStp=0):
	gID = "Step#%d"%(nStp)
	oH5 = h5py.File(fOut,'w')
	Nt,Np,Nd = xyz.shape
	if (vIDs is not None):
		Nv = len(vIDs)
	else:
		Nv = 0

	#Create step data
	oH5.create_group(gID)
	oH5[gID].attrs.create("time",nStp)

	for m in range(Np):
		lID = "Line#%d"%(m)

		#Add connectivity
		ijL = np.zeros((Nt-1,2),dtype=np.int)
		ijL[:,0] = np.arange(Nt-1)
		ijL[:,1] = np.arange(1,Nt)
		#Save line data
		oH5[gID].create_group(lID)
		oH5[gID][lID].attrs.create("Np",Nt)
		oH5[gID][lID].create_dataset("xyz",data=np.single(xyz[:,m,:]))
		oH5[gID][lID].create_dataset("LCon",data=ijL)

		for i in range(Nv):
			oH5[gID][lID].create_dataset(vIDs[i],data=np.single(V[:,m,i]))

	oH5.close()