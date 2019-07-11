#Various routines to deal with LFM-style data
import numpy as np
import kaipy.gamera.gamGrids as gg

clight = 2.9979e+10 #Speed of light [cm/s]
Mp = 1.6726219e-24 #g
gamma = (5.0/3)

def lfm2gg(fIn,fOut,doEarth=True,doJupiter=False):
	from pyhdf.SD import SD, SDC
	import h5py

	#Choose scaling
	if (doEarth):
		xScl = gg.Re
		print("Using Earth scaling ...")
	elif (doJupiter):
		xScl = gg.Rj
		print("Using Jovian scaling ...")
	iScl = 1/xScl
	hdffile = SD(fIn)
	#Grab x/y/z arrays from HDF file.  Scale by Re/Rj/Rs
	#LFM is k,j,i ordering
	x3 = iScl*np.double(hdffile.select('X_grid').get())
	y3 = iScl*np.double(hdffile.select('Y_grid').get())
	z3 = iScl*np.double(hdffile.select('Z_grid').get())
	lfmNc = x3.shape #Number of corners (k,j,i)
	nk = x3.shape[0]-1
	nj = x3.shape[1]-1
	ni = x3.shape[2]-1

	print("Reading LFM grid from %s, size (%d,%d,%d)"%(fIn,ni,nj,nk))

	with h5py.File(fOut,'w') as hf:
		hf.create_dataset("X",data=x3)
		hf.create_dataset("Y",data=y3)
		hf.create_dataset("Z",data=z3)

#Get LFM times
def lfmTimes(hdfs):
	from pyhdf.SD import SD, SDC
	Ts = [ SD(fIn).attributes().get('time') for fIn in hdfs]

	return np.array(Ts)

#Get LFM fields
def lfmFields(fIn):
	from pyhdf.SD import SD, SDC

	hdffile = SD(fIn)
	#Get cell-centered fields
	Bx3cc,By3cc,Bz3cc = getHDFVec(hdffile,'b')
	Vx3cc,Vy3cc,Vz3cc = getHDFVec(hdffile,'v')

	# #Get E field, E= - (VxB)/c
	# #Do this *BEFORE* subtracting off Earth dipole
	# scl = -1.0/clight
	# Ex3cc = scl*( Vy3cc*Bz3cc - Vz3cc*By3cc)
	# Ey3cc = scl*( Vz3cc*Bx3cc - Vx3cc*Bz3cc)
	# Ez3cc = scl*( Vx3cc*By3cc - Vy3cc*Bx3cc)


	return Vx3cc,Vy3cc,Vz3cc,Bx3cc,By3cc,Bz3cc

#Get LFM MHD variables
#Convert (D/Cs) -> (n,P)
#Returns units (#/cm3) and (nPa)

def lfmFlow(fIn):
	from pyhdf.SD import SD, SDC

	hdffile = SD(fIn)

	#Get soundspeed [km/s]
	C3 = getHDFScl(hdffile,"c",Scl=1.0e-5)
	#Get rho [g/cm3]
	D3 = getHDFScl(hdffile,"rho")

	#Conversion to MKS for P in Pascals
	D_mks = (D3*1.0e-3)*( (1.0e+2)**3.0 ) #kg/m3
	C_mks = C3*1.0e+3 #m/s

	P3 = 1.0e+9*D_mks*C_mks*C_mks/gamma #nPa
	n3 = D3/Mp #Number density, #/cm3

	return n3,P3

#Get data from HDF-4 file
def getHDFVec(hdffile,qi,Scl=1.0):
	qxi = qi+'x_'
	qyi = qi+'y_'
	qzi = qi+'z_'

	Qx3 = hdffile.select(qxi).get()
	Qy3 = hdffile.select(qyi).get()
	Qz3 = hdffile.select(qzi).get()

	#These are too big, corner-sized but corners are poison
	#Chop out corners

	Qx3cc = Scl*Qx3[:-1,:-1,:-1]
	Qy3cc = Scl*Qy3[:-1,:-1,:-1]
	Qz3cc = Scl*Qz3[:-1,:-1,:-1]
	return Qx3cc,Qy3cc,Qz3cc

def getHDFScl(hdffile,q,Scl=1.0):
	qi = q+"_"
	Q3 = np.double(hdffile.select(qi).get())
	#These are too big, corner-sized but corners are poison
	#Chop out corners

	Q3cc = Scl*Q3[:-1,:-1,:-1]
	return Q3cc	