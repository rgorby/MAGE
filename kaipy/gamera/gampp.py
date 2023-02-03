#Gamera post-processing routines
#Get data from serial/MPI gamera output
import glob
import numpy as np
from kaipy.kaiTools import MJD2UT
import itertools
import kaipy.kdefs as kdefs

#Object to use to pull data from HDF5 structure (serial or mpi)

#Initialize,
#gIn = gampp.GameraPipe(fdir,ftag)
#Set grid from data
#gIn.GetGrid()
#V = gIn.GetVar("D",stepnum)
#doFast=True skips various data scraping
idStr = "_0000_0000_0000.gam.h5"

class GameraPipe(object):

	#Initialize GP object
	#fdir = directory to h5 files
	#ftag = stub of h5 files
	def __init__(self,fdir,ftag,doFast=False,doVerbose=True):

		self.fdir = fdir
		self.ftag = ftag

		#Here just create variables
		#-----------
		#Ranks
		self.isMPI = False
		self.Nr = 0
		self.Ri = 1 ; self.Rj = 1 ; self.Rk = 1
		#Global/local cells
		self.is2D = False
		self.Ni = 0 ; self.Nj = 0 ; self.Nk = 0
		self.dNi= 0 ; self.dNj= 0 ; self.dNk= 0

		#Variables/slices
		self.Nv = 0 ; self.Nt = 0 ; self.Nv0 = 0
		self.T = [] ; self.vIDs = [] ; self.v0IDs = []
		self.s0 = 0 ; self.sFin = 0
		self.sids = np.array([])

		#Grids
		self.X = [] ; self.Y = [] ; self.Z = []

		self.doFast = doFast
		self.UnitsID = "NONE"

		#Example file data
		self.f0 = []

		#Stubs for MJD stuff
		self.hasMJD = False
		self.MJDs = []

		#Stub for UT date time object
		self.hasUT = False
		self.UT = []

		#Scrape data from directory
		self.OpenPipe(doVerbose)

	def OpenPipe(self,doVerbose=True):
		import kaipy.kaiH5 as kh5
		if (doVerbose):
			print("Opening pipe: %s : %s"%(self.fdir,self.ftag))

		#Test for serial (either old or new naming convention)
		fOld    = "%s/%s.h5"%(self.fdir,self.ftag)
		fNew    = "%s/%s.gam.h5"%(self.fdir,self.ftag)

		if ( len(glob.glob(fOld)) == 1):
			#Found old-style serial
			self.isMPI = False
			if (doVerbose):
				print("Found serial database")
			f0 = fOld
		elif ( len(glob.glob(fNew)) == 1):
			#Found new-style serial
			self.isMPI = False
			self.ftag = self.ftag + ".gam" #Add .gam to tag
			if (doVerbose):
				print("Found serial database")
			f0 = fNew
		else:
			print("%s not found, looking for MPI database"%(fOld))
			self.isMPI = True
			sStr = "%s/%s_*%s"%(self.fdir,self.ftag,idStr)

			fIns = glob.glob(sStr)
			if (len(fIns)>1):
				print("This shouldn't happen, bailing ...")
			if (len(fIns) == 0):
				print("No MPI database found, all out of options, bailing ...")
				quit()
			f0 = fIns[0]
			Ns = [int(s) for s in f0.split('_') if s.isdigit()]

			self.Ri = Ns[-5]
			self.Rj = Ns[-4]
			self.Rk = Ns[-3]
			self.Nr = self.Ri*self.Rj*self.Rk
			if (doVerbose):
				print("\tFound %d = (%d,%d,%d) ranks"%(self.Nr,self.Ri,self.Rj,self.Rk))

		#In either case, f0 is defined.  use it to get per file stuff
		self.Nt,sids = kh5.cntSteps(f0)
		if (self.doFast):
			self.T = np.zeros(self.Nt)
		else:
			self.T = kh5.getTs(f0,sids)

		self.sids = sids
		self.s0 = sids.min()
		self.sFin = sids.max()
		if (doVerbose):
			print("Found %d timesteps\n\tTime = [%f,%f]"%(self.Nt,self.T.min(),self.T.max()))
			print("\tSteps = [%d,%d]"%(sids.min(),sids.max()))

		#Get MJD if present
		MJDs = kh5.getTs(f0,sids,"MJD",-np.inf)
		if (MJDs.max()>0):
			self.hasMJD = True
			self.MJDs = MJDs
			self.hasUT = True
			self.UT = MJD2UT(MJDs)


		#Get grid stuff
		Dims = kh5.getDims(f0)
		Nd = len(Dims)
		self.dNi = Dims[0]-1
		self.dNj = Dims[1]-1
		if (Nd>2):
			self.dNk = Dims[2]-1

		else:
			self.dNk = 0
			self.Rk = 1
		self.Ni = self.dNi*self.Ri
		self.Nj = self.dNj*self.Rj
		self.Nk = self.dNk*self.Rk


		#Variables
		self.v0IDs = kh5.getRootVars(f0)
		self.vIDs  = kh5.getVars(f0,sids.min())
		self.Nv0 = len(self.v0IDs)
		self.Nv  = len(self.vIDs)
		if (Nd>2):
			self.is2D = False
		else:
			self.is2D = True
		if (doVerbose):
			if (Nd>2):
				nCells = self.Ni*self.Nj*self.Nk
			else:
				nCells = self.Ni*self.Nj
			print("Grid size = (%d,%d,%d)"%(self.Ni,self.Nj,self.Nk))
			print("\tCells = %e"%(nCells))
			print("Variables (Root/Step) = (%d,%d)"%(self.Nv0,self.Nv))
			print("\tRoot: %s"%(self.v0IDs))
			print("\tStep: %s"%(self.vIDs))

		self.SetUnits(f0)
		if (doVerbose):
			print("Units Type = %s"%(self.UnitsID))
			#print("Pulling grid ...")
		self.GetGrid(doVerbose)
		self.f0 = f0

	def SetUnits(self,f0):
		import h5py
		with h5py.File(f0,'r') as hf:
			uID = hf.attrs.get("UnitsID","CODE")
		if (not isinstance(uID,str)):
			self.UnitsID = uID.decode('utf-8')

	def GetGrid(self,doVerbose):
		import kaipy.kaiH5 as kh5
		from alive_progress import alive_bar
		if (self.is2D):
			self.X = np.zeros((self.Ni+1,self.Nj+1))
			self.Y = np.zeros((self.Ni+1,self.Nj+1))
		else:
			self.X = np.zeros((self.Ni+1,self.Nj+1,self.Nk+1))
			self.Y = np.zeros((self.Ni+1,self.Nj+1,self.Nk+1))
			self.Z = np.zeros((self.Ni+1,self.Nj+1,self.Nk+1))
		if (doVerbose):
			#print("Del = (%d,%d,%d)"%(self.dNi,self.dNj,self.dNk))
			titStr = "%s/Grid"%(self.ftag)
		else:
			titStr = None
		NrX = max(self.Nr,1)
		with alive_bar(NrX,title=titStr,length=kdefs.barLen) as bar:
			for (i,j,k) in itertools.product(range(self.Ri),range(self.Rj),range(self.Rk)):
				iS = i *self.dNi
				jS = j *self.dNj
				kS = k *self.dNk
				iE = iS+self.dNi
				jE = jS+self.dNj
				kE = kS+self.dNk
				#print("Bounds = (%d,%d,%d,%d,%d,%d)"%(iS,iE,jS,jE,kS,kE))
				if (self.isMPI):
					fIn = self.fdir + "/" + kh5.genName(self.ftag,i,j,k,self.Ri,self.Rj,self.Rk)
				else:
					fIn = self.fdir + "/" + self.ftag + ".h5"
				if (self.is2D):
					self.X[iS:iE+1,jS:jE+1] = kh5.PullVar(fIn,"X")
					self.Y[iS:iE+1,jS:jE+1] = kh5.PullVar(fIn,"Y")
				else:
					self.X[iS:iE+1,jS:jE+1,kS:kE+1] = kh5.PullVar(fIn,"X")
					self.Y[iS:iE+1,jS:jE+1,kS:kE+1] = kh5.PullVar(fIn,"Y")
					self.Z[iS:iE+1,jS:jE+1,kS:kE+1] = kh5.PullVar(fIn,"Z")
				bar()
	#Get 3D variable "vID" from Step# sID
	def GetVar(self,vID,sID=None,vScl=None,doVerb=True):
		import kaipy.kaiH5 as kh5
		from alive_progress import alive_bar
		if (doVerb):
			if (sID is None):
				titStr = "%s/%s"%(self.ftag,vID)
				
			else:
				titStr = "%s/Step#%d/%s"%(self.ftag,sID,vID)
		else:
			titStr = ''

		if (self.is2D):
			V = np.zeros((self.Ni,self.Nj))
		else:
			V = np.zeros((self.Ni,self.Nj,self.Nk))

		NrX = max(self.Nr,1)
		
		with alive_bar(NrX,title=titStr.ljust(kdefs.barLab),length=kdefs.barLen) as bar:
			for (i,j,k) in itertools.product(range(self.Ri),range(self.Rj),range(self.Rk)):

				iS = i*self.dNi
				jS = j*self.dNj
				kS = k*self.dNk
				iE = iS+self.dNi
				jE = jS+self.dNj
				kE = kS+self.dNk
				#print("Bounds = (%d,%d,%d,%d,%d,%d)"%(iS,iE,jS,jE,kS,kE))
				if (self.isMPI):
					fIn = self.fdir + "/" + kh5.genName(self.ftag,i,j,k,self.Ri,self.Rj,self.Rk)
				else:
					fIn = self.fdir + "/" + self.ftag + ".h5"

				if (self.is2D):
					V[iS:iE,jS:jE] = kh5.PullVar(fIn,vID,sID)

				else:
					V[iS:iE,jS:jE,kS:kE] = kh5.PullVar(fIn,vID,sID)
				bar()
		if (vScl is not None):
			V = vScl*V
		return V

	#Get variable slice of constant i,j,k
	#Directions = idir/jdir/kdir strings
	#Indexing = (1,Nijk)
	#FIXME: Currently pulling whole 3D array and then slicing, lazy
	def GetSlice(self,vID,sID,ijkdir='idir',n=1,vScl=None,doVerb=True):
		sDirs = ['IDIR','JDIR','KDIR']
		ijkdir = ijkdir.upper()

		if (not (ijkdir in sDirs)):
			print("Invalid slice direction, defaulting to I")
			ijkdir = 'IDIR'
		if (sID is None):
			cStr = "Reading %s/%s"%(self.ftag,vID)
		else:
			cStr = "Reading %s/Step#%d/%s"%(self.ftag,sID,vID)
		V = self.GetVar(vID,sID,vScl,doVerb=False)
		#Now slice
		np = n-1 #Convert from Fortran to Python indexing
		if (ijkdir == "IDIR"):
			Vs = V[np,:,:]
			sStr = "(%d,:,:)"%(n)
		elif (ijkdir == "JDIR"):
			Vs = V[:,np,:]
			sStr = "(:,%d,:)"%(n)
		elif (ijkdir == "KDIR"):
			Vs = V[:,:,np]
			sStr = "(:,:,%d)"%(n)
		if (doVerb):
			print(cStr+sStr)
		return Vs

	#Wrappers for root variables (or just set sID to "None")
	def GetRootVar(self,vID,vScl=None,doVerb=True):
		V = self.GetVar(vID,sID=None,vScl=vScl,doVerb=doVerb)
		return V
	def GetRootSlice(self,vID,ijkdir='idir',n=1,vScl=None,doVerb=True):
		Vs = self.GetSlice(vID,sID=None,ijkdir=ijkdir,n=n,vScl=vScl,doVerb=doVerb)
		return Vs
