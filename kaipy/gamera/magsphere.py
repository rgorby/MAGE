#Various tools to post-process and analyze Gamera magnetosphere runs
from kaipy.kdefs import *
import gampp
from gampp import GameraPipe
import numpy as np
import glob
#Object to pull from MPI/Serial magnetosphere runs (H5 data), extends base

ffam =  "monospace"
dLabC = "black" #Default label color
dLabFS = "medium" #Default label size
dBoxC = "lightgrey" #Default box color
TINY = 1.0e-8
rmStr = "mixtest"

#Assuming LFM/EGG type grid
class GamsphPipe(GameraPipe):
	#Initialize object, rely on base class, take optional unit identifier
	def __init__(self,fdir,ftag,doFast=False,uID="Earth"):

		print("Initializing %s magnetosphere"%(uID))
		#TODO: Add different unit/planet options here
		self.MagM = -EarthM0g*1.0e+5
		self.bScl = 4.58    #->nT
		self.pScl = 1.67e-2 #->nPa
		self.vScl = 1.0e+2  #-> km/s
		self.tScl = 63.8    #->seconds
		self.tRMScl = 63.8 #->seconds, for Remix time scaling

		self.hasRemix = False #Remix data present
		self.Nrm = 0 #Number of remix outputs
		#2D equatorial grid, stretched polar (Ni,Nj*2+1)
		self.xxi = [] ; self.yyi = []
		self.xxc = [] ; self.yyc = []

		GameraPipe.__init__(self,fdir,ftag,doFast=doFast)

		self.Rin = self.xxi[0,0]
		
	def OpenPipe(self):
		GameraPipe.OpenPipe(self)

		#Now do magnetosphere specific things
		if (self.UnitsID == "EARTH"):
			self.bScl   = 1.0  #->nT
			self.pScl   = 1.0  #->nPa
			self.vScl   = 1.0  #-> km/s
			self.tScl   = 1.0 #->Seconds
			self.tRMScl = 63.8 #->Seconds

		#Rescale time
		self.T = self.tScl*self.T

		#Create warped polar slice grid
		Nr = self.Ni
		Np = 2*(self.Nj)
		self.xxi = np.zeros((Nr+1,Np+1))
		self.yyi = np.zeros((Nr+1,Np+1))
		self.xxc = np.zeros((Nr  ,Np  ))
		self.yyc = np.zeros((Nr  ,Np  ))
		self.BzD = np.zeros((Nr  ,Np  ))
		#Create corners for stretched polar grid
		for i in range(Nr+1):
			for j in range(self.Nj):
				self.xxi[i,j] = self.X[i,j,0]
				self.yyi[i,j] = self.Y[i,j,0]
			for j in range(self.Nj,Np+1):
				jp = Np-j
				self.xxi[i,j] =  self.X[i,jp,0]
				self.yyi[i,j] = -self.Y[i,jp,0]
		#Get centers for stretched polar grid & BzD
		for i in range(Nr):
			for j in range(Np):
				self.xxc[i,j] = 0.25*(self.xxi[i,j]+self.xxi[i+1,j]+self.xxi[i,j+1]+self.xxi[i+1,j+1])
				self.yyc[i,j] = 0.25*(self.yyi[i,j]+self.yyi[i+1,j]+self.yyi[i,j+1]+self.yyi[i+1,j+1])
				r = np.sqrt(self.xxc[i,j]**2.0+self.yyc[i,j]**2.0)
				rm5 = r**(-5.0)
				self.BzD[i,j] = -r*r*self.MagM*rm5

		#Do remix data things
		rmOStr = "%s/%s*.h5"%(self.fdir,rmStr)
		rmOuts = glob.glob(rmOStr)
		Nrm = len(rmOuts)
		print("Found %d ReMIX outputs"%(Nrm))
		if (Nrm>0):
			import h5py
			self.hasRemix = True
			self.Nrm = Nrm
			self.tRm = np.zeros(Nrm)
			self.nCPCP = np.zeros(Nrm)
			self.sCPCP = np.zeros(Nrm)
			self.rmOuts = rmOuts

			if (not self.doFast):
				for i in range(Nrm):
					fMix = rmOuts[i]
					with h5py.File(fMix,'r') as hf:
						Atts = hf.attrs.keys()
						if ('t' in Atts):
							self.tRm[i] = self.tRMScl*hf.attrs['t']
						if ('nCPCP' in Atts):
							self.nCPCP[i] = hf.attrs['nCPCP']
							self.sCPCP[i] = hf.attrs['sCPCP']
				print("\tTime (Min/Max) = %f/%f"%(self.tRm.min(),self.tRm.max()))
				#Sort into time ordering
				I = self.tRm.argsort()
				self.tRm = self.tRm[I]
				self.nCPCP = self.nCPCP[I]
				self.sCPCP = self.sCPCP[I]
				self.rmOuts = [rmOuts[i] for i in I]

	#Get "egg" slice, variable matched to stretched polar grid
	#Either equatorial or meridional
	def EggSlice(self,vID,sID=None,vScl=None,doEq=True,doVerb=True):
		#Get full 3D variable first
		Q = self.GetVar(vID,sID,vScl,doVerb)

		#For upper/lower half planes, average above/below
		Nk2 = self.Nk/2
		Nk4 = self.Nk/4
		if (doEq):
			ku = -1
			kl = Nk2-1
		else:
			ku = Nk4 - 1
			kl = 3*Nk4-1
		Nr = self.Ni
		Np = 2*self.Nj
		Qk = np.zeros((Nr,Np))
		for j in range(Np):
			if (j>=self.Nj):
				#Lower half plane
				jp = Np-j-1
				Qk[:,j] = 0.5*( Q[:,jp,kl] + Q[:,jp,kl+1] )
			else:
				jp = j
				Qk[:,j] = 0.5*( Q[:,jp,ku] + Q[:,jp,ku+1] )

		return Qk
	#Standard equatorial dbz (in nT)
	def DelBz(self,s0=0):

		Bz = self.EggSlice("Bz",s0) #Unscaled

		dbz = Bz*self.bScl - self.BzD
		return dbz
		
	#Equatorial magnitude of field (in nT)
	def eqMagB(self,s0=0):
		Bx = self.EggSlice("Bx",s0) #Unscaled
		By = self.EggSlice("By",s0) #Unscaled
		Bz = self.EggSlice("Bz",s0) #Unscaled
		Beq = self.bScl*np.sqrt(Bx**2.0+By**2.0+Bz**2.0)
		return Beq

	#Return data for meridional 2D field lines
	#Need to use Cartesian grid
	def bStream(self,s0=0,xyBds=[-35,25,-25,25],dx=0.05):
		
		#Get field data
		U = self.bScl*self.EggSlice("Bx",s0,doEq=False)
		V = self.bScl*self.EggSlice("Bz",s0,doEq=False)

		x1,y1,gu,gv,gM = self.doStream(U,V,xyBds,dx)
		return x1,y1,gu,gv,gM

	def vStream(self,s0=0,xyBds=[-35,25,-25,25],dx=0.05):
		#Get field data
		U = self.vScl*self.EggSlice("Vx",s0,doEq=True)
		V = self.vScl*self.EggSlice("Vy",s0,doEq=True)

		x1,y1,gu,gv,gM = self.doStream(U,V,xyBds,dx)
		return x1,y1,gu,gv,gM

	#Map from Gamera step to remix file
	def Gam2Remix(self,n):
		tGam = self.T[n-self.s0]
		if (self.hasRemix):
			#Find nearest time slice
			i0 = np.abs(self.tRm-tGam).argmin()
			fMix = self.rmOuts[i0]
		else:
			fMix = ""
		return fMix
	def GetCPCP(self,n):
		tGam = self.T[n-self.s0]
		cpcp = [0.0,0.0]
		if (self.hasRemix):
			#Find nearest time slice
			i0 = np.abs(self.tRm-tGam).argmin()
			cpcp = [self.nCPCP[i0],self.sCPCP[i0]]
		return cpcp

	#Add time label, xy is position in axis (not data) coords
	def AddTime(self,n,Ax,xy=[0.9,0.95],cLab=dLabC,fs=dLabFS,T0=0.0,doBox=True,BoxC=dBoxC):
		ffam = "monospace"
		#Get time in seconds
		t = self.T[n-self.s0] - T0
		Nm = np.int( (t-T0)/60.0 ) #Minutes, integer
		Hr = Nm/60
		Min = np.mod(Nm,60)
		Sec = np.mod(np.int(t),60)

		tStr = "Time %02d:%02d:%02d"%(Hr,Min,Sec)
		if (doBox):
			Ax.text(xy[0],xy[1],tStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam,bbox=dict(boxstyle="round",fc=dBoxC))

	def AddSW(self,n,Ax,xy=[0.725,0.025],cLab=dLabC,fs=dLabFS,T0=0.0,doBox=True,BoxC=dBoxC,doAll=True):
		import kaipy.kaiH5 as kh5
		#Start by getting SW data
		vIDs = ["D","P","Vx","Bx","By","Bz"]
		Nv = len(vIDs)
		qSW = np.zeros(Nv)

		fSW = self.ChunkName(self.Ri-1,0,0)
		for i in range(Nv):
			Q = kh5.PullVar(fSW,vIDs[i],n)
			qSW[i] = Q[-1,0,0]
		D = qSW[0] ; P = qSW[1] ; Vx = qSW[2] ; Bx = qSW[3] ; By = qSW[4] ; Bz = qSW[5]
		SWStr = "Solar Wind\n"
		MagB = self.bScl*np.sqrt(Bx**2.0+By**2.0+Bz**2.0)
		#Clock = atan(by/bz), cone = acos(Bx/B)
		r2deg = 180.0/np.pi
		if (MagB>TINY):
			clk  = r2deg*np.arctan2(By,Bz)
			cone = r2deg*np.arccos(self.bScl*Bx/MagB)
		else:
			clk = 0.0
			cone = 0.0
		if (clk<0):
			clk = clk+360.0
		Deg = r"$\degree$"
		SWStr = "Solar Wind\nIMF: %4.1f [nT], %5.1f"%(MagB,clk) + Deg# + ", %5.2f"%(cone) + Deg
		if (doAll):
			SWStr = SWStr + "\nDensity: %5.1f [#/cc] \nSpeed:  %6.1f [km/s] "%(D,self.vScl*np.abs(Vx))

		if (doBox):
			Ax.text(xy[0],xy[1],SWStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam,bbox=dict(boxstyle="round",fc=dBoxC))
		else:
			Ax.text(xy[0],xy[1],SWStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam,bbox=dict(boxstyle="round",fc=dBoxC))
	def AddCPCP(self,n,Ax,xy=[0.9,0.95],cLab=dLabC,fs=dLabFS,doBox=True,BoxC=dBoxC):
		cpcp = self.GetCPCP(n)
		tStr = "CPCP   (North/South)\n%6.2f / %6.2f [kV]"%(cpcp[0],cpcp[1])
		if (doBox):
			Ax.text(xy[0],xy[1],tStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam,bbox=dict(boxstyle="round",fc=dBoxC))
		else:
			Ax.text(xy[0],xy[1],tStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam)

	def doStream(self,U,V,xyBds=[-35,25,-25,25],dx=0.05):
		from scipy.interpolate import griddata
		N1 = np.int( (xyBds[1]-xyBds[0])/dx )
		N2 = np.int( (xyBds[3]-xyBds[2])/dx )

		#Create matching Cartesian grid
		x1 = np.linspace(xyBds[0],xyBds[1],N1)
		y1 = np.linspace(xyBds[2],xyBds[3],N2)

		xx1,yy1 = np.meshgrid(x1,y1)
		r1 = np.sqrt(xx1**2.0+yy1**2.0)
		#Flatten and interpolate
		px = self.xxc.flatten()
		py = self.yyc.flatten()
		pu = U.flatten()
		pv = V.flatten()
		pMag = np.sqrt(U**2.0+V**2.0).flatten()

		gu = griddata(zip(px,py),pu  ,(xx1,yy1),method='linear',fill_value=0.0)
		gv = griddata(zip(px,py),pv  ,(xx1,yy1),method='linear',fill_value=0.0)
		gM = griddata(zip(px,py),pMag,(xx1,yy1),method='linear',fill_value=0.0)

		kOut = (r1<=self.Rin*1.01)
		gu[kOut] = 0.0
		gv[kOut] = 0.0
		gM[kOut] = 0.0

		return x1,y1,gu,gv,gM
