#Various tools to post-process and analyze Gamera magnetosphere runs
from kaipy.kdefs import *
import kaipy.gamera.gampp
from kaipy.gamera.gampp import GameraPipe
import numpy as np
import glob
import kaipy.kaiH5 as kh5
import timeit

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

		#New mix stuff
		self.hasRemix  = False #Remix data present (new format)
		self.mixPipe = None

		#Old mix stuff
		self.tRMScl = 63.8 #->seconds, for Remix time scaling
		
		self.hasRemixO = False #Remix data present (old format)
		self.Nrm = 0 #Number of remix outputs

		#2D equatorial grid, stretched polar (Ni,Nj*2+1)
		self.xxi = [] ; self.yyi = []
		self.xxc = [] ; self.yyc = []

		GameraPipe.__init__(self,fdir,ftag,doFast=doFast)

		self.Rin = self.xxi[0,0]
		
	def OpenPipe(self,doVerbose=True):
		GameraPipe.OpenPipe(self,doVerbose)
		self.GetM0()

		#Now do magnetosphere specific things
		if (self.UnitsID != "CODE"):
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
		#Upper half plane
		for j in range(self.Nj):
			self.xxi[:,j] = self.X[:,j,0]
			self.yyi[:,j] = self.Y[:,j,0]
		#Lower half plane
		for j in range(self.Nj,Np+1):
			jp = Np-j
			self.xxi[:,j] =  self.X[:,jp,0]
			self.yyi[:,j] = -self.Y[:,jp,0]

		#Get centers for stretched polar grid & BzD
		self.xxc = 0.25*(self.xxi[:-1,:-1] + self.xxi[1:,:-1] + self.xxi[:-1,1:] + self.xxi[1:,1:])
		self.yyc = 0.25*(self.yyi[:-1,:-1] + self.yyi[1:,:-1] + self.yyi[:-1,1:] + self.yyi[1:,1:])
		r = np.sqrt(self.xxc**2.0 + self.yyc**2.0)
		rm5 = r**(-5.0)
		self.BzD = -r*r*self.MagM*rm5

		if (self.hasMJD):
			print("Found MJD data")
			print("\tTime (Min/Max) = %f/%f"%(self.MJDs.min(),self.MJDs.max()))

	#Do remix data things
		#Old-style
		rmOStr = "%s/%s*.h5"%(self.fdir,rmStr)
		rmOuts = glob.glob(rmOStr)
		Nrm = len(rmOuts)
		
		if (Nrm>0):
			import h5py
			print("Found %d ReMIX outputs"%(Nrm))
			self.hasRemixO = True
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
		#New-style
		Stubs = self.ftag.split(".")
		rmOStr = "%s/%s.mix.h5"%(self.fdir,Stubs[0])
		rmOuts = glob.glob(rmOStr)
		Nrm = len(rmOuts)
		if (Nrm>0):
			mixtag = Stubs[0]+".mix"
			self.hasRemix = True
			print("Found ReMIX data, reading ...")
			self.mixPipe = GameraPipe(self.fdir,mixtag,doVerbose=False)
			self.nCPCP = kh5.getTs(rmOStr,sIds=self.sids,aID="nCPCP")
			self.sCPCP = kh5.getTs(rmOStr,sIds=self.sids,aID="sCPCP")
	#Get magnetic moment from file
	def GetM0(self):
		import h5py
		with h5py.File(self.f0,'r') as hf:
			M0 = hf.attrs.get("MagM0",self.MagM)
		#Store value
		self.MagM = M0
		
	#Get "egg" slice, variable matched to stretched polar grid
	#Either equatorial or meridional
	def EggSlice(self,vID,sID=None,vScl=None,doEq=True,doVerb=True):
		#Get full 3D variable first
		Q = self.GetVar(vID,sID,vScl,doVerb)

		#For upper/lower half planes, average above/below
		Nk2 = self.Nk//2
		Nk4 = self.Nk//4
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
		if (self.hasRemixO):
			#Find nearest time slice
			i0 = np.abs(self.tRm-tGam).argmin()
			fMix = self.rmOuts[i0]
		else:
			fMix = None
		return fMix
		
	#Get CPCP @ gamera step #n
	def GetCPCP(self,n):
		n0 = n-self.s0
		tGam = self.T[n0]
		if (self.hasRemixO): #Old remix style
			#Find nearest time slice
			i0 = np.abs(self.tRm-tGam).argmin()
			cpcp = [self.nCPCP[i0],self.sCPCP[i0]]
		elif (self.hasRemix):
			cpcp = [self.nCPCP[n0],self.sCPCP[n0]]
		else:
			cpcp = [0.0,0.0]
		return cpcp

	#Add time label, xy is position in axis (not data) coords
	def AddTime(self,n,Ax,xy=[0.9,0.95],cLab=dLabC,fs=dLabFS,T0=0.0,doBox=True,BoxC=dBoxC):
		ffam = "monospace"
		HUGE = 1.0e+8
		#Decide whether to do UT or elapsed
		if (self.hasMJD):
			minMJD = self.MJDs[n-self.s0]
		else:
			minMJD = -HUGE
		if (self.hasMJD and minMJD>TINY):
			from astropy.time import Time
			dtObj = Time(self.MJDs[n-self.s0],format='mjd').datetime
			tStr = "  " + dtObj.strftime("%H:%M:%S") + "\n" + dtObj.strftime("%m/%d/%Y")

		else:	
			#Get time in seconds
			t = self.T[n-self.s0] - T0
			Nm = np.int( (t-T0)/60.0 ) #Minutes, integer
			Hr = Nm/60
			Min = np.mod(Nm,60)
			Sec = np.mod(np.int(t),60)

			tStr = "Elapsed Time\n  %02d:%02d:%02d"%(Hr,Min,Sec)
		if (doBox):
			Ax.text(xy[0],xy[1],tStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam,bbox=dict(boxstyle="round",fc=dBoxC))

	def AddSW(self,n,Ax,xy=[0.725,0.025],cLab=dLabC,fs=dLabFS,T0=0.0,doBox=True,BoxC=dBoxC,doAll=True):
		import kaipy.kaiH5 as kh5
		#Start by getting SW data
		vIDs = ["D","P","Vx","Bx","By","Bz"]
		Nv = len(vIDs)
		qSW = np.zeros(Nv)
		if (self.isMPI):
			fSW = self.fdir + "/" + kh5.genName(self.ftag,self.Ri-1,0,0,self.Ri,self.Rj,self.Rk)
		else:
			fSW = self.fdir + "/" + self.ftag + ".h5"

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

	#Replacement for remixpp adding inset remix plots
	#Becomes a wrapper to remixpp.CMIViz
	def CMIViz(self,AxM=None,nStp=0,doNorth=True,loc="upper left",dxy=[20,20]):
		from kaipy.kaiH5 import CheckOrDie
		import h5py
		import kaipy.gamera.remixpp as remixpp
		if (doNorth):
			pID = "Potential NORTH"
			cID = "Field-aligned current NORTH"
			gID = "NORTH"
		else:
			pID = "Potential SOUTH"
			cID = "Field-aligned current SOUTH"
			gID = "SOUTH"

		if (self.hasRemixO): #Old remix style
			fMix = self.Gam2Remix(nStp)
			CheckOrDie(fMix)
			with h5py.File(fMix,'r') as hf:
				P = hf[gID]['Potential'][()]
				C = hf[gID]['Field-aligned current'][()]

		elif (self.hasRemix): #New style remix
			P  = self.mixPipe.GetVar(pID,nStp,doVerb=False).T
			C  = self.mixPipe.GetVar(cID,nStp,doVerb=False).T
			
		RIn = self.xxi[0,0]
		llBC = np.arcsin(np.sqrt(1.0/RIn))*180.0/np.pi
		nLat,nLon = C.shape
		#Now call remixpp.CMIViz
		remixpp.CMIPic(nLat,nLon,llBC,P,C,AxM,doNorth,loc,dxy)
