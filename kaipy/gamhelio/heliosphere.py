#Various tools to post-process and analyze Gamera heliosphere runs
from kaipy.kdefs import *
import kaipy.gamera.gampp
from kaipy.gamera.gampp import GameraPipe
import numpy as np
import glob
import kaipy.kaiH5 as kh5
import timeit

#Object to pull from MPI/Serial heliosphere runs (H5 data), extends base

ffam =  "monospace"
dLabC = "black" #Default label color
dLabFS = "medium" #Default label size
dBoxC = "lightgrey" #Default box color
TINY = 1.0e-8
rmStr = "mixtest"

#Assuming LFM/EGG type grid -need to adapt to helio grid
class GamsphPipe(GameraPipe):
	#Initialize object, rely on base class, take optional unit identifier
	def __init__(self,fdir,ftag,doFast=False,uID="Inner"):

		print("Initializing %s heliosphere"%(uID))
		
		#units for helio
		#self.MagM = -EarthM0g*1.0e+5
		self.bScl = 100.    #->nT
		#self.pScl = 1.67e-2 #->nPa
		self.vScl = 150.  #-> km/s
		self.tScl = 4637.    #->seconds
		self.dScl = 200. #cm-3
		self.TScl = 1.e-6/4/np.pi/200./1.38e-16/2./1.e6 #in K

		#2D equatorial grid, stretched polar (Ni,Nj*2+1)
		#??????
		#2D equatorial grid
		self.xxi = [] ; self.yyi = [] #corners
		self.xxc = [] ; self.yyc = [] #centers

		#2D merid grid
		self.xxi_m = []; self.zzi = [] #corners
		self.xxc_m = []; self.zzc = []

		#base class, will use OpenPipe below
		GameraPipe.__init__(self,fdir,ftag,doFast=doFast)

		self.R0 = self.xxc[0,0]
		
	def OpenPipe(self,doVerbose=True):
		GameraPipe.OpenPipe(self,doVerbose)
		
		if (self.UnitsID != "CODE"):
			self.bScl   = 1.0  #->nT
			#self.pScl   = 1.0  #->nPa
			self.vScl   = 1.0  #-> km/s
			self.tScl   = 1.0 #->Seconds
			# [EP] added
			self.dScl = 1.0
			self.TScl = 1.0 
			#self.tRMScl = 63.8 #->Seconds

		#Rescale time
		self.T = self.tScl*self.T

		Neq_a = self.Nj//2 #cell above eq plane
		print (Neq_a)

		Nr = self.Ni
		Np = self.Nk

		#corners in eq XY plane
		self.xxi = np.zeros((Nr+1,Np+1)) 
		self.yyi = np.zeros((Nr+1,Np+1))
		#centers
		self.xxc = np.zeros((Nr  ,Np  ))
		self.yyc = np.zeros((Nr  ,Np  ))
		

		#equatorial plane
		#corners i,k in eq plane, j index is Neq_a
		self.xxi[:,:] = self.X[:,Neq_a,:]
		self.yyi[:,:] = self.Y[:,Neq_a,:]

		#centers i,k
		self.xxc = 0.25*(self.xxi[:-1,:-1] + self.xxi[1:,:-1] + self.xxi[:-1,1:] + self.xxi[1:,1:])
		self.yyc = 0.25*(self.yyi[:-1,:-1] + self.yyi[1:,:-1] + self.yyi[:-1,1:] + self.yyi[1:,1:])
		r = np.sqrt(self.xxc**2.0 + self.yyc**2.0)

		#!!!for helio meridional plane grid would be different; should be done like for magnetosphere


		#Upper half plane
		#for j in range(self.Nj):
		#	self.xxi[:,j] = self.X[:,j,0] #k=0
		#	self.yyi[:,j] = self.Y[:,j,0]
		#Lower half plane
		#for j in range(self.Nj,Np+1):
		#	jp = Np-j
		#	self.xxi[:,j] =  self.X[:,jp,0]
		#	self.yyi[:,j] = -self.Y[:,jp,0]

		if (self.hasMJD):
			print("Found MJD data")
			print("\tTime (Min/Max) = %f/%f"%(self.MJDs.min(),self.MJDs.max()))

	def EqSlice(self,vID,sID=None,vScl=None,doEq=True,doVerb=True):
		#Get full 3D variable first
		Q = self.GetVar(vID,sID,vScl,doVerb)

		Nj2 = self.Nj//2 

		#above and below the eq plane
		ja = Nj2 - 1
		jb = ja + 1

		Nr = self.Ni
		Np = self.Nk

		#equatorial j-slice of var 
		Qj = np.zeros((Nr,Np))
		#taking average above/below eq plane
		Qj[:,:] = 0.5*( Q[:,ja,:] + Q[:,jb,:] )
		return Qj

	#merid plane Y=0
	def MeridGrid(self):
		#Get Grid
		self.GetGrid(doVerbose=True)

		Nk2 = self.Nk//2
		Nt = self.Nj
		
		#kooking from -Y to XZ plane
		xright = self.X[:,:,0] #corners
		xleft = self.X [:,:,Nk2]

		zright = self.Z[:,:,0] #corners
		zleft = self.Z[:,:,Nk2]

		#stack right and left together
		xmer = np.hstack( (xright, xleft[:,::-1]) ) #reverse j 
		zmer = np.hstack( (zright, zleft[:,::-1]) ) #reverse j

		#cell centers
		xmer_c = 0.25*( xmer[:-1,:-1]+xmer[:-1,1:]+xmer[1:,:-1]+xmer[1:,1:] )
		xmer_c = np.delete(xmer_c, Nt, axis = 1)
		zmer_c = 0.25*( zmer[:-1,:-1]+zmer[:-1,1:]+zmer[1:,:-1]+zmer[1:,1:] )
		zmer_c = np.delete(zmer_c, Nt, axis = 1)
		return xmer_c, zmer_c

	#merid plane Y=0
	def MeridGridHalfs(self):
		self.GetGrid(doVerbose=True)

		Nk2 = self.Nk//2
		Nt = self.Nj

		#looking from -Y to XZ plane
		xright = self.X[:,:,0] #corners
		zright = self.Z[:,:,0] #corners

		xleft = self.X [:,:,Nk2]
		zleft = self.Z[:,:,Nk2]

		xright_c = 0.25*( xright[:-1,:-1]+xright[:-1,1:]+xright[1:,:-1]+xright[1:,1:] )
		zright_c = 0.25*( zright[:-1,:-1]+zright[:-1,1:]+zright[1:,:-1]+zright[1:,1:] )
		r = np.sqrt(xright_c**2 + zright_c**2)

		#centers: right plane, left plane, radius
		return xright, zright, xleft, zleft, r

	def iSliceGrid(self):
		#Get Grid
		self.GetGrid(doVerbose=True)

		rxy = np.sqrt(self.X**2 + self.Y**2)
		theta = np.arctan2(rxy,self.Z)
		phi = np.arctan2(self.Y,self.X)

		theta = 90. - theta*180./np.pi
		phi [phi < 0] += 2*np.pi
		phi = phi*180./np.pi

		#last i-index == face of the last cell
		lat = theta[-1,:,:]
		lon = phi[-1,:,:]

		return lat, lon

	def MeridSlice(self,vID,sID=None,vScl=None,doVerb=True):
		#Get full 3D variable first
		Q = self.GetVar(vID,sID,vScl,doVerb)
		
		Nk2 = self.Nk//2
		Np = self.Nk
		
		#Nr = self.Ni
		#Nt = 2*self.Nj
		#XZ meridional slice (k=0) of var 
		#Qj = np.zeros((Nr,Nt))
		
		Qright = 0.5*( Q[:,:,0] + Q[:,:,Np-1] ) 
		Qleft  = 0.5*( Q[:,:,Nk2-1] + Q[:,:,Nk2] )
		#print (Qright.shape, Qleft.shape)
		#Qj = np.hstack( (Qright, Qleft[:,::-1]) ) #reverse in j
		#print (Qj.shape)
		return Qright, Qleft

	def iSliceVar(self,vID,sID=None,vScl=None,doVerb=True):
		#Get full 3D variable first
		Q = self.GetVar(vID,sID,vScl,doVerb)

		#cell centered valies from the last cell
		Qi = Q[-1,:,:]
		return Qi

	def iSliceMagV(self,s0=0):
		Vx = self.iSliceVar("Vx",s0) #Unscaled
		Vy = self.iSliceVar("Vy",s0) #Unscaled
		Vz = self.iSliceVar("Vz",s0) #Unscaled
		Vi = self.vScl*np.sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
		return Vi

	def iSliceD(self,s0=0):
		Di = self.iSliceVar("D",s0) #Unscaled
		Di = Di*self.dScl
		return Di

	#Equatorial speed (in km/s) in eq plane
	def eqMagV(self,s0=0):
		Vx = self.EqSlice("Vx",s0) #Unscaled
		Vy = self.EqSlice("Vy",s0) #Unscaled
		Vz = self.EqSlice("Vz",s0) #Unscaled
		Veq = self.vScl*np.sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
		return Veq

	#Normalized density (D*r*r/21.5/21.5 in cm-3) in eq plane
	def eqNormD (self,s0=0):
		#self.GetGrid(self,doVerb=True)

		D = self.EqSlice("D",s0) #Unscaled
		#radius of inner boudnary (center of the first cell)
		R0 = self.xxc[0,0]
		print ("R0 = %f"%R0)

		#calculate normalized density
		#r = np.sqrt(self.xxc**2.0 + self.yyc**2.0)
		Norm = (self.xxc**2.0 + self.yyc**2.0)/R0/R0
		NormDeq = self.dScl*D*Norm
		return NormDeq

	#Normalized Br (Br*r*r/21.5/21.5) in eq plane
	def eqNormBr (self,s0=0):
		Bx = self.EqSlice("Bx",s0) #Unscaled
		By = self.EqSlice("By",s0) #Unscaled
		Bz = self.EqSlice("Bz",s0) #Unscaled
		#Br = self.EqSlice("Br",s0) #Unscaled
		#R0 = self.xxc[0,0]

		#calculate normalized Br
		#r = np.sqrt(self.xxc**2.0 + self.yyc**2.0)
		Br = (Bx*self.xxc + By*self.yyc)*np.sqrt(self.xxc**2.0 + self.yyc**2.0)/self.R0/self.R0
		
		NormBreq = self.bScl*Br
		return NormBreq

	#Temperature
	def eqTemp (self,s0=0):
		Pres = self.EqSlice("P",s0)
		D = self.EqSlice("D",s0)

		Temp = Pres/D*self.TScl
		
		return Temp
	
	#Meridional speed (in km/s)
	def MerMagV(self,s0=0):
		Vxr, Vxl = self.MeridSlice("Vx",s0) #Unscaled
		Vyr, Vyl = self.MeridSlice("Vy",s0) #Unscaled
		Vzr, Vzl = self.MeridSlice("Vz",s0) #Unscaled
		MagVr = self.vScl*np.sqrt(Vxr**2.0+Vyr**2.0+Vzr**2.0)
		MagVl = self.vScl*np.sqrt(Vxl**2.0+Vyl**2.0+Vzl**2.0)
		return MagVr, MagVl

	def MerDNrm(self,s0=0):
		xr, zr, xl, zl, r = self.MeridGridHalfs()
		print (r.shape)
		Dr, Dl = self.MeridSlice("D",s0) #Unscaled
		Drn = Dr*self.dScl*r*r/self.R0/self.R0
		Dln = Dl*self.dScl*r*r/self.R0/self.R0
		return Drn, Dln


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

	#def AddSW(self,n,Ax,xy=[0.725,0.025],cLab=dLabC,fs=dLabFS,T0=0.0,doBox=True,BoxC=dBoxC,doAll=True):
	#	import kaipy.kaiH5 as kh5
	#	#Start by getting SW data
	#	vIDs = ["D","P","Vx","Bx","By","Bz"]
	#	Nv = len(vIDs)
	#	qSW = np.zeros(Nv)
	#	if (self.isMPI):
	#		fSW = self.fdir + "/" + kh5.genName(self.ftag,self.Ri-1,0,0,self.Ri,self.Rj,self.Rk)
	#	else:
	#		fSW = self.fdir + "/" + self.ftag + ".h5"
	#
	#	for i in range(Nv):
	#		Q = kh5.PullVar(fSW,vIDs[i],n)
	#		qSW[i] = Q[-1,0,0]
	#	D = qSW[0] ; P = qSW[1] ; Vx = qSW[2] ; Bx = qSW[3] ; By = qSW[4] ; Bz = qSW[5]
	#	SWStr = "Solar Wind\n"
	#	MagB = self.bScl*np.sqrt(Bx**2.0+By**2.0+Bz**2.0)
	#	#Clock = atan(by/bz), cone = acos(Bx/B)
	#	r2deg = 180.0/np.pi
	#	if (MagB>TINY):
	#		clk  = r2deg*np.arctan2(By,Bz)
	#		cone = r2deg*np.arccos(self.bScl*Bx/MagB)
	#	else:
	#		clk = 0.0
	#		cone = 0.0
	#	if (clk<0):
	#		clk = clk+360.0
	#	Deg = r"$\degree$"
	#	SWStr = "Solar Wind\nIMF: %4.1f [nT], %5.1f"%(MagB,clk) + Deg# + ", %5.2f"%(cone) + Deg
	#	if (doAll):
	#		SWStr = SWStr + "\nDensity: %5.1f [#/cc] \nSpeed:  %6.1f [km/s] "%(D,self.vScl*np.abs(Vx))

	#	if (doBox):
	#		Ax.text(xy[0],xy[1],SWStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam,bbox=dict(boxstyle="round",fc=dBoxC))
	#	else:
	#		Ax.text(xy[0],xy[1],SWStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam,bbox=dict(boxstyle="round",fc=dBoxC))
	#def AddCPCP(self,n,Ax,xy=[0.9,0.95],cLab=dLabC,fs=dLabFS,doBox=True,BoxC=dBoxC):
	#	cpcp = self.GetCPCP(n)
	#	tStr = "CPCP   (North/South)\n%6.2f / %6.2f [kV]"%(cpcp[0],cpcp[1])
	#	if (doBox):
	#		Ax.text(xy[0],xy[1],tStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam,bbox=dict(boxstyle="round",fc=dBoxC))
	#	else:
	#		Ax.text(xy[0],xy[1],tStr,color=cLab,fontsize=fs,transform=Ax.transAxes,family=ffam)

	
