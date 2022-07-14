#Various tools to post-process and analyze Gamera heliosphere runs
from kaipy.kdefs import *
import kaipy.gamera.gampp
from kaipy.gamera.gampp import GameraPipe
import numpy as np
import glob
import kaipy.kaiH5 as kh5
import timeit

#Object to pull from MPI/Serial heliosphere runs (H5 data), extends base

ffam   =  "monospace"
dLabC  = "black" #Default label color
dLabFS = "medium" #Default label size
dBoxC  = "lightgrey" #Default box color
TINY   = 1.0e-8
MK     = 1.e6 #MegaKelvin

#Adapted to helio grid
class GamsphPipe(GameraPipe):
	#Initialize object, rely on base class, take optional unit identifier
	def __init__(self,fdir,ftag,doFast=False,uID="Inner"):

		print("Initializing %s heliosphere"%(uID))
		
		#units for inner helio
		self.bScl = 100.    #->nT
		self.vScl = 150.  #-> km/s
		self.tScl = 4637.    #->seconds
		self.dScl = 200. #cm-3
		self.TScl = 1.e-6/4/np.pi/200./kbltz/MK #in MK
    
        # units for OHelio
		#self.bScl = 5.    #->nT
		#self.vScl = 34.5  #-> km/s
		#self.tScl = 1.4e8/34.5
		#self.dScl = 10. #cm-3
		#self.TScl = 0.144 #in MK

		#2D equatorial grid
		self.xxi = [] ; self.yyi = [] #corners
		self.xxc = [] ; self.yyc = [] #centers

		#base class, will use OpenPipe below
		GameraPipe.__init__(self,fdir,ftag,doFast=doFast)

		#inner boundary distance
		self.R0 = self.xxc[0,0]

		#j and k for radial profile
		self.jRad = self.Nj//2
		self.kRad = self.Nk//4

	def OpenPipe(self,doVerbose=True):
		GameraPipe.OpenPipe(self,doVerbose)
		
		if (self.UnitsID != "CODE"):
			self.bScl   = 1.0          #->nT
			self.vScl   = 1.0          #-> km/s
			self.tScl   = 1.0          #-> Seconds
			self.dScl   = 1.0          #-> cm-3
			self.TScl   = 1.0/kbltz/MK #-> MKelvin

		#Rescale time
		self.T = self.tScl*self.T

		Neq_a = self.Nj//2 #cell above eq plane

		Nr = self.Ni
		Np = self.Nk

		#corners in eq XY plane
		self.xxi = np.zeros((Nr+1,Np+1)) 
		self.yyi = np.zeros((Nr+1,Np+1))
		#centers
		self.xxc = np.zeros((Nr  ,Np  ))
		self.yyc = np.zeros((Nr  ,Np  ))
		
		#Grid for equatorial plane. Should probably be done as a separate function
		#equatorial plane
		#corners i,k in eq plane, j index is Neq_a
		self.xxi[:,:] = self.X[:,Neq_a,:]
		self.yyi[:,:] = self.Y[:,Neq_a,:]

		#centers i,k
		self.xxc = 0.25*(self.xxi[:-1,:-1] + self.xxi[1:,:-1] + self.xxi[:-1,1:] + self.xxi[1:,1:])
		self.yyc = 0.25*(self.yyi[:-1,:-1] + self.yyi[1:,:-1] + self.yyi[:-1,1:] + self.yyi[1:,1:])
		r = np.sqrt(self.xxc**2.0 + self.yyc**2.0)

		if (self.hasMJD):
			print("Found MJD data")
			print("\tTime (Min/Max) = %f/%f"%(self.MJDs.min(),self.MJDs.max()))

	#Var eq slice
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

	#Radial profile thru cell centers
	def RadialProfileGrid(self):
		self.GetGrid(doVerbose=True)
		#cell corners
		x = self.X [:,:,:]
		y = self.Y [:,:,:]
		z = self.Z [:,:,:]
		#cell centers
		x_c = 0.125*(x[:-1,:-1,:-1]+x[:-1,:-1,1:]+x[:-1,1:,:-1]+x[:-1,1:,1:]+
		x[1:,:-1,:-1]+x[1:,:-1,1:]+x[1:,1:,:-1]+x[1:,1:,1:])
		y_c = 0.125*(y[:-1,:-1,:-1]+y[:-1,:-1,1:]+y[:-1,1:,:-1]+y[:-1,1:,1:]+
		y[1:,:-1,:-1]+y[1:,:-1,1:]+y[1:,1:,:-1]+y[1:,1:,1:])
		z_c = 0.125*(z[:-1,:-1,:-1]+z[:-1,:-1,1:]+z[:-1,1:,:-1]+z[:-1,1:,1:]+
		z[1:,:-1,:-1]+z[1:,:-1,1:]+z[1:,1:,:-1]+z[1:,1:,1:])
		#radius of cell centers
		jR = self.jRad
		kR = self.kRad
		r = np.sqrt(x_c[:,jR,kR]**2.0 + y_c[:,jR,kR]**2.0 + z_c[:,jR,kR]**2.)

		return r	

	#NOT USED merid plane Y=0
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

	#merid plane Y=0 from two halfs
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

	#Grid at 1 AU lat lon
	def iSliceGrid(self):
		#Get Grid
		self.GetGrid(doVerbose=True)

		rxy = np.sqrt(self.X**2 + self.Y**2)
		theta = np.arctan2(rxy,self.Z)
		phi = np.arctan2(self.Y,self.X)

		#theta [theta < 0] += np.pi/2.
		theta += -np.pi/2.
		theta = theta*180./np.pi
		phi [phi < 0] += 2*np.pi
		phi = phi*180./np.pi

		#last i-index == face of the last cell
		lat = theta[-1,::-1,:]
		lon = phi[-1,:,:]
		#these are corners

		return lat, lon

	#Vars at Y=0
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

	#Var at 1 AU
	def iSliceVar(self,vID,sID=None,vScl=None,doVerb=True):
		#Get full 3D variable first
		Q = self.GetVar(vID,sID,vScl,doVerb)

		#cell centered values from the last cell
		Qi = Q[-1,:,:]
                #cell centered values from the first cell
		#Qi = Q[0,:,:]
		#jd_c = self.MJDs[sID]
		#print ('jd_c = ', jd_c)
		return Qi

	#Var along 1D radial line
	def RadialProfileVar(self,vID,sID=None,vScl=None,doVerb=True):
		#Get full 3D variable first
		Q = self.GetVar(vID,sID,vScl,doVerb)

		#set j and k for a radial profile
		jR = self.jRad
		kR = self.kRad
		Nr = self.Ni
              
		Qi = np.zeros(Nr)
        #variable in a cell center
		Qi = Q[:,jR,kR] 
	
		return Qi

	#Radial Profile: Normalized Density
	def RadProfDen(self,s0=0):
		D = self.RadialProfileVar("D", s0)
		r = self.RadialProfileGrid()
		Norm = r**2./r[0]/r[0]
		
		D = D*Norm*self.dScl
		return D

	#Radial Profile: Speed
	def RadProfSpeed(self,s0=0):
		Vx = self.RadialProfileVar("Vx", s0)
		Vy = self.RadialProfileVar("Vy", s0)
		Vz = self.RadialProfileVar("Vz", s0)

		MagV = self.vScl*np.sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
		return MagV

	#Radial Profile: Normalized Flux rho*V*r^2
	def RadProfFlux(self,s0=0):
		D = self.RadialProfileVar("D", s0)
		Vx = self.RadialProfileVar("Vx", s0)
		Vy = self.RadialProfileVar("Vy", s0)
		Vz = self.RadialProfileVar("Vz", s0)
		r = self.RadialProfileGrid()
		
		Norm = r[:]**2./r[0]/r[0]

		Flux = D*Norm*self.dScl*self.vScl*np.sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
		return Flux

	#Speed at 1 AU
	def iSliceMagV(self,s0=0):
		Vx = self.iSliceVar("Vx",s0) #Unscaled
		Vy = self.iSliceVar("Vy",s0) #Unscaled
		Vz = self.iSliceVar("Vz",s0) #Unscaled
		Vi = self.vScl*np.sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
		return Vi

	#Density at 1 AU
	def iSliceD(self,s0=0):
		Di = self.iSliceVar("D",s0) #Unscaled
		Di = Di*self.dScl
		return Di

	#Br at 1 AU
	def iSliceBr(self,s0=0):
		Bx = self.iSliceVar("Bx",s0) #Unscaled
		By = self.iSliceVar("By",s0) #Unscaled
		Bz = self.iSliceVar("Bz",s0) #Unscaled

		self.GetGrid(doVerbose=True)
		x = self.X[-1,:,:]
		y = self.Y[-1,:,:]
		z = self.Z[-1,:,:]
		#centers
		x_c = 0.25*( x[:-1,:-1]+x[:-1,1:]+x[1:,:-1]+x[1:,1:] )
		y_c = 0.25*( y[:-1,:-1]+y[:-1,1:]+y[1:,:-1]+y[1:,1:] )
		z_c = 0.25*( z[:-1,:-1]+z[:-1,1:]+z[1:,:-1]+z[1:,1:] )
		Br = self.bScl*(Bx*x_c + By*y_c + Bz*z_c)/np.sqrt(x_c**2.+y_c**2.+z_c**2.)
		return Br

	#Br at first cell
	def iSliceBrBound(self,s0=0):
		Bx = self.iSliceVar("Bx",s0) #Unscaled
		By = self.iSliceVar("By",s0) #Unscaled
		Bz = self.iSliceVar("Bz",s0) #Unscaled

		self.GetGrid(doVerbose=True)
		x = self.X[0,:,:]
		y = self.Y[0,:,:]
		z = self.Z[0,:,:]
		#centers
		x_c = 0.25*( x[:-1,:-1]+x[:-1,1:]+x[1:,:-1]+x[1:,1:] )
		y_c = 0.25*( y[:-1,:-1]+y[:-1,1:]+y[1:,:-1]+y[1:,1:] )
		z_c = 0.25*( z[:-1,:-1]+z[:-1,1:]+z[1:,:-1]+z[1:,1:] )
		Br = self.bScl*(Bx*x_c + By*y_c + Bz*z_c)/np.sqrt(x_c**2.+y_c**2.+z_c**2.)

		return Br


	#temperature at 1 AU
	def iSliceT(self,s0=0):
		Pi = self.iSliceVar("P",s0) #Unscaled
		Di = self.iSliceVar("D",s0) #Unscaled

		Temp = Pi/Di*self.TScl
		return Temp
		

	#Equatorial speed (in km/s) in eq plane
	def eqMagV(self,s0=0):
		Vx = self.EqSlice("Vx",s0) #Unscaled
		Vy = self.EqSlice("Vy",s0) #Unscaled
		Vz = self.EqSlice("Vz",s0) #Unscaled
		Veq = self.vScl*np.sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
		return Veq

	#Normalized density (D*r*r/21.5/21.5 in cm-3) in eq plane
	def eqNormD (self,s0=0):

		D = self.EqSlice("D",s0) #Unscaled

		Norm = (self.xxc**2.0 + self.yyc**2.0)/self.R0/self.R0
		NormDeq = self.dScl*D*Norm
		return NormDeq

	#Normalized Br (Br*r*r/21.5/21.5) in eq plane
	def eqNormBr (self,s0=0):
		Bx = self.EqSlice("Bx",s0) #Unscaled
		By = self.EqSlice("By",s0) #Unscaled
		Bz = self.EqSlice("Bz",s0) #Unscaled

		Br = (Bx*self.xxc + By*self.yyc)*np.sqrt(self.xxc**2.0 + self.yyc**2.0)/self.R0/self.R0
		
		NormBreq = self.bScl*Br
		return NormBreq

	#Temperature T(r/r0) in eq plane
	def eqTemp (self,s0=0):
		Pres = self.EqSlice("P",s0)
		D = self.EqSlice("D",s0)
		
		#T(r/r0)
		Temp = Pres/D*self.TScl*np.sqrt(self.xxc**2.0 + self.yyc**2.0)/self.R0
		
		return Temp
	
	#Meridional speed (in km/s) in Y=0 plane
	def MerMagV(self,s0=0):
		Vxr, Vxl = self.MeridSlice("Vx",s0) #Unscaled
		Vyr, Vyl = self.MeridSlice("Vy",s0) #Unscaled
		Vzr, Vzl = self.MeridSlice("Vz",s0) #Unscaled
		MagVr = self.vScl*np.sqrt(Vxr**2.0+Vyr**2.0+Vzr**2.0)
		MagVl = self.vScl*np.sqrt(Vxl**2.0+Vyl**2.0+Vzl**2.0)
		return MagVr, MagVl

	#Normalized D in Y=0 plane
	def MerDNrm(self,s0=0):
		xr, zr, xl, zl, r = self.MeridGridHalfs()
		Dr, Dl = self.MeridSlice("D",s0) #Unscaled
		Drn = Dr*self.dScl*r*r/self.R0/self.R0
		Dln = Dl*self.dScl*r*r/self.R0/self.R0
		return Drn, Dln

	#Mormalized Br in Y=0 plane
	def MerBrNrm(self,s0=0):
		xr, zr, xl, zl, r = self.MeridGridHalfs()
		Bxr, Bxl = self.MeridSlice("Bx",s0) #Unscaled
		Bzr, Bzl = self.MeridSlice("Bz",s0) #Unscaled

		#cell centers to calculate Br
		xr_c = 0.25*( xr[:-1,:-1]+xr[:-1,1:]+xr[1:,:-1]+xr[1:,1:] )
		zr_c = 0.25*( zr[:-1,:-1]+zr[:-1,1:]+zr[1:,:-1]+zr[1:,1:] )
		
		xl_c = 0.25*( xl[:-1,:-1]+xl[:-1,1:]+xl[1:,:-1]+xl[1:,1:] )
		zl_c = 0.25*( zl[:-1,:-1]+zl[:-1,1:]+zl[1:,:-1]+zl[1:,1:] )

		#calculating Br
		Br_r = (Bxr*xr_c + Bzr*zr_c)*r*self.bScl/self.R0/self.R0
		Br_l = (Bxl*xl_c + Bzl*zl_c)*r*self.bScl/self.R0/self.R0 
		return Br_r, Br_l

	#Normalized Temp in Y=0 plane 
	def MerTemp(self,s0=0):
		xr, zr, xl, zl, r = self.MeridGridHalfs()

		Pr, Pl = self.MeridSlice("P",s0) #Unscaled
		Dr, Dl = self.MeridSlice("D",s0) #Unscaled

		Tempr = Pr/Dr*self.TScl*r/self.R0
		Templ = Pl/Dl*self.TScl*r/self.R0
		return Tempr, Templ

	#Not used for helio as of now
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

	
