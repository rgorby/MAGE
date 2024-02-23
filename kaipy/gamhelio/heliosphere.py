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
    def __init__(self,fdir,ftag,doFast=False,uID="Inner",doParallel=False,nWorkers=4):

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
        GameraPipe.__init__(self,fdir,ftag,doFast=doFast,doParallel=doParallel,nWorkers=nWorkers)

        #inner boundary distance
        self.R0 = np.sqrt(self.X[0, 0, 0]**2 + self.Y[0, 0, 0]**2 + self.Z[0, 0, 0]**2)

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
    
    #Var theta slice
    def jSlice(self,vID,sID=None,vScl=None,doEq=True,doVerb=True,jidx=-1):
        #Get full 3D variable first
        Q = self.GetVar(vID,sID,vScl,doVerb)

        if(jidx == -1):
            Nj2 = self.Nj//2 
        else:
            Nj2 = jidx
        #above and below the j plane
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

    #merid plane from kidx from two halfs
    def MeridGridHalfs(self,kidx=None,phi=None):
        self.GetGrid(doVerbose=True)

        if(kidx):
            Nk1 = kidx
            Nk2 = self.Nk//2 + Nk1
        elif(phi):
            Nk1 = int(phi/(2*np.pi)*self.Nk)
            Nk2 = self.Nk//2 + Nk1
        else: 
            Nk1 = 0
            Nk2 = self.Nk//2

        #looking from -Y to XZ plane
        xright = self.X[:,:,Nk1] #corners
        yright = self.Y[:,:,Nk1] #corners
        zright = self.Z[:,:,Nk1] #corners

        xleft = self.X[:,:,Nk2]
        yleft = self.Y[:,:,Nk2]
        zleft = self.Z[:,:,Nk2]
        rleft = -np.sqrt(xleft**2. + yleft**2)

        xright_c = 0.25*( xright[:-1,:-1]+xright[:-1,1:]+xright[1:,:-1]+xright[1:,1:] )
        yright_c = 0.25*( yright[:-1,:-1]+yright[:-1,1:]+yright[1:,:-1]+yright[1:,1:] )
        zright_c = 0.25*( zright[:-1,:-1]+zright[:-1,1:]+zright[1:,:-1]+zright[1:,1:] )
        r = np.sqrt(xright_c**2 + zright_c**2 + yright_c**2)

        #centers: right plane, left plane, radius
        return xright, yright, zright, xleft, yleft, zleft, r

    #Grid at 1 AU lat lon
    def iSliceGrid(self,idx=-1):
        #Get Grid
        self.GetGrid(doVerbose=True)

        rxy = np.sqrt(self.X**2 + self.Y**2)
        theta = np.arctan2(rxy,self.Z)
        phi = np.arctan2(self.Y,self.X)

        #theta [theta < 0] += np.pi/2.
        theta = theta - np.pi/2
        theta = theta*180./np.pi
        phi [phi < 0] += 2*np.pi
        phi = phi*180./np.pi

        #last i-index == face of the last cell
        lat = theta[idx,::-1,:]
        lon = phi[idx,:,:]
        #these are corners

        return lat, lon

    #Vars at Y=0
    def MeridSlice(self,vID,sID=None,vScl=None,doVerb=True,indx=(None,None)):
        #Get full 3D variable first
        Q = self.GetVar(vID,sID,vScl,doVerb)
        
        Nk2 = self.Nk//2
        kidx, phi = indx
        
        if(kidx):
            Nk1 = kidx
            Nk2 = self.Nk//2 + Nk1
            Np = Nk1 - 1
        elif(phi):
            Nk1 = int(phi/(2*np.pi)*self.Nk)
            Nk2 = self.Nk//2 + Nk1
            Np =  Nk1 - 1
        else: 
            Nk1 = 0
            Nk2 = self.Nk//2
            Np = self.Nk
        
        #Nr = self.Ni
        #Nt = 2*self.Nj
        #XZ meridional slice (k=0) of var 
        #Qj = np.zeros((Nr,Nt))
        
        Qright = 0.5*( Q[:,:,Nk1] + Q[:,:,Np-1] ) 
        Qleft  = 0.5*( Q[:,:,Nk2-1] + Q[:,:,Nk2] )
        #print (Qright.shape, Qleft.shape)
        #Qj = np.hstack( (Qright, Qleft[:,::-1]) ) #reverse in j
        #print (Qj.shape)
        return Qright, Qleft

    #Var at 1 AU
    def iSliceVar(self,vID,sID=None,vScl=None,doVerb=True,idx=-1):
        #Get full 3D variable first
        Q = self.GetVar(vID,sID,vScl,doVerb)

        #cell centered values from the last cell
        Qi = Q[idx,:,:]
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
    def iSliceMagV(self,s0=0,idx=-1):
        Vx = self.iSliceVar("Vx",s0,idx=idx) #Unscaled
        Vy = self.iSliceVar("Vy",s0,idx=idx) #Unscaled
        Vz = self.iSliceVar("Vz",s0,idx=idx) #Unscaled
        Vi = self.vScl*np.sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
        return Vi

    #Density at 1 AU
    def iSliceD(self,s0=0,idx=-1):
        Di = self.iSliceVar("D",s0,idx=idx) #Unscaled
        Di = Di*self.dScl
        return Di

    #Br at 1 AU
    def iSliceBr(self,s0=0,idx=-1):
        Bx = self.iSliceVar("Bx",s0,idx=idx) #Unscaled
        By = self.iSliceVar("By",s0,idx=idx) #Unscaled
        Bz = self.iSliceVar("Bz",s0,idx=idx) #Unscaled

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

    #Bx at 1AU
    def iSliceBx(self,s0=0,idx=-1):
        Bx = self.iSliceVar("Bx",s0,idx=idx) #Unscaled

        self.GetGrid(doVerbose=True)
        x = self.X[-1,:,:]

        #centers
        x_c = 0.25*( x[:-1,:-1]+x[:-1,1:]+x[1:,:-1]+x[1:,1:] )

        BxScl = self.bScl*Bx*x_c
        return BxScl

    #By at 1AU
    def iSliceBy(self,s0=0,idx=-1):
        By = self.iSliceVar("By",s0,idx=idx) #Unscaled

        self.GetGrid(doVerbose=True)
        y = self.Y[-1,:,:]
        #centers
        y_c = 0.25*( y[:-1,:-1]+y[:-1,1:]+y[1:,:-1]+y[1:,1:] )
        ByScl = self.bScl*By*y_c
        return ByScl
    
    #Bz at 1AU
    def iSliceBz(self,s0=0,idx=-1):
        Bz = self.iSliceVar("Bz",s0,idx=idx) #Unscaled

        self.GetGrid(doVerbose=True)
        z = self.Z[-1,:,:]
        #centers

        z_c = 0.25*( z[:-1,:-1]+z[:-1,1:]+z[1:,:-1]+z[1:,1:] )
        BzScl = self.bScl*Bz*z_c
        return BzScl
    
    #Br at first cell
    def iSliceBrBound(self,s0=0,idx=-1):
        Bx = self.iSliceVar("Bx",s0,idx=idx) #Unscaled
        By = self.iSliceVar("By",s0,idx=idx) #Unscaled
        Bz = self.iSliceVar("Bz",s0,idx=idx) #Unscaled

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
    def iSliceT(self,s0=0,idx=-1):
        Pi = self.iSliceVar("P",s0,idx=idx) #Unscaled
        Di = self.iSliceVar("D",s0,idx=idx) #Unscaled

        Temp = Pi/Di*self.TScl
        return Temp
        
    #Equatorial speed (in km/s) in eq plane
    def eqMagV(self,s0=0):
        Vx = self.EqSlice("Vx",s0) #Unscaled
        Vy = self.EqSlice("Vy",s0) #Unscaled
        Vz = self.EqSlice("Vz",s0) #Unscaled
        Veq = self.vScl*np.sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
        return Veq
    
    #Equatorial speed (in km/s) in eq plane
    def jMagV(self,s0=0,jidx=-1):
        Vx = self.jSlice("Vx",s0,jidx=jidx) #Unscaled
        Vy = self.jSlice("Vy",s0,jidx=jidx) #Unscaled
        Vz = self.jSlice("Vz",s0,jidx=jidx) #Unscaled
        Veq = self.vScl*np.sqrt(Vx**2.0+Vy**2.0+Vz**2.0)
        return Veq

    def eqNormD (self, s0=0):
        """Compute the normalized number density in the equatorial plane.

        Compute the normalized number density in the equatorial plane.
        
        The number density is normalized by the factor (r/r0)**2, where r0 is
        the radius of the inner edge of the grid (should be 21.5 Rsun).

        Parameters
        ----------
        self : GamsphPipe
            This object
        s0 : int
            Simulation step number to fetch

        Returns
        -------
        NormDeq : np.ndarray, shape same as self.xxc
            Normalized number density in equatorial plane

        Raises
        ------
        None
        """
        # Fetch the unscaled data.
        D = self.EqSlice("D", s0)

        # Compute the normalization factor for each data point.
        Norm = (self.xxc**2.0 + self.yyc**2.0)/self.R0/self.R0

        # Convert the number density from code units (dimensionless) to
        # physical units (cm**-3), then normalize the data using the scale
        # factor.
        NormDeq = self.dScl*D*Norm

        # Return the data.
        return NormDeq
    

    #Normalized density (D*r*r/21.5/21.5 in cm-3) in eq plane
    def jNormD (self,s0=0,jidx=-1):

        D = self.jSlice("D",s0,jidx=jidx) #Unscaled

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
    
    #Normalized Br (Br*r*r/21.5/21.5) in eq plane
    def jNormBr (self,s0=0,jidx=-1):
        Bx = self.jSlice("Bx",s0,jidx=jidx) #Unscaled
        By = self.jSlice("By",s0,jidx=jidx) #Unscaled
        Bz = self.jSlice("Bz",s0,jidx=jidx) #Unscaled

        Br = (Bx*self.xxc + By*self.yyc)*np.sqrt(self.xxc**2.0 + self.yyc**2.0)/self.R0/self.R0
        
        NormBreq = self.bScl*Br
        return NormBreq
    
    #Normalized Br (Br*r*r/21.5/21.5) in eq plane
    def eqBx (self,s0=0):
        Bx = self.EqSlice("Bx",s0) #Unscaled

        BxScl = (Bx*self.xxc + Bx*self.yyc)*np.sqrt(self.xxc**2.0 + self.yyc**2.0)/self.R0/self.R0
        
        BxScl = self.bScl*BxScl
        return BxScl

    #Normalized Br (Br*r*r/21.5/21.5) in eq plane
    def eqBy (self,s0=0):
        By = self.EqSlice("By",s0) #Unscaled

        ByScl = (By*self.xxc + By*self.yyc)*np.sqrt(self.xxc**2.0 + self.yyc**2.0)/self.R0/self.R0
        
        ByScl= self.bScl*ByScl
        return ByScl

    #Normalized Br (Br*r*r/21.5/21.5) in eq plane
    def eqBz (self,s0=0):
        Bz = self.EqSlice("Bz",s0) #Unscaled

        BzScl = Bz
        
        BzScl = self.bScl*BzScl
        return BzScl


    #Temperature T(r/r0) in eq plane
    def eqTemp (self,s0=0):
        Pres = self.EqSlice("P",s0)
        D = self.EqSlice("D",s0)
        
        #T(r/r0)
        Temp = Pres/D*self.TScl*np.sqrt(self.xxc**2.0 + self.yyc**2.0)/self.R0
        
        return Temp

    #Temperature T(r/r0) in eq plane
    def jTemp (self,s0=0,jidx=-1):
        Pres = self.jSlice("P",s0,jidx=jidx)
        D = self.jSlice("D",s0,jidx=jidx)
        
        #T(r/r0)
        Temp = Pres/D*self.TScl*np.sqrt(self.xxc**2.0 + self.yyc**2.0)/self.R0
        
        return Temp
    
    #Meridional speed (in km/s) in Y=indx plane
    def MerMagV(self,s0=0,indx=(None,None)):
        Vxr, Vxl = self.MeridSlice("Vx",s0,indx=indx) #Unscaled
        Vyr, Vyl = self.MeridSlice("Vy",s0,indx=indx) #Unscaled
        Vzr, Vzl = self.MeridSlice("Vz",s0,indx=indx) #Unscaled
        MagVr = self.vScl*np.sqrt(Vxr**2.0+Vyr**2.0+Vzr**2.0)
        MagVl = self.vScl*np.sqrt(Vxl**2.0+Vyl**2.0+Vzl**2.0)
        return MagVr, MagVl

    #Normalized D in Y=indx plane
    def MerDNrm(self,s0=0,indx=(None,None)):
        xr, yr, zr, xl, yl, zl, r = self.MeridGridHalfs(*indx)
        Dr, Dl = self.MeridSlice("D",s0,indx=indx) #Unscaled
        Drn = Dr*self.dScl*r*r/self.R0/self.R0
        Dln = Dl*self.dScl*r*r/self.R0/self.R0
        return Drn, Dln

    #Mormalized Br in Y=indx plane
    def MerBrNrm(self,s0=0,indx=(None,None)):
        xr, yr, zr, xl, yl, zl, r = self.MeridGridHalfs(*indx)
        #rxyr = np.sqrt(xr**2. + yr**2)
        #rxyl = -np.sqrt(xl**2. + yl**2)
        Bxr, Bxl = self.MeridSlice("Bx",s0,indx=indx) #Unscaled
        Byr, Byl = self.MeridSlice("By",s0,indx=indx) #Unscaled
        Bzr, Bzl = self.MeridSlice("Bz",s0,indx=indx) #Unscaled

        #cell centers to calculate Br
        xr_c = 0.25*( xr[:-1,:-1]+xr[:-1,1:]+xr[1:,:-1]+xr[1:,1:] )
        yr_c = 0.25*( yr[:-1,:-1]+yr[:-1,1:]+yr[1:,:-1]+yr[1:,1:] )
        zr_c = 0.25*( zr[:-1,:-1]+zr[:-1,1:]+zr[1:,:-1]+zr[1:,1:] )
        
        xl_c = 0.25*( xl[:-1,:-1]+xl[:-1,1:]+xl[1:,:-1]+xl[1:,1:] )
        yl_c = 0.25*( yl[:-1,:-1]+yl[:-1,1:]+yl[1:,:-1]+yl[1:,1:] )
        zl_c = 0.25*( zl[:-1,:-1]+zl[:-1,1:]+zl[1:,:-1]+zl[1:,1:] )

        #calculating Br
        Br_r = (Bxr*xr_c + Byr*yr_c + Bzr*zr_c)*r*self.bScl/self.R0/self.R0
        Br_l = (Bxl*xl_c + Byl*yl_c + Bzl*zl_c)*r*self.bScl/self.R0/self.R0 
        return Br_r, Br_l

    #Normalized Temp in Y=indx plane 
    def MerTemp(self,s0=0,indx=(None,None)):
        xr, yr, zr, xl, yl, zl, r = self.MeridGridHalfs(*indx)

        Pr, Pl = self.MeridSlice("P",s0,indx=indx) #Unscaled
        Dr, Dl = self.MeridSlice("D",s0,indx=indx) #Unscaled

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
            Nm = int( (t-T0)/60.0 ) #Minutes, integer
            Hr = Nm/60
            Min = np.mod(Nm,60)
            Sec = np.mod(int(t),60)

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

    
