# Gamera post-processing routines
# Get data from serial/MPI gamera output
import glob
import numpy as np

# Object to use to pull data from HDF5 structure (serial or mpi)

# Initialize,
# gIn = gampp.GameraPipe(fdir,ftag)
# Set grid from data
# gIn.GetGrid()
# V = gIn.GetVar("D",stepnum)
# doFast=True skips various data scraping
idStr = "_0000_0000_0000.gam.h5"
idStr_old = "_0000_0000_0000_000000000000.h5"


class GameraPipe(object):

        # Initialize GP object
        # fdir = directory to h5 files
        # ftag = stub of h5 files
        def __init__(self, fdir, ftag, doFast=False, doVerbose=True):

                self.fdir = fdir
                self.ftag = ftag

                # Here just create variables
                # -----------
                self.isOld = False
                # Ranks
                self.isMPI = False
                self.Nr = 0
                self.Ri = 1; self.Rj = 1; self.Rk = 1
                # Global/local cells
                self.is2D = False
                self.Ni = 0; self.Nj = 0; self.Nk = 0
                self.dNi = [] ; self.dNj = [] ; self.dNk = []

                # Variables/slices
                self.Nv = 0; self.Nt = 0; self.Nv0 = 0
                self.T = []; self.vIDs = []; self.v0IDs = []
                self.s0 = 0; self.sFin = 0

                # Grids
                self.X = []; self.Y = []; self.Z = []

                self.doFast = doFast
                self.UnitsID = "NONE"

                # Stubs for MJD stuff
                self.hasMJD = False
                self.MJDs = []

                # Scrape data from directory
                self.doVerbose = True
                self.OpenPipe(doVerbose)

        def OpenPipe(self, doVerbose=True):
                import kaipy.kaiH5 as kh5
                if (doVerbose):
                        print("Opening pipe: %s : %s" % (self.fdir, self.ftag))

                # Test for serial (either old or new naming convention)
                fOld = "%s/%s.h5" % (self.fdir, self.ftag)
                fNew = "%s/%s.gam.h5" % (self.fdir, self.ftag)
                print("fNew ",fNew)

                if (len(glob.glob(fOld)) == 1):
                        # Found old-style serial
                        self.isMPI = False
                        if (doVerbose):
                                print("Found serial database")
                        f0 = fOld
                elif (len(glob.glob(fNew)) == 1):
                        # Found new-style serial
                        self.isMPI = False
                        self.ftag = self.ftag + ".gam"  # Add .gam to tag
                        if (doVerbose):
                                print("Found serial database")
                        f0 = fNew
                else:
                        print("%s not found, looking for MPI database" % (fNew))
                        self.isMPI = True
                        self.isOld = False
                        sStr = "%s/%s_*%s" % (self.fdir, self.ftag, idStr)

                        fIns = glob.glob(sStr)
                        if (len(fIns) > 1):
                                print("This shouldn't happen, bailing ...")
                        if (len(fIns) == 0):
                                print("No new MPI database found, checking old ...")
                                self.isMPI = True
                                self.isOld = True
                                sStr = "%s/%s_*%s" % (self.fdir, self.ftag, idStr_old)

                                fIns = glob.glob(sStr)
                                if (len(fIns) > 1):
                                        print("This shouldn't happen, bailing ...")
                                if (len(fIns) == 0):
                                        print("No old MPI database found, out of options, bailing ...")
                                        quit()
                
                        f0 = fIns[0]
                        print("fIns ",fIns)
                        Ns = [int(s) for s in f0.split('_') if s.isdigit()]

                        self.Ri = Ns[0]
                        self.Rj = Ns[1]
                        self.Rk = Ns[2]
                        self.Nr = self.Ri*self.Rj*self.Rk
                        if (doVerbose):
                                print("\tFound %d = (%d,%d,%d) ranks" %
                                      (self.Nr, self.Ri, self.Rj, self.Rk))
                        self.dNi = np.zeros(self.Ri,dtype=int)
                        self.dNj = np.zeros(self.Rj,dtype=int)
                        self.dNk = np.zeros(self.Rk,dtype=int)

                # In either case, f0 is defined.  use it to get per file stuff
                self.Nt, sids = kh5.cntSteps(f0)
                if (self.doFast):
                        self.T = np.zeros(self.Nt)
                else:
                        self.T = kh5.getTs(f0, sids)

                self.s0 = sids.min()
                self.sFin = sids.max()
                if (doVerbose):
                        print("Found %d timesteps\n\tTime = [%f,%f]" % (
                            self.Nt, self.T.min(), self.T.max()))
                        print("\tSteps = [%d,%d]" % (sids.min(), sids.max()))

                # Get MJD if present
#                MJDs = kh5.getTs(f0, sids, "MJD", -np.inf)
#                if (MJDs.max() > 0):
#                        self.hasMJD = True
 #                       self.MJDs = MJDs

                # Get grid stuff

                Dims = kh5.getDims(f0)
                print("DIMS ",Dims)
                Nd=len(Dims)
                """
                self.dNi=Dims[0]-1
                self.dNj=Dims[1]-1
                if (Nd > 2):
                        self.dNk=Dims[2]-1

                else:
                        self.dNk=0
                        self.Rk=1
                self.Ni=self.dNi*self.Ri
                self.Nj=self.dNj*self.Rj
                self.Nk=self.dNk*self.Rk
"""
                self.Ni = 0
                for i in range(self.Ri):
                        if (self.isOld):
                                joff = i*self.Rj 
                                fi = "%s/%s_%04d_%04d_%04d_%04d_0000_0000_%012d.h5" % (self.fdir, self.ftag, self.Ri, self.Rj, self.Rk,i,joff)
                        else:
                                fi = "%s/%s_%04d_%04d_%04d_%04d_0000_0000.gam.h5" % (self.fdir, self.ftag, self.Ri, self.Rj, self.Rk,i)
                        print(" fi is ",fi)
                        dims = kh5.getDims(fi)
                        self.dNi[i] = dims[0]-1
                        self.Ni = self.Ni + self.dNi[i]
                        print(" Ni, dNi ",self.Ni,self.dNi)
                self.Nj = 0
                for j in range(self.Rj):
                        if (self.isOld):
                                fj = "%s/%s_%04d_%04d_%04d_0000_%04d_0000_%012d.h5" % (self.fdir, self.ftag, self.Ri, self.Rj, self.Rk,j,j)
                        else:
                                fj = "%s/%s_%04d_%04d_%04d_0000_%04d_0000.gam.h5" % (self.fdir, self.ftag, self.Ri, self.Rj, self.Rk,j)
                        print(fj)
                        dims = kh5.getDims(fj)
                        self.dNj[j] = dims[1]-1
                        self.Nj = self.Nj + self.dNj[j]
                        print(" Nj, dNj ",self.Nj,self.dNj)
                if (Nd >2 ):
                        for k in range(self.Rk):
                                if (self.isOld):
                                        koff = self.Ri*self.Rj*k
                                        fk = "%s/%s_%04d_%04d_%04d_0000_0000_%04d_%012d.h5" % (self.fdir, self.ftag, self.Ri, self.Rj, self.Rk,k,koff)
                                else:
                                        fk = "%s/%s_%04d_%04d_%04d_0000_0000_%04d.gam.h5" % (self.fdir, self.ftag, self.Ri, self.Rj, self.Rk,k)
                                print(fk)
                                dims = kh5.getDims(fk)
                                self.dNk[k] = dims[2]-1
                                self.Nk = self.Nk + self.dNk[k]
                                print(" Nk, dNk ",self.Nk,self.dNk)

                # Variables
                self.v0IDs=kh5.getRootVars(f0)
                self.vIDs=kh5.getVars(f0, sids.min())
                self.Nv0=len(self.v0IDs)
                self.Nv=len(self.vIDs)
                if (Nd > 2):
                        self.is2D=False
                else:
                        self.is2D=True
                if (doVerbose):
                        if (Nd > 2):
                                nCells=self.Ni*self.Nj*self.Nk
                        else:
                                nCells=self.Ni*self.Nj
                        print("Grid size = (%d,%d,%d)" % (self.Ni, self.Nj, self.Nk))
                        print("\tCells = %e" % (nCells))
                        print("Variables (Root/Step) = (%d,%d)" % (self.Nv0, self.Nv))
                        print("\tRoot: %s" % (self.v0IDs))
                        print("\tStep: %s" % (self.vIDs))

                self.SetUnits(f0)
                if (doVerbose):
                        print("Units Type = %s" % (self.UnitsID))
                        print("Pulling grid ...")
                self.GetGrid(doVerbose)

        def SetUnits(self, f0):
                import h5py
                with h5py.File(f0, 'r') as hf:
                        uID=hf.attrs.get("UnitsID", "CODE")
                if (not isinstance(uID, str)):
                        self.UnitsID=uID.decode('utf-8')

        # Chunks go from 0,Ri-1 (Python iteration)
        def ChunkName(self, i, j, k):
                print( "ChunkName %s,%d,%d,%d "%(self.isOld,i,j,k))
                if (self.isMPI):
                        fIn="blah"
                        n=j + i*self.Rj + k*self.Ri*self.Rj
                        if (self.isOld):
                                fRn="%04d_%04d_%04d_%04d_%04d_%04d_%012d" % (self.Ri, self.Rj, self.Rk, i, j, k,n)
                                fIn=self.fdir+"/"+self.ftag+"_"+fRn+".h5"
                        else:
                                fRn="%04d_%04d_%04d_%04d_%04d_%04d" % (self.Ri, self.Rj, self.Rk, i, j, k)
                                fIn=self.fdir+"/"+self.ftag+"_"+fRn+".gam.h5"
                else:
                        fIn=self.fdir+"/"+self.ftag+".h5"
                return fIn
                print(fIn)
        def GetGrid(self, doVerbose):
                import kaipy.kaiH5 as kh5
                if (self.is2D):
                        self.X=np.zeros((self.Ni+1, self.Nj+1))
                        self.Y=np.zeros((self.Ni+1, self.Nj+1))
                else:
                        self.X=np.zeros((self.Ni+1, self.Nj+1, self.Nk+1))
                        self.Y=np.zeros((self.Ni+1, self.Nj+1, self.Nk+1))
                        self.Z=np.zeros((self.Ni+1, self.Nj+1, self.Nk+1))
                if (doVerbose):
                        print("Deli = " , self.dNi)
                        print("Delj = " , self.dNj)
                        print("Delk = " , self.dNk)
                iS = 0
                for i in range(self.Ri):
                        iE=iS+self.dNi[i]
                        jS = 0
                        for j in range(self.Rj):
                                jE=jS+self.dNj[j]
                                kS = 0
                                for k in range(self.Rk):
                                        kE=kS+self.dNk[k]
                                        print("Bounds = (%d,%d,%d,%d,%d,%d)" % (iS, iE, jS, jE, kS, kE))
                                        #print(fIn)
                                        # print("Bounds = (%d,%d,%d,%d,%d,%d)"%fIn.GetDims)
                                        fIn=self.ChunkName(i, j, k)
                                        if (self.is2D):
                                                self.X[iS:iE+1, jS:jE+1]=kh5.PullVar(fIn, "X")
                                                self.Y[iS:iE+1, jS:jE+1]=kh5.PullVar(fIn, "Y")
                                        else:
                                                self.X[iS:iE+1, jS:jE+1, kS:kE+1]=kh5.PullVar(fIn, "X")
                                                self.Y[iS:iE+1, jS:jE+1, kS:kE+1]=kh5.PullVar(fIn, "Y")
                                                self.Z[iS:iE+1, jS:jE+1, kS:kE+1]=kh5.PullVar(fIn, "Z")

                                        kS=k*self.dNk[k]                                        
                                jS=jS + self.dNj[j]
                        iS=iS + self.dNi[i]
        # Get 3D variable "vID" from Step# sID
        def GetVar(self, vID, sID, vScl=None, doVerb=True):
                import kaipy.kaiH5 as kh5

                if (doVerb):
                        if (sID is None):
                                print("Reading %s/%s" % (self.ftag, vID))
                        else:
                                print("Reading %s/Step#%d/%s" % (self.ftag, sID, vID))

                if (self.is2D):
                        V=np.zeros((self.Ni, self.Nj))
                else:
                        V=np.zeros((self.Ni, self.Nj, self.Nk))

                for i in range(self.Ri):
                        for j in range(self.Rj):
                                for k in range(self.Rk):
                                        iS=i*self.dNi
                                        jS=j*self.dNj
                                        kS=k*self.dNk
                                        iE=iS+self.dNi
                                        jE=jS+self.dNj
                                        kE=kS+self.dNk
                                        # print("Bounds = (%d,%d,%d,%d,%d,%d)"%(iS,iE,jS,jE,kS,kE))
                                        fIn=self.ChunkName(i, j, k)
                                        if (self.is2D):
                                                V[iS:iE, jS:jE]=kh5.PullVar(fIn, vID, sID)

                                        else:
                                                V[iS:iE, jS:jE, kS:kE]=kh5.PullVar(fIn, vID, sID)
                if (vScl is not None):
                        V=vScl*V
                return V

        # Get variable slice of constant i,j,k
        # Directions = idir/jdir/kdir strings
        # Indexing = (1,Nijk)
        # FIXME: Currently pulling whole 3D array and then slicing, lazy
        def GetSlice(self, vID, sID, ijkdir='idir', n=1, vScl=None, doVerb=True):
                sDirs=['IDIR', 'JDIR', 'KDIR']
                ijkdir=ijkdir.upper()

                if (not (ijkdir in sDirs)):
                        print("Invalid slice direction, defaulting to I")
                        ijkdir='IDIR'
                if (sID is None):
                        cStr="Reading %s/%s" % (self.ftag, vID)
                else:
                        cStr="Reading %s/Step#%d/%s" % (self.ftag, sID, vID)
                V=self.GetVar(vID, sID, vScl, doVerb=False)
                # Now slice
                np=n-1  # Convert from Fortran to Python indexing
                if (ijkdir == "IDIR"):
                        Vs=V[np, :, :]
                        sStr="(%d,:,:)" % (n)
                elif (ijkdir == "JDIR"):
                        Vs=V[:, np, :]
                        sStr="(:,%d,:)" % (n)
                elif (ijkdir == "KDIR"):
                        Vs=V[:, :, np]
                        sStr="(:,:,%d)" % (n)
                if (doVerb):
                        print(cStr+sStr)
                return Vs

        # Wrappers for root variables (or just set sID to "None")
        def GetRootVar(self, vID, vScl=None, doVerb=True):
                V=self.GetVar(vID, sID=None, vScl=vScl, doVerb=doVerb)

        def GetRootSlice(self, vID, ijkdir='idir', n=1, vScl=None, doVerb=True):
                Vs=self.GetSlice(vID, sID=None, ijkdir=ijkdir, n=n, vScl=vScl, doVerb=doVerb)
                return Vs
