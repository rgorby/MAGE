#!/usr/bin/env python

#python wsa2TDgamera.py 

import os,sys,glob
import scipy
from scipy import interpolate
from scipy.optimize import newton_krylov,anderson
import h5py
import numpy as np
import matplotlib.pyplot as plt

import time
import kaipy.gamhelio.wsa2TDgamera.params as params
import kaipy.gamhelio.lib.wsa as wsa
import kaipy.gamhelio.lib.poisson as poisson
import kaipy.gamera.gamGrids as gg

#plotting function for debug
def plot(wsa_file, var_wsa, var_wsa_rolled):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
 
    #fig=plt.figure(figsize=(16,12))
    fig=plt.figure(figsize=(10,6.))
    gs = gridspec.GridSpec(2,1,height_ratios=[20,1])
    ax1 = fig.add_subplot(gs[0,0])
    axc = fig.add_subplot(gs[1,0])
    p1=ax1.pcolormesh(var_wsa_rolled[::-1,:]*1.e5,cmap='RdBu_r',vmin=-150.,vmax=150.)
    ax1.contour(var_wsa_rolled[::-1,:],[0.],colors='white')

    plt.colorbar(p1,cax=axc, orientation = 'horizontal').set_label('Br [nT]')

    ax1.set_xlim((0,var_wsa.shape[1]))
    ax1.set_ylim((0,var_wsa.shape[0]))
    ax1.set_aspect("equal")

    ##in rotating system of coordinates
    #ax2 = plt.subplot(212,sharex=ax1)
    ##p2=ax2.pcolormesh(var_wsa_rolled)
    #p2=ax2.pcolormesh(var_wsa_rolled[::-1,:],cmap='RdBu_r',vmin=var_wsa_rolled.min(),vmax=-var_wsa_rolled.min())
    #plt.colorbar(p2,ax=ax2).set_label('V')
    #ax2.set_xlim((0,var_wsa_rolled.shape[1]))
    #ax2.set_ylim((0,var_wsa_rolled.shape[0]))

    fig.suptitle(wsaFile)
    plt.savefig(wsaFile[:-4]+'png')


# [EP] function to plot boundary conditions in rotating frame to make a movie
def plotBc(wsa_file, phi, theta, var1, var2, var3, var4):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    fig=plt.figure(figsize=(12,9))
    gs = gridspec.GridSpec(2,2,wspace=0.2, hspace =0.1)

    phi = phi*180./np.pi
    theta = (np.pi/2.-theta)*180./np.pi

    var4 = var4/1.e6 #temp in MK
    
    ax1 = plt.subplot(gs[0,0], aspect='equal')
    p1=ax1.pcolormesh(phi, theta, var1.T, shading = 'auto', cmap = 'rainbow', vmin = 300, vmax = 850)
    plt.colorbar(p1,ax=ax1,aspect = 15, orientation = 'horizontal').set_label(r'$V_r$, km/s')

    ax2 = plt.subplot(gs[0,1],sharex=ax1, aspect='equal')
    p2=ax2.pcolormesh(phi, theta,var2.T, shading = 'auto', cmap = 'RdBu_r', vmin = -150, vmax = 150)
    plt.colorbar(p2,ax=ax2,aspect = 15, orientation = 'horizontal').set_label(r'$B_r, nT$')

    ax3 = plt.subplot(gs[1,0],sharex=ax1, aspect='equal')
    p3=ax3.pcolormesh(phi, theta,var3.T, shading = 'auto', cmap = 'copper_r', vmin = 300, vmax = 1200)
    plt.colorbar(p3,ax=ax3,aspect = 15, orientation = 'horizontal').set_label(r'$Rho, cm^{-3}$')

    ax4 = plt.subplot(gs[1,1],sharex=ax1, aspect='equal')
    p4=ax4.pcolormesh(phi, theta, var4.T,shading = 'auto', cmap = 'copper', vmin = 0.5, vmax = 2.5)
    plt.colorbar(p4,ax=ax4,aspect = 15, orientation = 'horizontal').set_label('Temperature, K')
    
    date = wsa_file.split('/')[-1][4:12]
    year = date[0:4] 
    month = date[4:6]
    day = date[6:8]

    plt.suptitle(year + ':' + month + ':' + day, y=0.85)
    plt.savefig(wsaFile[:-5]+'_bc.png', bbox_inches='tight')

#----------- PARSE ARGUMENTS ---------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('ConfigFileName',help='The name of the configuration file to use',default='startup.config')
args = parser.parse_args()
#----------- PARSE ARGUMENTS ---------#

# Read params from config file
prm = params.params(args.ConfigFileName)
(ni,nj,nk) = (prm.Ni,prm.Nj,prm.Nk)
Ng = prm.NO2

#grid parameters
tMin = prm.tMin
tMax = prm.tMax
Rin = prm.Rin
Rout = prm.Rout
Ni = prm.Ni
Nj = prm.Nj
Nk = prm.Nk

#----------GENERATE HELIO GRID------

print("Generating gamera-helio grid ...")

X3,Y3,Z3 = gg.GenKSph(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)

#to generate non-uniform grid for GL cme (more fine in region 0.1-0.3 AU) 
#X3,Y3,Z3 = gg.GenKSphNonUGL(Ni=Ni,Nj=Nj,Nk=Nk,Rin=Rin,Rout=Rout,tMin=tMin,tMax=tMax)
gg.WriteGrid(X3,Y3,Z3,fOut=os.path.join(prm.GridDir,prm.gameraGridFile))

print("Gamera-helio grid ready!")

#----------GENERATE HELIO GRID------


# [EP] sorted list of WSA files
wsaFiles = sorted(glob.glob(os.path.join(prm.adaptdir,prm.adaptWildCard)))

print(wsaFiles)

# [EP] electric fields on edges
#+2 in j directions, two ghost cells at start and end
et_save = np.zeros( (nj+2,nk+1) )
ep_save = np.zeros( (nj+1+2,nk) )

#Normalization
Vnorm = 1.e5 #cm/s => km/s
Bnorm = 1.e-5 #Gs => nT
mp = 1.67e-24
kblts = 1.38e-16

#open innerbcTD.h5 for ouput
with h5py.File(os.path.join(prm.IbcDir,prm.gameraIbcFile),'w') as hf:

#[EP] going through the list of WSA files
    for (fcount,wsaFile) in enumerate(wsaFiles):
        #print(fcount)
    	############### WSA STUFF #####################
        isFirstFile = (wsaFile == wsaFiles[0]) 
       	#[EP] reading WSA file
        jd_c,phi_wsa_v,theta_wsa_v,phi_wsa_c,theta_wsa_c,bi_wsa,v_wsa,n_wsa,T_wsa = wsa.read(wsaFile,prm.densTempInfile,prm.normalized, verbose = isFirstFile)
        #bi_wsa in Gs CGS units
	    #v_wsa in cm/s
	    #n_wsa in g/cm-3
	    #T_wsa in K


        #convert julian date from wsa fits into modified julian date
        mjd_c = jd_c - 2400000.5

        if isFirstFile:
            #take JD from the first wsa file
            jd0 = jd_c

            # GAMERA GRID
            # read GAMERA grid from innerbc.h5
            
            print ('reading heliogrid.h5 ...')
            f = h5py.File(os.path.join(prm.GridDir,prm.gameraGridFile), 'r')
            #Nphi, Nth, Nr = np.shape(f['X'])
            #corners
            x = f['X'][:]
            y = f['Y'][:]
            z = f['Z'][:]

            #centers
            xc = 0.125*(f['X'][:-1,:-1,:-1]+f['X'][:-1,:-1,1:]+f['X'][:-1,1:,:-1]+f['X'][:-1,1:,1:]+
                 f['X'][1:,:-1,:-1]+f['X'][1:,:-1,1:]+f['X'][1:,1:,:-1]+f['X'][1:,1:,1:])
            yc = 0.125*(f['Y'][:-1,:-1,:-1]+f['Y'][:-1,:-1,1:]+f['Y'][:-1,1:,:-1]+f['Y'][:-1,1:,1:]+
                 f['Y'][1:,:-1,:-1]+f['Y'][1:,:-1,1:]+f['Y'][1:,1:,:-1]+f['Y'][1:,1:,1:])
            zc = 0.125*(f['Z'][:-1,:-1,:-1]+f['Z'][:-1,:-1,1:]+f['Z'][:-1,1:,:-1]+f['Z'][:-1,1:,1:]+
                 f['Z'][1:,:-1,:-1]+f['Z'][1:,:-1,1:]+f['Z'][1:,1:,:-1]+f['Z'][1:,1:,1:]) 

            #radius of inner boundary. Index order [k,j,i]
            R0 = np.sqrt(x[0,0,Ng]**2+y[0,0,Ng]**2+z[0,0,Ng]**2)

            #[EP for testing]
            #cell corners including ghost cells
            r = np.sqrt(x[:]**2+y[:]**2+z[:]**2)
            rxy = np.sqrt(x[:]**2+y[:]**2)
          
            # remove the ghosts from angular dimensions (corners)
            P = np.arctan2(y[Ng:-Ng,Ng:-Ng,:],x[Ng:-Ng,Ng:-Ng,:])
            P [ P < 0] += 2*np.pi
            T = np.arccos(z[Ng:-Ng,Ng:-Ng,:]/r[Ng:-Ng,Ng:-Ng,:])

            #grid for output into innerbc.h5
            P_out = P[:,:,0:Ng+1]
            T_out = T[:,:,0:Ng+1]
            R_out = r[Ng:-Ng,Ng:-Ng,0:Ng+1]
            print ("shapes of output phi and theta ", P_out.shape, T_out.shape, R_out.shape)

            #centers spherical grid excluding ghosts in angular directions
            
            #Rc = np.sqrt(xc[Ng:-Ng, Ng:-Ng,:]**2 + yc[Ng:-Ng, Ng:-Ng,:]**2 + zc[Ng:-Ng, Ng:-Ng,:]**2)
            #Pc = np.arctan2(yc[Ng:-Ng, Ng:-Ng,:], xc[Ng:-Ng, Ng:-Ng,:])
            #Tc = np.arccos(zc[Ng:-Ng,Ng:-Ng,:]/Rc)

            #include one extra cell in j direction at start and end
            Pg = Ng-1
            Rc = np.sqrt(xc[Ng:-Ng, Pg:-Pg,:]**2 + yc[Ng:-Ng, Pg:-Pg,:]**2 + zc[Ng:-Ng, Pg:-Pg,:]**2)
            Pc = np.arctan2(yc[Ng:-Ng, Pg:-Pg,:], xc[Ng:-Ng, Pg:-Pg,:])
            Tc = np.arccos(zc[Ng:-Ng,Pg:-Pg,:]/Rc)
            

            Pc [Pc < 0] += 2*np.pi
            #GAMERA grid centers at the inner boundary, 1D array
            phi = Pc[:,0,0]
            theta = Tc[0,:,0]

            #debug
            #print (phi)
            #print (theta)

            # what exactly does this do???
            pois = poisson.poisson(theta,phi)

        #time from the 1st wsa map in seconds
        time_sec = (jd_c - jd0)*24.*60.*60.

        omega=2*np.pi/prm.Tsolar*(25.38/27.27)
        #shift of wsa maps needed if wsa solutions are provided in inertial frame (folder UPDATED)
        #if wsa solutions are provided in rotating carrington frame (folder CARR), no need to shift.
        #shift phi coordinate in wsa data according to the shift of the wsa map relative to the first one 
        #wsa maps move to the right with cadence 1 day

        phi_prime=(phi_wsa_c-omega*prm.adaptCadence*fcount)%(2*np.pi)
        #looking for index of shift
        if np.where(np.ediff1d(phi_prime)<0)[0].size!=0: #for the first map size =0, for other maps size=1
            ind0=np.where(np.ediff1d(phi_prime)<0)[0][0]+1
            #print 'ind = ', ind0
        else:
            ind0=0 # this is for the first map
        
        #shifting phi_prime to the left
        phi_prime=np.roll(phi_prime,-ind0)
        bi_wsa_rolled=np.roll(bi_wsa,-ind0,axis=1)
        v_wsa_rolled=np.roll(v_wsa,-ind0,axis=1)
        n_wsa_rolled=np.roll(n_wsa,-ind0,axis=1)
        T_wsa_rolled=np.roll(T_wsa,-ind0,axis=1)

        #plot br from original wsa map (top plot) and shifted to the origin map(bottom plot)
        #changes in time in the bottom plot are purely due to time-dependent variations of B_r (rotation is eliminted) 
        
        plot(wsaFile, bi_wsa, bi_wsa_rolled)
        ##plot(wsaFile, v_wsa, v_wsa_rolled)

        ###INTERPOLATION OF ROLLED WSA MAPS TO GAMERA GRID phi-theta####

        # bivariate spline approximation over a rectangular mesh
        fbi = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,bi_wsa_rolled.T,kx=1,ky=1)  
        # interpolation to Gamera grid
        br = fbi(phi,theta)

        #Next Slava used SMOOTHING for br for the paper, we do not need it for now

        fv = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,v_wsa_rolled.T,kx=1,ky=1)  
        vr = fv(phi,theta)

        f = interpolate.RectBivariateSpline(phi_wsa_c,theta_wsa_c,n_wsa_rolled.T,kx=1,ky=1)  
        rho = f(phi,theta)

        #not interpolating temperature, calculating sound speed cs
        #assuming uniform total pressure Rho_max*k*T0 = p+Br^2/8pi

        #TODO: Check Temp calculation
        T0 = 0.9e6 
        Rho0 = 1100.*mp       #density in the HCS
        #cs = np.sqrt(prm.gamma/rho*(rho.max()*1.38e-16*T0/1.67e-24-br**2/8/np.pi))
        Temp = mp/rho/kblts*(Rho0*kblts*T0/mp-br**2/8./np.pi/2.)

        # Poisson solver after interpolation onto GAMERA grid
        if fcount>0:
            print ('fcount = ', fcount)
            #right-hand side of laplacian equation
            pois.setRHS( (br-br_save).T) #after transponding it becomes (nj,nk)
            guess=np.zeros_like(br.T)
            #electric field potential psi: (Laplacian(Psi) = dB_r/dt)
            Psi = newton_krylov(pois.residual,guess, method='lgmres',verbose=True,iter=100)#,f_rtol=1.e-6) #iter=100
            print('Residual: %g' % abs(pois.residual(Psi)).max())

            print ('Psi.shape = ', Psi.shape) # (nj, nk) =(128, 256)
            #Psi is defined in cell centers

            #calculate electric field componenet
            #suffix  _a denotes that this is adapt field
            #E_theta = dPsi/dphi/sin(theta)
            #E_phi = -dPsi/dtheta/r

            et_a = np.zeros( (Psi.shape[0],Psi.shape[1]+1) ) #(128, 257)
            et_a[:,1:-1] = np.diff(Psi,axis=1)/np.diff(phi) #except first and last cells in k
            et_a[:,0] = (Psi[:,0] - Psi[:,-1])/(phi[0]-phi[-1]+2*np.pi) #k=0
            et_a[:,-1]=et_a[:,0] #k=Nk+1
            et_a /= np.sin(theta[:,None]) 
            print ('E_theta.shape = ', et_a.shape) 
            """
            note, here we assume theta constant along phi and same theta on the
            boundary and in the center of the GAMERA cell - Elena: we should probably fix that 
            """
            ep_a=np.zeros((Psi.shape[0]+1,Psi.shape[1])) #(129, 256)
            ep_a[1:-1,:] = -np.diff(Psi,axis=0)/np.diff(theta)[:,None] #except first and last cells in j
            #for j=0 and j=N_j we set nearby values
            ep_a[0,:]=ep_a[1,:]   # used to set these to zero, but more appropriate to repeat from next theta, since Ephi does not depend on theta at the pole.
            ep_a[-1,:]=ep_a[-2,:] 
            print ('E_phi.shape = ', ep_a.shape)
            #[EP]: two lines above are from LFM where theta went from pole to pole


            # I do not understand all that business with interpolation of electric fields in time (see adapt2lfm.py)
            # Convert to CGS. FIX ME!!! UNITS HARD CODED
            et_a*= prm.Rin*prm.scale/prm.adaptCadence/24./3600.
            ep_a*= prm.Rin*prm.scale/prm.adaptCadence/24./3600.

            et_save = et_a
            ep_save = ep_a
            #[EP] for debug
            dbr = br - br_save
            #et_save and ep_save are defined at times of adapt and on cell edges

        br_save = br

        """
        After we obtained et_save ep_save at cell edges we calculate B_theta and B_phi
        at cell centers and faces
        """
        vrt = vr.T #(nj,nk) in cell centers
        bp_a = np.zeros_like(vrt) #B_phi in cell centers
        bt_a = np.zeros_like(vrt) #B_theta in cell centers
        bp_kface_a = np.zeros( (vrt.shape[0],vrt.shape[1]+1) ) #(nj,nk+1)
        bt_jface_a = np.zeros( (vrt.shape[0]+1,vrt.shape[1]) ) #(nj+1,nk)
        vrt_kface  = np.zeros( (vrt.shape[0],vrt.shape[1]+1) )
        vrt_jface  = np.zeros( (vrt.shape[0]+1,vrt.shape[1]) )

        if fcount >0:
            # B_phi and B_theta defined at cell centers
            bp_a = 0.5*(et_save[:,:-1]+et_save[:,1:])/vrt
            bt_a = -0.5*(ep_save[:-1,:]+ep_save[1:,:])/vrt

            # the above are at cell centers, also need at the
            # corresponding faces, see below
            
            # First interpolate velocity to faces
            vrt_kface[:,1:-1] = 0.5*(vrt[:,:-1]+vrt[:,1:]); vrt_kface[:,0] = 0.5*(vrt[:,-1]+vrt[:,0]); vrt_kface[:,-1] = vrt_kface[:,0]
            vrt_jface[1:-1,:] = 0.5*(vrt[1:,:]+vrt[:-1,:]) ; vrt_jface[0,:]=vrt[1,:].mean(); vrt_jface[-1,:]=vrt[-2,:].mean(); 
            
            #B_phi and B_theta at faces
            bp_kface_a = et_save/vrt_kface
            bt_jface_a = -ep_save/vrt_jface

        #transponse again to agree with GAMERA indexing nk,nj,ni
        # Note, these are defined at cell centers on the boundary (at rmin)
        bp_a = bp_a.T  
        bt_a = bt_a.T
        #in kaiju we do not to save B-components at cell centers, so we do not need bp_a and bt_a

        # at faces; change shapes to match order in gamera nk, nj
        bt_jface_a = bt_jface_a.T 
        bp_kface_a = bp_kface_a.T

        et_save = et_save.T
        ep_save = ep_save.T


        # Scale inside ghost region
        #print(rho.shape)
        (vr,rho,Temp,br,bp_kface_a,bt_jface_a,et_save,ep_save) = [np.dstack(prm.NO2*[var]) for var in (vr,rho,Temp,br,bp_kface_a,bt_jface_a,et_save,ep_save)]
        rho*=(R0/Rc[0,0,:Ng])**2
        Temp*=(R0/Rc[0,0,:Ng])
        br*=(R0/Rc[0,0,:Ng])**2
        bp_kface_a*=(R0/Rc[0,0,:Ng])
        et_save*=(R0/Rc[0,0,:Ng])

        #tangential velocities are set to zero
        #vp = zeros_like(vr)
        #vt = zeros_like(vr)

        #print vr.shape, rho.shape, cs.shape, br.shape, bt_jface_a.shape, bp_kface_a.shape
        #print et_save.shape, ep_save.shape

        #Agreement. For innerbcTD.h5 the output units are V[km/s], Rho[cm-3], T[K], B[nT]
        #v_wsa /= Vnorm
        #n_wsa /= mp
        #bi_wsa /= Bnorm 

        print (wsaFile)
        #removing two bounding cells in theta and normalizing
        vrp = vr[:,1:-1,:]/Vnorm
        vp = np.zeros_like(vrp)/Vnorm
        vt = np.zeros_like(vrp)/Vnorm
        rhop = rho[:,1:-1,:]/mp
        Tempp = Temp[:,1:-1,:]
        brp = br[:,1:-1,:]/Bnorm
        bt_jface_a_p = bt_jface_a[:,1:-1,:]/Bnorm

        bp_kface_a_p = bp_kface_a[:,1:-1,:]/Bnorm
        et_save_p = et_save[:,1:-1,:]
        ep_save_p = ep_save[:,1:-1,:]

        #print vrp.shape, rhop.shape, csp.shape, brp.shape, bt_jface_a_p.shape, bp_kface_a_p.shape
        #print et_save_p.shape, ep_save_p.shape
        
	#V in cm/s B in Gs n in gcm-3
        print (fcount, time_sec, mjd_c)

        if prm.dumpBC:
            if fcount == 0:
                #write out phi and th coords of corners at inner boundary grid
                hf.create_dataset("X", data=P_out)
                hf.create_dataset("Y", data=T_out)
                hf.create_dataset("Z", data=R_out)
            grname = "Step#"+str(fcount)
            grp = hf.create_group(grname)
            grp.attrs.create("time", time_sec)
            grp.attrs.create("MJD", mjd_c)
            grp.create_dataset("vr",data=vrp) #cc
            grp.create_dataset("vp",data=vp) #cc !zeros
            grp.create_dataset("vt",data=vt) #cc !zeros
            #hf.create_dataset("vr_kface",data=vr_kface) #kface
            grp.create_dataset("rho",data=rhop) #cc
            grp.create_dataset("T",data=Tempp) #cc
            grp.create_dataset("br",data=brp) #cc
            #hf.create_dataset("br_kface",data=br_kface) #kface 
            #hf.create_dataset("bp",data=bp_a) #cc
            #hf.create_dataset("bt",data=bt_a) #cc
            grp.create_dataset("bt_jface",data=bt_jface_a_p) #jface
            grp.create_dataset("bp_kface",data=bp_kface_a_p) #kface
            grp.create_dataset("et",data=et_save_p) #k-edges
            grp.create_dataset("ep",data=ep_save_p) #j-edges

        plotBc(wsaFile,phi, theta[1:-1], vrp[:,:,Ng-1], brp[:,:,Ng-1], rhop[:,:,Ng-1], Tempp[:,:,Ng-1])
        

        # [EP] test if calculated tengential electric fields give Br from wsa
        if fcount > 30:
            print ('Elena debug')
            print (fcount)

            dphi = phi[2]-phi[1]
            dtheta = theta[2]-theta[1]

            dbrp = dbr[:,1:-1]

            
            #cell edges dphi and dtheta at inner boundary face
            dlp = dphi*rxy[Ng:-Ng-1,Ng:-Ng,Ng] #(256,129) dlp change with nj
            #dlp = dphi*r[Ng:-Ng-1,Ng:-Ng,Ng]*sin(T[:,:])
            dlt = dtheta*R0 #dlt is same for all cells

            et_use = et_save_p.T #(257,128)
            ep_use = ep_save_p.T #(256,129)

            #rotE of the cell face
            circE = zeros((nk,nj))

            #circE1 = - (ep_use[:,:-1]*dlp + et_use[1:,:]*dlt - ep_use[:,1:]*dlp - et_use[:-1,:]*dlt)

            for k in range(256):
                for j in range(128):
                    circE[k,j] = - (ep_use[k,j+1]*dlp[k,j+1] - ep_use[k,j]*dlp[k,j] + et_use[k,j]*dlt - et_use[k+1,j]*dlt)

            dt = prm.adaptCadence*24.*3600.
            dbrdt = dlp[:,:-1]*dlt*prm.scale*dbrp/dt
           
            resid = dbrdt - circE
            
            fig1 = plt.figure(); plt.pcolormesh(np.log10(np.abs(resid.T))); plt.colorbar()
            fig1.suptitle(wsaFile)
            plt.savefig(wsaFile[:-5]+'_testFL.png')
         


        #if fcount==2:
         #   sys.exit("passed two wsa files")
            





    



        
