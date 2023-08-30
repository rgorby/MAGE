import numpy as np
import h5py as h5
from scipy.interpolate import RectBivariateSpline
from kaipy.rcm.wmutils.wmData import wmParams
# 
def genWM(params, useWM=True):

        import os
      
        fInChorus = 'DWang_chorus_lifetime.h5'
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

        fInChorus = os.path.join(__location__,fInChorus)

        print("Reading %s"%fInChorus)
	
        if useWM:
                return readWM(params,fInChorus)
        else:
                return toyWM(params)
      


# Add wpi-induced electron lifetime model to input file and create an output file
# Writes arrays to file in rcmconfig.h5 format
def genh5(fIn, fOut, inputParams, useWM=True):

        if fIn != fOut: 
               oH5 = h5.File(fOut, 'w')
               iH5 = h5.File(fIn,'r') 
               for Q in iH5.keys():
                     sQ = str(Q)
                     oH5.create_dataset(sQ, data=iH5[sQ])
        else:
               oH5 = h5.File(fOut, 'r+')
        
        if not ('Taui' in oH5.keys()):
               kpi, mlti, li, eki, taui = genWM(inputParams, useWM = useWM)   
               attrs = inputParams.getAttrs()

               oH5.create_dataset('Kpi', data=kpi)
               oH5.create_dataset('MLTi', data=mlti)
               oH5.create_dataset('Li', data=li)
               oH5.create_dataset('Eki', data=eki)
               oH5.create_dataset('Taui', data=taui)
               for key in attrs.keys():
                       oH5.attrs[key] = attrs[key]
        oH5.close()

def ReSample(cL,cK,Qp,s):
        Q = np.log10(Qp)
        SmthQ = RectBivariateSpline(cL,cK,Q.T,s=s)
        Qs = SmthQ(cL,cK).T
        return 10.0**(Qs)

def readWM(params,fInChorus):
        
        # add electron lifetime for the chorus wave loss
        f5 = h5.File(fInChorus, 'r')
        kpi=f5['Kp_1D'][:][0]
        mlti=np.append(f5['MLT_1D'][:][0],24.) # 0, 1,...,24
        li=f5['L_1D'][:][0]
        eki=10.**(f5['E_1D'][:][0]) # in MeV
        taui=(10.**(f5['Tau2_4D'][:]))*24.*3600. # in second, use method 2 data
        #print ("kpi",kpi,"mlti",mlti,"li",li,"eki",eki)
        f5.close()
        nk,nm,nl,ne = taui.shape
        #expand mlt from 0:23 to 0:24
        tauEi = np.array([np.append(taui[0,:,:,:],np.array([taui[0,0,:,:]]),0)])
        for i in range(1,nk):
              tauEi=np.append(tauEi,np.array([np.append(taui[i,:,:,:],np.array([taui[i,0,:,:]]),0)]),0)	
        tauEi = tauEi.T
        #Smooth the expanded tau table  
        ne1,nl1,nm1,nk1 = tauEi.shape
        tauSi = np.zeros((ne1,nl1,nm1,nk1)) #create new smoothed tau
        s0 = 250 #Magic smoothing number
        for m in range(nm1):
            for k in range(nk1):
                MLT0 = mlti[m]
                Kp0  = kpi[k]
                #print("MLT,Kp = %f,%f"%(MLT0,Kp0))
                m0 = np.abs(mlti-MLT0).argmin()
                k0 = np.abs(kpi-Kp0).argmin()
                tauMK = tauEi[:,:,m0,k0]
                tauMKs = ReSample(li,eki*1e3,tauMK,s0)
                tauSi[:,:,m,k] = tauMKs 

        return kpi,mlti,li,eki,tauSi

def toyWM(params):
        nKpi = params.nKp
        nMLTi = params.nMLT
        nLi = params.nL
        nEki = params.nEk
        
        kpi = np.linspace(1,7,nKpi)
        mlti = np.linspace(0,24,nMLTi) #Note the dimension of MLT is 25
        li = np.linspace(3.,7.,nLi)
        eki = np.exp(np.linspace(-3,0.1,nEki)) #in MeV
        taui = kpi[:,None,None,None]*mlti[None,:,None,None]*li[None,None,:,None]*eki[None,None,None,:]
        taui = taui.T
 
        return kpi,mlti,li,eki,taui

 
