import numpy as np
import h5py as h5
from kaipy.rcm.wmutils.wmData import wmParams
# 
def genWM(params, useWMh5=True):

        import os
      
        fIn = 'DWang_chorus_lifetime.h5'
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

        fIn = os.path.join(__location__,fIn)

        print("Reading %s"%(fIn))
	
        if useWMh5:
                return readWMh5(params,fIn)
        else:
                return toyWM(params)


# Add wpi-induced electron loss to rcmconfig.h5
# Writes arrays to file in rcmconfig.h5 format
def genh5(fname, inputParams, useWMh5=True):

        kpi, mlti, li, eki, tau1i, tau2i = genWM(inputParams, useWMh5 = useWMh5)
        attrs = inputParams.getAttrs()

        f5 = h5.File(fname, 'r+')
        f5.create_dataset('Kpi', data=kpi)
        f5.create_dataset('MLTi', data=mlti)
        f5.create_dataset('Li', data=li)
        f5.create_dataset('Eki', data=eki)
        f5.create_dataset('Tau1i', data=tau1i)
        f5.create_dataset('Tau2i', data=tau2i)
        for key in attrs.keys():
                f5.attrs[key] = attrs[key]
        f5.close()

def readWMh5(params,fIn):

        f5 = h5.File(fIn, 'r')
        kpi=f5['Kp_1D'][:][0]
        mlti=np.append(f5['MLT_1D'][:][0],24.)
        li=f5['L_1D'][:][0]
        eki=10.**(f5['E_1D'][:][0]) # in MeV
        tau1i=(10.**(f5['Tau1_4D'][:]))*24.*3600. # in second 
        tau2i=(10.**(f5['Tau2_4D'][:]))*24.*3600.
        #print ("kpi",kpi,"mlti",mlti,"li",li,"eki",eki)
        nk,nm,nl,ne = tau1i.shape
        #expand mlt from 0:23 to 0:24
        tau1ai = np.array([np.append(tau1i[0,:,:,:],np.array([tau1i[0,0,:,:]]),0)])
        tau2ai = np.array([np.append(tau2i[0,:,:,:],np.array([tau2i[0,0,:,:]]),0)])
        for i in range(1,7):
              tau1ai=np.append(tau1ai,np.array([np.append(tau1i[i,:,:,:],np.array([tau1i[i,0,:,:]]),0)]),0)
              tau2ai=np.append(tau2ai,np.array([np.append(tau2i[i,:,:,:],np.array([tau2i[i,0,:,:]]),0)]),0)	
        tau1ai = tau1ai.T
        tau2ai = tau2ai.T
        return kpi,mlti,li,eki,tau1ai,tau2ai


def toyWM(params):
        nKpi = params.nKp
        nMLTi = params.nMLT
        nLi = params.nL
        nEki = params.nEk
        
        kpi = np.linspace(1,7,nKpi)
        mlti = np.linspace(0,24,nMLTi) #Note the dimension of MLT is 25
        li = np.linspace(3.,7.,nLi)
        eki = np.exp(np.linspace(-3,0.1,nEki)) #in MeV
        #print ("kpi",kpi,"mlti",mlti,"li",li,"eki",eki) 
        tau1i = np.zeros((nKpi,nMLTi,nLi,nEki))
        tau2i = np.zeros((nKpi,nMLTi,nLi,nEki)).T 
        tau1i = kpi[:,None,None,None]*mlti[None,:,None,None]*li[None,None,:,None]*eki[None,None,None,:]
        tau1i = tau1i.T
   
        return kpi,mlti,li,eki,tau1i,tau2i
      

 
