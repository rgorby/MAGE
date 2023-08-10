import numpy as np
import h5py as h5
from kaipy.rcm.wmutils.wmData import wmParams
# 
def genWM(params, useWM=True):

        import os
      
        fInChorus = 'DWang_chorus_lifetime.h5'
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

        fInChorus = os.path.join(__location__,fInChorus)

        fInTDS = 'tauTDS.txt'
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

        fInTDS = os.path.join(__location__,fInTDS)

        print("Reading %s and %s"%(fInChorus,fInTDS))
	
        if useWM:
                return readWM(params,fInChorus,fInTDS)
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
        
        if not ('Tau1i' in oH5.keys()):
               kpi, mlti, li, eki, tau1i, tau2i, ekTDSi, tauTDSi = genWM(inputParams, useWM = useWM)
               attrs = inputParams.getAttrs()

               oH5.create_dataset('Kpi', data=kpi)
               oH5.create_dataset('MLTi', data=mlti)
               oH5.create_dataset('Li', data=li)
               oH5.create_dataset('Eki', data=eki)
               oH5.create_dataset('Tau1i', data=tau1i)
               oH5.create_dataset('Tau2i', data=tau2i)
               oH5.create_dataset('EkTDSi', data=ekTDSi)
               oH5.create_dataset('TauTDSi', data=tauTDSi)
               for key in attrs.keys():
                       oH5.attrs[key] = attrs[key]
        oH5.close()

def readWM(params,fInChorus,fInTDS):
        
        # add electron lifetime for the chorus wave loss
        f5 = h5.File(fInChorus, 'r')
        kpi=f5['Kp_1D'][:][0]
        mlti=np.append(f5['MLT_1D'][:][0],24.)
        li=f5['L_1D'][:][0]
        eki=10.**(f5['E_1D'][:][0]) # in MeV
        print('shape of eki:',eki.shape)
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
        f5.close()

        #add electron lifetime for the Time Domain Structure loss
        tdmArrays=np.loadtxt(fInTDS)
        #print(tdm_arrays)
        ekTDSi = tdmArrays[:,0].T/1.e6 #in MeV
        print ('shape of ekiTDS',ekTDSi.shape)
        tauTDSi = tdmArrays[:,2].T*24.*3600. #in second, read in the electron lifetime against TDS with Ew= 4mV/m
        print ('tauTDSi[0]',tauTDSi[0])           
        return kpi,mlti,li,eki,tau1ai,tau2ai,ekTDSi,tauTDSi



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
      

 
