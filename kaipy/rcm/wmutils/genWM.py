import numpy as np
import h5py as h5
from kaipy.rcm.wmutils.wmData import wmParams
# 
def genWM(params, useWMh5=True):
	
        if useWMh5:
                return readWMh5(params,'DWang_chorus_lifetime.h5')
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
        kpi=f5['Kp_1D']
        mlti=f5['MLT_1D']
        li=f5['L_1D']
        eki=f5['E_1D']
        tau1i=f5['Tau1_4D']  
        tau2i=f5['Tau2_4D']

	#check with params

        return kpi,mlti,li,eki,tau1i,tau2i


def toyWM(params):
        nKpi = params.nKp
        nMLTi = params.nMLT
        nLi = params.nL
        nEki = params.nEk
        
        kpi = np.linspace(1,7,nKpi)
        mlti = np.linspace(0,23,nMLTi)
        li = np.linspace(3.,7.,nLi)
        eki = np.exp(np.linspace(-3,0.1,nEki))*1.e6 #in eV
        
        tau1i = np.zeros((nKpi,nMLTi,nLi,nEki))
        tau2i = np.zeros((nKpi,nMLTi,nLi,nEki))  
        for nk in range(nKpi):
            for nm in range(nMLTi):
                for nl in range(nLi):
                    for ne in range(nEki):
                         tau1i[nk,nm,nl,ne] = kpi[nk]*mlti[nm]*li[nl]*eki[ne] 
                         tau2i[nk,nm,nl,ne] = tau1i[nk,nm,nl,ne]
        
        return kpi,mlti,li,eki,tau1i,tau2i
      

 
