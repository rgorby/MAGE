import numpy as np
import h5py as h5
from scipy.interpolate import RectBivariateSpline
from kaipy.rcm.wmutils.wmData import wmParams

def genWM(params, useWM=True):

        import os

        fInChorus = 'chorus_polynomial.txt'
        __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

        fInChorus = os.path.join(__location__,fInChorus)

        print("Reading %s"%fInChorus)

        if useWM:
                return genChorus(params,fInChorus)
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

#read parameters of the polynomial fit, Wang+,2023
def readPoly(fIn):
    table = []
    with open(fIn, 'r') as file:
        # Skip the first row
        next(file)
        for line in file:
            row = line.strip().split('\t')[2:-1]  # Discard the first two elements of each row
            row = [float(x) for x in row]  # Convert the strings to float
            rowLen = len(row)
            table.append(np.array(row))
    return (rowLen,np.array(table))

#Chorus polynomial fit for the electron lifetime
def ChorusPoly(Lpoly,Kpoly,polyTgt):

        c0 = polyTgt[0]#Intercept
        c1 = polyTgt[1] #L              
        c2 = polyTgt[2]#log10(E)        
        c3 = polyTgt[3]# L^2            
        c4 = polyTgt[4]#log10(E)^2      
        c5 = polyTgt[5]#L^3             
        c6 = polyTgt[6]#log10(E)^3      
        c7 = polyTgt[7]#log10(E)*L      
        c8 = polyTgt[8]#log10(E)*L^2    
        c9 = polyTgt[9]#log10(E)^2*L

        lenL = len(Lpoly)
        lenK = len(Kpoly)

        tau = np.ones((lenL,lenK))
        # Duplicating the array in columns
        L = np.tile(Lpoly, (lenK, 1)).T

        # Duplicating the array in rows
        K = np.tile(Kpoly, (lenL, 1))

        tau =   c0*tau+\
        c1*L+c2*K+\
        c3*np.power(L,2)+\
        c4*np.power(K,2)+\
        c5*np.power(L,3)+\
        c6*np.power(K,3)+\
        c7*L*K+\
        c8*np.power(L,2)*K+\
        c9*L*np.power(K,2) #in log10(days)

        tau = 10.0**tau*(60.*60.*24.) #in seconds

        return tau

def ReSample(L,MLT,Qp,xMLT):
        Nr,Np = Qp.shape
        #Add ghosts in MLT to handle periodic boundary
        Ng = 2
        Npg = Np+Ng*2
        gMLT = np.arange(0-Ng,24+Ng+1)
        Qpg = np.zeros((Nr,Npg))
        #Set center and then left/right strips
        Qpg[:,2:-2] = Qp
        Qpg[:,1] = Qp[:,-1]
        Qpg[:,0] = Qp[:,-2]
        Qpg[:,-1] = Qp[:,0]
        Qpg[:,-2] = Qp[:,1]

        Q = np.log10(Qpg)
        upQ = RectBivariateSpline(L,gMLT,Q,s=10)

        Qu = upQ(L,xMLT)
        xQp = 10.0**(Qu)
        #Enforce equality at overlap point
        tauP = 0.5*(xQp[:,0]+xQp[:,-1])
        xQp[:, 0] = tauP
        xQp[:,-1] = tauP

        return xQp

def genChorus(params,fInChorus):
        rowLen,paramArray = readPoly(fInChorus)
        polyArray = paramArray.reshape(24,7,rowLen) #Dim MLT: 24, Dim Kp: 7
        lenMLT = 24
        #Kpi
        startValue = 1.0 #in Re 
        endValue = 7.0
        lenKp = 7
        Kpi = np.linspace(startValue, endValue, num=lenKp) 
        #Eki
        startValue = 1.0e-3 #in MeV
        endValue = 2.0  
        lenEk = 155  
        Eki = np.linspace(np.log10(startValue), np.log10(endValue), lenEk) #in log10(MeV)
        #Li
        startValue = 3.0 #in Re 
        endValue = 7.0
        lenL = 41  
        Li = np.linspace(startValue, endValue, num=lenL) 
        #Tau from polynomial fit
        tauP = np.zeros((lenKp,lenMLT,lenL,lenEk))  
        for k in range(0,lenKp): # Kp: 1,2,...,7
            for m in range(0,lenMLT): # MLT: 0,1,...,23
                #print("xMLT,xKp",xMLT,xKp)
                polyKM = polyArray[m,k,:]
                #print('polyTgt',polyTgt)
                tauPolyKM = ChorusPoly(Li,Eki,polyKM)
                tauP[k,m,:,:] = tauPolyKM[:,:]
        #expand MLT from 0-23 to 0-24
        tauE = np.array([np.append(tauP[0,:,:,:],np.array([tauP[0,0,:,:]]),0)])
        for i in range(1,lenKp):
              tauE=np.append(tauE,np.array([np.append(tauP[i,:,:,:],np.array([tauP[i,0,:,:]]),0)]),0)
        tauE = tauE.T
        #Interpolation in the MLT dimesion
        xFac = 4
        lenMLTx = lenMLT*xFac+1 # 97
        MLTi = np.linspace(0,24,lenMLT+1)
        xMLTi = np.linspace(0,24,lenMLTx) 
        tauX = np.zeros((lenEk,lenL,lenMLTx,lenKp))
        for i in range(lenEk):
            for j in range(lenKp):
                Q = tauE[i,:,:,j]
                tauX[i,:,:,j] = ReSample(Li,MLTi,Q,xMLTi)
        Eki = 10.0**Eki #in MeV
        return Kpi,xMLTi,Li,Eki,tauX



        









