import h5py as h5
import numpy as np
import hashlib

class AlamParams:

    #All energies in eV
    def __init__(self, distType = 'wolf', 
                    num_e      = 39   , num_p     = 120,
                    alamMin_e  = -1.0 , alamMin_p = 10, 
                    ktMax      = 15000, L_kt      = 10, tiote=4.0, 
                    p1         = 3.0  , p2        = 1.0,
                    addPsphere = True):
        self.distType     = distType
        self.num_e        = num_e
        self.num_p        = num_p
        self.aMin_e       = alamMin_e
        self.aMin_p       = alamMin_p
        self.ktMax        = ktMax
        self.L_kt         = L_kt
        self.tiote        = tiote
        self.p1           = p1
        self.p2           = p2
        self.doAddPsphere = addPsphere

    def getAttrs(self):
        return {
                'hash': self.getHash(),
                'distType': self.distType,
                'num_e': self.num_e,
                'num_p': self.num_p,
                'aMin_e': self.aMin_e,
                'aMin_p': self.aMin_p,
                'ktMax': self.ktMax,
                'L_kt': self.L_kt,
                'tiote': self.tiote,
                'p1': self.p1,
                'p2': self.p2,
                'doAddPsphere': self.doAddPsphere
        }
    # Generate an identifier for the current settings
    def getHash(self):
        l = [self.num_e, self.num_p, self.aMin_e, self.aMin_p, self.ktMax, 
                    self.L_kt, self.tiote, self.p1, self.p2, self.doAddPsphere]
        s = '-'.join(["{:1.2f}".format(i) for i in l]).encode('utf-8')
        m = hashlib.md5()
        m.update(s)
        return m.hexdigest()

class AlamData:

    def __init__(self, filename=None, doUsePsphere=False, doReadNow=True, alamdict=None, alamParams=None):
        self.filename = filename
        self.doUsePsphere = doUsePsphere
        
        if doReadNow and filename is not None:
            self.alams, self.amins, self.amaxs = self.readFile(filename)
        
        if alamdict is not None:
            self.buildFromData(alamdict)
        
        if alamParams is not None:
            self.params = alamParams
        else:
            self.params = AlamParams()

    def readFile(self, fname):
        if 'h5' in fname:
            return self.readrcmconfig(fname)

    def readrcmconfig(self, fname):
        f5 = h5.File(fname, 'r')

        alamdict = {}

        #Initialize
        numSpec = max(f5['ikflavc'])
        for i in range(1,numSpec+1):
            alamdict['spec'+str(i)] = np.array([])
        
        #Populate
        for i in range(len(f5['alamc'])):
            specStr = 'spec' + str(f5['ikflavc'][i])
            alamdict[specStr] = np.append(alamdict[specStr], f5['alamc'][i])

        #If doUsePsphere == False, remove 0-channel from electron species
        if not self.doUsePsphere:
            for k in alamdict.keys():
                if alamdict[k][1] < 0:  # Find channel with negative alam, i.e. electron channel
                    alamdict[k] = alamdict[k][1:]
                    break

        #print(alamdict)
        #alamdict['spec2'][0] = 5.645

        amindict, amaxdict = self.getAlamMinMax(alamdict)
        
        return alamdict, amindict, amaxdict


    def buildFromData(self, alamdict):
        if self.doUsePsphere == False and alamdict['spec1'][0] == 0:
            alamdict['spec1'] = alamdict['spec1'][1:]
        self.alams = alamdict
        self.amins, self.amaxs = self.getAlamMinMax(self.alams)

#Helper functions, broken out so others can use them too
    def getAlamMinMax(self, alamdict):
        amindict = {}
        amaxdict = {}
        for k in alamdict.keys():
            amin = np.array([])
            amax = np.array([])
            alams = alamdict[k]
            for i in range(len(alams)):
                if i == 0:
                    amax = np.append(amax, 0.5*(np.abs(alams[i])+np.abs(alams[i+1])))
                    amin = np.append(amin, 0)
                elif i == len(alams)-1:
                    amax = np.append(amax, 1.5*np.abs(alams[i]) - 0.5*np.abs(alams[i-1]))
                    amin = np.append(amin, 0.5*(np.abs(alams[i])+np.abs(alams[i-1])))
                else:
                    amax = np.append(amax, 0.5*(np.abs(alams[i])+np.abs(alams[i+1])))
                    amin = np.append(amin, 0.5*(np.abs(alams[i])+np.abs(alams[i-1])))
            amindict[k] = amin
            amaxdict[k] = amax

        return amindict, amaxdict