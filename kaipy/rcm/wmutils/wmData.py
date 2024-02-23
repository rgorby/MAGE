import h5py as h5
import numpy as np

class wmParams:

    #All energies in eV
    def __init__(self, dim = 4, nKp = 6, nMLT = 97, nL = 41, nEk = 155):
        self.dim        = dim
        self.nKp        = nKp
        self.nMLT       = nMLT
        self.nL         = nL
        self.nEk        = nEk
    def getAttrs(self):
        return {
                'tauDim': self.dim,
                'nKp': self.nKp,
                'nMLT': self.nMLT,
                'nL': self.nL,
                'nEk': self.nEk,
        }

