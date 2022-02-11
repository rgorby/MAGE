""" This script is responsible for taking AlamParam objects and turning it into lambdas and associated attributes for AlamData objects
"""
import numpy as np

import kaipy.rcm.lambdautils.AlamParams as aP
import kaipy.rcm.lambdautils.AlamData as aD
import kaipy.rcm.lambdautils.DistTypes as dT

#------
# Main functions
#------

def getAlamMinMax(alams):
    amin = np.array([])
    amax = np.array([])
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

    return amin.tolist(), amax.tolist()

def genSpeciesFromParams(specParams):
    """ Takes a SpecParams object and generates a new Species object
    """
    alams = specParams.genAlams()
    n = len(alams)  # !!!Bad workaround. Since genAlams might change n, it should set it in its specParams container itself
    amins,amaxs = getAlamMinMax(alams)
    flav = specParams.flav
    fudge = specParams.fudge
    name = specParams.name
    return aD.Species(n, alams, amins, amaxs, flav, fudge, params=specParams, name=name)

def genAlamDataFromParams(alamParams):
    doUsePsphere = alamParams.doUsePsphere
    specList = [genSpeciesFromParams(sP) for sP in alamParams.specParams]
    if doUsePsphere:
        specList.insert(0,genPsphereSpecies())
    return aD.AlamData(doUsePsphere, specList, params=alamParams)


#------
# Helpers
#------

def genPsphereSpecies():
    params = dT.DT_Manual(name='Plasmasphere')
    return aD.Species(1, [0], [0], [0], 1, 0, params=params, name='Plasmasphere')
