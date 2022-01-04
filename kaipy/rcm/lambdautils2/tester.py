import kaipy.rcm.lambdautils2.AlamData as aD
import kaipy.rcm.lambdautils2.AlamParams as aP
import kaipy.rcm.lambdautils2.DistTypes as dT

import kaipy.rcm.lambdautils2.genAlam as genAlam

from dataclasses import asdict as dc_asdict

import h5py as h5
import fileIO

dtWolf = dT.DT_Wolf(p1=2,p2=3)
#dtWolf = getDistTypeObj()
specParams = aP.SpecParams(100,10,10000,dtWolf,2, 1/3)
alamParams = aP.AlamParams(True,[specParams])

spec = genAlam.genSpeciesFromParams(specParams)
alamData = genAlam.genAlamDataFromParams(alamParams)

fileIO.saveRCMConfig(alamData,params=alamParams)