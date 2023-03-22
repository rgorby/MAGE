import os
import subprocess
import numpy as np
import json
import h5py as h5
from dataclasses import asdict as dc_asdict

import kaipy.kaijson as kj
import kaipy.sif.lambdautils.AlamParams as aP
import kaipy.sif.lambdautils.DistTypes as dT

def saveAlamParams(f5name, alamData):
    with h5.File(f5name, 'a') as f5:
        f5.attrs['AlamParams'] = [kj.dumps(dc_asdict(p),noIndent=True) for p in alamData.params]

def saveAlamData(f5name, alamData, doPrint=False):
    """ Takes an AlamData object, formats it to sifconfig.h5 style, and saves it
    """

    with h5.File(f5name, 'a') as f5:

        gSpec = f5.create_group("Species")

        for i, spec in enumerate(alamData.species):
            gFlav = gSpec.create_group(str(i))

            gFlav.create_dataset('alami', data=spec.alami)
            gFlav.attrs['Name' ] = spec.params.name
            gFlav.attrs['N'   ] = spec.params.n
            gFlav.attrs['flav' ] = spec.params.flav
            gFlav.attrs['amu' ] = spec.params.amu
            gFlav.attrs['fudge'] = spec.params.fudge
            #lambdas = np.append(lambdas, spec.alami)
            #flavs  = np.append(flavs , [spec.params.flav  for f in range(spec.params.n)])
            #fudges = np.append(fudges, [spec.params.fudge for f in range(spec.params.n)])
        """
        if doPrint:
            print(lambdas)
            print(flavs)
            print(fudges)

        f5.create_dataset('Species/alami', data=lambdas)
        f5.create_dataset('ikflavc', data=flavs)
        f5.create_dataset('fudgec', data=fudges)
        """