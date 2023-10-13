import h5py as h5
from dataclasses import asdict as dc_asdict

import kaipy.kaijson as kj
import kaipy.raiju.lambdautils.AlamData as aD

def saveAlamParams(f5name: str, alamData: aD.AlamData):
    with h5.File(f5name, 'a') as f5:
        print("Adding params used to generate lambda distribution as root attribute")
        f5.attrs['AlamParams'] = [kj.dumps(dc_asdict(p),noIndent=True) for p in alamData.params]


def saveAlamData(f5name: str, alamData: aD.AlamData):
    """ Takes an AlamData object, formats it to raijuconfig.h5 style, and saves it
        Will also call saveAlamParams above if it did write AlamData to file
    """

    with h5.File(f5name, 'a') as f5:

        # Only add species if one doesn't already exist
        if "Species" in f5.keys():
            print("'Species' group already in {}, not writing new one".format(f5name))
            return
        
        # If it doens't exist yet, we add
        print("Adding Species to",f5name)
        gSpec = f5.create_group("Species")

        for i, spec in enumerate(alamData.species):
            gFlav = gSpec.create_group(str(i))

            gFlav.create_dataset('alami', data=spec.alami)
            gFlav['alami'].attrs['units'] = "eV * (Rp/nT)^(2/3)"
            gFlav.attrs['Name'     ] = spec.params.name
            gFlav.attrs['N'        ] = spec.params.n
            gFlav.attrs['flav'     ] = spec.params.flav
            gFlav.attrs['numNuc_p' ] = spec.params.numNuc_p
            gFlav.attrs['numNuc_n' ] = spec.params.numNuc_n
            gFlav.attrs['q'        ] = spec.params.q
            gFlav.attrs['fudge'    ] = spec.params.fudge

    
    # If we successfully added AlamData, save the params as well
    saveAlamParams(f5name, alamData)