from dataclasses import dataclass
#from dataclasses_json import dataclass_json
from dataclasses import asdict as dc_asdict
from typing import Optional, List

import kaipy.rcm.lambdautils.AlamParams as aP

# dataclasses_json isn't a default package. Since its only used for reading, don't want to make it a requirement for everyone
try:
    from dataclasses_json import dataclass_json
    dataclasses_json_module_imported = True
except ModuleNotFoundError:
    dataclass_json = None
    dataclasses_json_module_imported = False

def conditional_decorator(dec, dataclasses_json_module_imported):
    def decorator(func):
        if not dataclasses_json_module_imported:
            # Return the function unchanged, not decorated.
            return func
        return dec(func)
    return decorator

@dataclass
class Species:
	n: int  # Number of channels
	alams: List[float]  # Lambda channel values
	amins: List[float]  # Lower lambda bounds for species
	amaxs: List[float]  # Upper lambda bound for species
	flav: int  # "Flavor", used to distinguish species types in RCM
			   # 1 = electrons, 2 = protons
	fudge: Optional[float] = 0  # "Fudge factor" loss ratio
	params: Optional[aP.SpecParams] = None  # Parameters used to generate this instance of Species
	name: Optional[str] = None

#@dataclass_json
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class AlamData:
	""" Main class that most things will interact with
	"""
	doUsePsphere: bool  # Whether or not this dataset includes a zero-energy plasmasphere channel
	specs: List[Species]  # List of Species objects
	params: Optional[aP.AlamParams] = None  # Parameters used to generate this instance of AlamData

