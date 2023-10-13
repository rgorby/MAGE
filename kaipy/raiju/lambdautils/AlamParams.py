# Structures to specify a species in raiju. Populate and give to genAlam.py to get what goes into raijuconfig.h5
from dataclasses import dataclass
from dataclasses import asdict as dc_asdict
from typing import Optional, List

# Import other things from this package space
import kaipy.raiju.lambdautils.DistTypes as dT


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


#------
# Actual classes/structures
#------
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class SpecParams:
	""" Defines a single species
		A full species parameter set is defined by the combination of the params listed here AND the ones in whatever DistType is chosen
	"""
	n: int  # Number of channels
	amin: float  # Lower lambda bound for species
	amax: float  # Upper lambda bound for species
	distType: dT.DistType  # DistType params used to generate final lambda distribution
	flav: int  # "Flavor", used to distinguish species types in RCM
			   # 0 = 0-channel plasmasphere, 1 = electrons, 2 = protons
	numNuc_p : int # Number of protons in nucleus
	numNuc_n : int # Number of neutrons in nucleus
	q : int # Net charge of species
	fudge: Optional[float] = 0  # "Fudge factor" loss ratio
	name: Optional[str] = None

	def genAlami(self):  # This will call the given DistType's 'required' function to generate alams based on its rules
		specData = self.distType.genAlami(self)
		return specData