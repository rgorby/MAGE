from dataclasses import dataclass
#from dataclasses_json import dataclass_json
from dataclasses import asdict as dc_asdict
from typing import Optional, List

#Import other things from this package space
import kaipy.rcm.lambdautils.DistTypes as dT

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


# Defines a single species
@dataclass
class SpecParams:
	""" Defines a single species
		A full species parameter set is defined by both the params listed here AND the ones in whatever DistType is chosen
	"""
	n: int  # Number of channels
	amin: float  # Lower lambda bound for species
	amax: float  # Upper lambda bound for species
	distType: dT.DistType  # DistType params used to generate final lambad distribution
	flav: int  # "Flavor", used to distinguish species types in RCM
			   # 1 = electrons, 2 = protons
	fudge: Optional[float] = 0  # "Fudge factor" loss ratio
	name: Optional[str] = None

	def genAlams(self):  # This will call the given DistType's 'required' function to generate alams based on its rules
		specData = self.distType.genAlamsFromSpecies(self)
		return specData

#@dataclass_json
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class AlamParams:
	doUsePsphere: bool  # Whether or not the resulting dataset will include a zero-energy plasmasphere channel
	specParams: List[SpecParams]  # List of all specParams to be included in dataset

	# These can help to determine some things for higher-level lambda generators
	#tiote: Optional[float] = None  # Ratio of ion temperature to electron temperature
	#ktMax: Optional[float] = None  # Energy in eV
	emine: Optional[float] = None  # Electron min energy in eV
	eminp: Optional[float] = None  # Proton min energy in eV
	emaxe: Optional[float] = None  # Electron max energy in eV
	emaxp: Optional[float] = None  # Proton max energy in eV
	L_kt : Optional[float] = None  # L value that e/p min/maxes correspond to


