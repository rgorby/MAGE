# Contains the datastructures that include the fully realized lambda distributions
from dataclasses import dataclass
from dataclasses import asdict as dc_asdict
from typing import Optional, List

import kaipy.raiju.lambdautils.AlamParams as aP

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
    """
        Initialize by giving a SpecParam object
        Species will generate alam dist using the params
    """
    params: aP.SpecParams  # Parameters used to generate this instance of Species

    def __post_init__(self):
        # len(n+1), because n in specParams is number of channels, and these are interfaces
        self.alami = self.params.genAlami()

#@dataclass_json
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class AlamData:
    """ 
        Initialize by giving a list of SpecParams
        AlamData will generate all of them and hold as a list of Species
    """
    params: List[aP.SpecParams]

    def __post_init__(self):
        self.species = [Species(p) for p in self.params]

