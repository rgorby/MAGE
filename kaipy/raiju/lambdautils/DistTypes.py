# Defines distribution types explicitly, so we can always take a config file and know exactly how the distribution was generated
from dataclasses import dataclass
from dataclasses import asdict as dc_asdict
from typing import Optional, List


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
# Parameters needed to determine lambda distribution
#------

@dataclass
class DistType:  # Empty class just so we can force the type in dataclasses below
    name: str = "Empty"

#------
# Specific implementations of DistType
#------

#@dataclass_json
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class DT_Manual(DistType):
    """
        If you want to completely ignore all this generator stuff, you could use this
        and do whatever you want. But if you were doing all that anyways, its probably
        just easier to manually edit the h5 dataset directly
    """
    def __post_init__(self):
        if self.name == "Empty": self.name = "Manual"

#@dataclass_json
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class DT_Single(DistType):
    """
        Just a single channel, pretty boring
    """
    def __post_init__(self):
        self.name = "Single"

    def genAlami(self, sP):
        return [sP.amin, sP.amax]

#@dataclass_json
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class DT_Wolf(DistType):
    """ Lambda channel spacing based on Wolf's notes 
            (ask Anthony Sciola or Frank Toffoletto for a copy)
        With the addition that there can be 2 p values for the start and end, and pStar transitions between them
        kmin, kmax let you set what part of the full distribution you actually generate
    """
    p1: float = None
    p2: float = None
    kmin: int = 0
    kmax: int = -1

    def __post_init__(self):
        self.name = "Wolf"

    def genAlami(self, sP):
        """
            Takes a filled out SpecParams and returns the lambda interfaces
        """
        # Get/set needed variables
        kmin = self.kmin
        kmax = sP.n+1 if self.kmax == -1  else self.kmax # Add 1 because sP.n is # channels and we are generating lambda interfaces
        amin = sP.amin
        amax = sP.amax


        # Do math
        alams = []
        for k in range(sP.n+1):
            kfrac = (k-kmin)/(kmax-kmin)  # How far through the channel range are we
            pstar = (1-kfrac)*self.p1 + kfrac*self.p2
            lammax = amax-amin
            lam = lammax*((k - kmin + 0.5)/(kmax-kmin + 0.5))**pstar + amin
            alams.append(lam)

        return alams


# TODO: There's another DistType in the previous RCM generator, should port over