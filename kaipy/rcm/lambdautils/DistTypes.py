from dataclasses import dataclass
#from dataclasses_json import dataclass_json
from dataclasses import asdict as dc_asdict
from typing import Optional, List

import numpy as np

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


def getDistTypeFromKwargs(**kwargs):
    """ This takes a set of kwargs and, using a present 'name' key, decides which DistType implementation they belong to
        And then returns an object of that DistType
    """
    if 'name' in kwargs.keys():    
        if kwargs['name'] == 'Wolf':
            return DT_Wolf.from_dict(kwargs)
        elif kwargs['name'] == 'ValueSpec':
            return DT_ValueSpec.from_dict(kwargs)
    else:
        return DistType.from_dict(kwargs)


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
class DT_Manual(DistType):  # Also pretty much empty, but can be used to allow user to add anything they want to save 
    def __post_init__(self):
        if self.name == "Empty": self.name = "Manual"

#@dataclass_json
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class DT_Wolf(DistType):
    """ Lambda channel spacing based on Wolf's notes 
            (ask Anthony Sciola or Frank Toffoletto for a copy)
        With the addition that there can be 2 p values for the start and end, and pStar transitions between them
    """
    p1: float = None
    p2: float = None

    def __post_init__(self):
        self.name = "Wolf"

    def genAlamsFromSpecies(self,sP):
        return self.genAlams(sP.n,sP.amin,sP.amax)

    def genAlams(self, n, amin, amax, kmin=0, kmax=-1):
        if kmax == -1: kmax = n

        alams = []
        for k in range(n):
            kfrac = (k-kmin)/(kmax-kmin)  # How far through the channel range are we
            pstar = (1-kfrac)*self.p1 + kfrac*self.p2
            lammax = amax-amin
            lam = lammax*((k - kmin + 0.5)/(kmax-kmin + 0.5))**pstar + amin
            alams.append(lam)
        return alams


#@dataclass_json
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class ValueSpec:
    start: float
    end: float
    scaleType: str  # See 'goodScaleTypes' below for valid strings
    n: Optional[int] = None  # Number of channels
    c: Optional[float] = None
    """ c has different meanings depending on scaleType
                lin: If given, will use set spacing c
                log: If given, will use log of base c
                spacing_lin: If given, will use as Sum_{k=1}^{N} c*k, where N is calculated based on start, end, and c
    """

    def __post_init__(self):
        goodScaleTypes = ['lin', 'log', 'spacing_lin']
        if self.scaleType not in goodScaleTypes:
            print("Error in ValueSpec, scaleType must be in {}, not {}. Defaulting to {}".format(goodScaleTypes, self.scaleType, goodScaleTypes[0]))
            self.scaleType = goodScaleTypes[0]
        if self.n is None and self.c is None:
            if self.scaleType == 'lin': self.c = 1
            print("Error in ValueSpec, must provide either (n) or (c). See source code to see what (c) does for each scaleType. Defaulting to " + str(self.c))
            

    def eval(self, doEnd):
        # Performs the appropriate operation given self's attributes and returns a list of values
        if self.scaleType == 'lin':
            line = np.linspace(sL.start, sL.end, sL.n, endpoint=doEnd)
        elif self.scaleType == 'log':
            lbase = self.c
            sign = 1 if self.start > 0 else -1
            start = np.log(np.abs(self.start))/np.log(lbase)
            end = np.log(np.abs(self.end))/np.log(lbase)
            line = np.logspace(start, end, self.n, base=lbase, endpoint=doEnd)
        elif self.scaleType == 'spacing_lin':
            diff = self.end-self.start
            if self.c is not None:
                #Set n based on c if needed
                #But also force n to be an integer
                self.n = int(0.5*(np.sqrt(8*diff/self.c + 1) + 1))
            #(Re)calculate c based on integer n
            self.c = 2*diff/(self.n**2 + self.n)
            print("Spacing_lin: n={}, c={}".format(self.n, self.c))

            spacings = np.array([self.c*k for k in range(self.n)])
            line = np.array([self.start + np.sum(spacings[:k]) for k in range(self.n)])

        return line

#@dataclass_json
@conditional_decorator(dataclass_json, dataclasses_json_module_imported)
@dataclass
class DT_ValueSpec(DistType):
    """ Lambda channel spacing based on a series of slope specifications
    """
    specList: List[ValueSpec] = None

    def __post_init__(self):
        self.name = "ValueSpec"
        #Check to see if all slopes are contiguous
        if len(self.specList) > 1:
            tol = 1E-4
            for i in range(len(self.specList)-1):
                if np.abs(self.specList[i].end - self.specList[i+1].start) > tol:
                    print("Error creating a DistType_SlopeSpec: SlopeSpec[{}].end ({}) != SlopeSpec[{}].start ({}). Undefined behavior"\
                        .format(i, self.specList[i].end, i+1, self.specList[i+1].start))

    def genAlamsFromSpecies(self, sP):
        #See if end points match up
        tol = 1E-4
        if np.abs(self.specList[0].start - sP.amin) > tol:
            print("SpecList[0].start={}, SpecParams.amin={}. Overwriting SpecParams.amin to SpecList[0].start".format(self.specList[0].start, sP.amin))
            sP.amin = self.specList[0].start
        if np.abs(self.specList[-1].end - sP.amax) > tol:
            print("SpecList[-1].end={}, SpecParams.amax={}. Overwriting SpecParams.amax to SpecList[-1].end".format(self.specList[0].start, sP.amin))
            sP.amax = self.specList[-1].end
        return self.genAlams(sP.n, sP.amin,sP.amax)

    def genAlams(self, n, amin, amax):
        nSL = len(self.specList)
        alams = np.array([])
        for i in range(nSL):
            doEnd = False if i < nSL-1 else True
            alams = np.append(alams, self.specList[i].eval(doEnd))
        return alams.tolist()







