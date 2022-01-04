from dataclasses import dataclass
from dataclasses_json import dataclass_json
from dataclasses import asdict as dc_asdict
from typing import Optional, List

import numpy as np


def getDistTypeFromKwargs(**kwargs):
    """ This takes a set of kwargs and, using a present 'name' key, decides which DistType implementation they belong to
        And then returns an object of that DistType
    """
    if 'name' in kwargs.keys():    
        if kwargs['name'] == 'Wolf':
            return DT_Wolf(**kwargs)
    else:
        return DistType(**kwargs)


#------
# Parameters needed to determine lambda distribution
#------

@dataclass
class DistType:  # Empty class just so we can force the type in dataclasses below
    name: str = "Empty"
    
#------
# Specific implementations of DistType
#------
@dataclass_json
@dataclass
class DT_Manual(DistType):  # Also pretty much empty, but can be used to allow user to add anything they want to save 
    def __post_init__(self):
        if self.name == "Empty": self.name = "Manual"

@dataclass_json
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


@dataclass_json
@dataclass
class SlopeSpec:
    n: int
    start: float
    end: float
    slopeType: str  # Accepts 'lin' and 'log'
    base: Optional[float] = 10

    def __post_init__(self):
        goodSlopeTypes = ['lin', 'log']
        if self.slopeType not in goodSlopeTypes:
            print("Error in SlopeSpec, slopeType must be in {}, not {}. Defaulting to {}".format(goodSlopeTypes, self.slopeType, goodSlopeTypes[0]))
            self.slopeType = goodSlopeTypes[0]

@dataclass_json
@dataclass
class DT_SlopeSpec:
    """ Lambda channel spacing based on a series of slope specifications
    """
    specList: List[SlopeSpec]

    def __post_init(self):
        self.name = "SlopeSpec"
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
            sL = self.specList[i]
            if sL.slopeType == 'lin':
                line = np.linspace(sL.start, sL.end, sL.n, endpoint=doEnd)
            if sL.slopeType == 'log':
                lbase = sL.base
                sign = 1 if sL.start > 0 else -1
                start = np.log(np.abs(sL.start))/np.log(lbase)
                end = np.log(np.abs(sL.end))/np.log(lbase)
                line = np.logspace(start, end, sL.n, base=lbase, endpoint=doEnd)
            alams = np.append(alams, line)
        return alams














