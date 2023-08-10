"""
NOTE: This module is stripped down to only stuff that is currently used in kaipy.
The few remaining routines only wrap spacepy

This module provides coordinate transformations relevant to geospace modeling.
Many/most of these are wrappers to external C or Fortran libraries, employing
the NumPy vectorize() method to handle ndarray inputs when appropriate.

#K: Ripped out everything but,

x,y,z = SMtoGSM(x,y,z,dateTime)   - convert from solar magnetic to geocentric
                                    solar magnetospheric coordinates
x,y,z = GSMtoSM(x,y,z,dateTime)   - convert from geocentric solar magnetospheric
                                    to solar magnetic coordinates
x,y,z = GSEtoGSM(x,y,z,dateTime)  - convert from geocentric solar ecliptic to
                                    geocentric magnetospheric coordinates

"""

import numpy as np
import datetime
from spacepy.coordinates import Coords
from spacepy.time import Ticktock

def SMtoGSM(x,y,z,ut):
    """
    >>> SMtoGSM(1,2,3, datetime.datetime(2009,1,27,0,0,0)) # doctest:+ELLIPSIS
    (-0.126..., 2.0, 3.159...)
    """
    #Adapting code from scutils:convertGameraVec
    fromSys  = 'SM'
    fromType = 'car'
    toSys    = 'GSM'
    toType   = 'car'

    invec = Coords(np.column_stack((x,y,z)),fromSys,fromType, use_irbem=False)
    invec.ticks = Ticktock(ut)
    outvec = invec.convert(toSys,toType)

    return outvec.x[0],outvec.y[0],outvec.z[0]


def GSMtoSM(x,y,z,ut):
    """
    >>> GSMtoSM(1,2,3, datetime.datetime(2009,1,27,0,0,0)) # doctest:+ELLIPSIS
    (1.997..., 2.0, 2.451...)
    """
    #Adapting code from scutils:convertGameraVec
    fromSys  = 'GSM'
    fromType = 'car'
    toSys    = 'SM'
    toType   = 'car'

    invec = Coords(np.column_stack((x,y,z)),fromSys,fromType, use_irbem=False)
    invec.ticks = Ticktock(ut)
    outvec = invec.convert(toSys,toType)

    return outvec.x[0],outvec.y[0],outvec.z[0]


def GSEtoGSM(x,y,z,ut):
    """
    >>> GSEtoGSM(1,2,3, datetime.datetime(2009,1,27,0,0,0)) # doctest:+ELLIPSIS
    (0.99..., 0.540..., 3.564...)
    """
    #Adapting code from scutils:convertGameraVec
    fromSys  = 'GSE'
    fromType = 'car'
    toSys    = 'GSM'
    toType   = 'car'

    invec = Coords(np.column_stack((x,y,z)),fromSys,fromType, use_irbem=False)
    invec.ticks = Ticktock(ut)
    outvec = invec.convert(toSys,toType)

    return outvec.x[0],outvec.y[0],outvec.z[0]


