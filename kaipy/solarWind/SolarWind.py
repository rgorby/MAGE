import datetime

import numpy

import kaipy.transform
from kaipy.solarWind.TimeSeries import TimeSeries 
from kaipy.solarWind.ols import ols

class SolarWind(object):
    """
    Abstract base class for Solar Wind processing
    """

    def __init__(self):
        """
        Derived classes must read in Solar Wind data and store the
        results in a standard way via the a kaipy.TimeSeries object.
        Variables that must be set are:

        key       |  name                         |  units
        ----------------------------------------------
        n         |  Density                      |  #/cc
        vx        |  Velocity vector (gsm)        |  km/s
        vy        |  Velocity vector (gsm)        |  km/s
        vz        |  Velocity vector (gsm)        |  km/s
        t         |  Temperature                  |  k Kelvin
        cs        |  Sound speed                  |  km/s
        bx        |  Magnetic field vector (gsm)  |  nT
        by        |  Magnetic field vector (gsm)  |  nT
        bz        |  Magnetic field vector (gsm)  |  nT
        b         |  Magnetic field vector        |  nT
        time_min  |  Elapsed time since start     |  minutes
        -----------------------------------------------
        meta     | Metadata about the run.  Must  |  n/a
                 | include ['data']['Start Date'] |
                 | for coordinate transforms.     |

        Derived classes should make a call to
        SolarWind._appendDerivedQuantities() to compute some
        additional solar wind variables.
        """
        # The TimeSeries object stores all the Solar Wind data.
        self.data = TimeSeries()
        
            
    def _getTiltAngle(self, dateTime):
        """
        Get the tilt angle for the current Date & Time
        """
        (x,y,z) = kaipy.transform.SMtoGSM(0,0,1, dateTime)

        return numpy.arctan2(x,z)

    def _gsm2sm(self, dateTime, x,y,z):
        """
        Convert from GSM to SM coordinates for the current Date & Time
        """        
        return kaipy.transform.GSMtoSM(x,y,z, dateTime)


    def bxFit(self):
        """
        Compute & return coefficients for a multiple linear regression fit of Bx to By & Bz.

        Get the fit by applying the linear regression fit:
        """
        # Before doing anything, convert to SM coordinates.
        bx_sm = []
        by_sm = []
        bz_sm = []

        for i,time in enumerate(self.data.getData('time_min')):
            b_sm = self._gsm2sm(self.data.getData('meta')['Start date']+datetime.timedelta(minutes=time),
                                self.data.getData('bx')[i],
                                self.data.getData('by')[i],
                                self.data.getData('bz')[i])
            bx_sm.append(b_sm[0])
            by_sm.append(b_sm[1])
            bz_sm.append(b_sm[2])

        bx_sm = numpy.array(bx_sm)
        by_sm = numpy.array(by_sm)
        bz_sm = numpy.array(bz_sm)

        # Now that we're in SM, do the fit!
        y = bx_sm
        x = numpy.array([by_sm, bz_sm])

        linearFit = kaipy.solarWind.ols.ols(y,x.T, y_varnm='bx_sm', x_varnm=['by_sm','bz_sm'])
        ##Obtain information about the fit via the summary:
        #linearFit.summary()

        coef = linearFit.b

        ## Compute the variance... See p
        #xbar=numpy.average(y)
        #xi = coef[0]+coef[1]*by_sm+coef[2]*bz_sm
        #n=0
        #variance = 0.0
        #for i in range(len(xi)):
        #    n += 1
        #    variance += (xi[i]-xbar)**2
        #variance /= (n)
        #print 'variance is', variance
        #print 'std dev is', numpy.sqrt(variance)
        #
        ## Compute chi^2... See p.667-669 of Numerical Recipes in C++
        #chisq=0
        #xi = self.data.getData('time_min')
        #for i in range(len(y)):
        #    chisq += (y[i]-coef[0]-coef[1]*by_sm[i]-coef[2]*bz_sm[i])**2
        #print 'chisquared is', chisq
                                
        return coef

    def _appendDerivedQuantities(self):
        """
        Calculate & append standard derived quantities to the data dictionary.
        
        Note: single '_' underscore so this function can be called by derived classes
        """

        # --- Magnetic Field magnitude
        if 'b' not in self.data:
            b = numpy.sqrt(self.data.getData('bx')**2 +
                           self.data.getData('by')**2 +
                           self.data.getData('bz')**2)
            self.data.append('b', 'Magnitude of Magnetic Field', r'$\mathrm{nT}$', b)
        else:
            b = self.data.getData('b')

        # --- Velocity Field magnitude
        if 'v' not in self.data:
            v = numpy.sqrt(self.data.getData('vx')**2 +
                           self.data.getData('vy')**2 +
                           self.data.getData('vz')**2)
            self.data.append('v', 'Magnitude of Velocity', r'$\mathrm{km/s}$', v)
        else:
            v = self.data.getData('v')

        # -- Sound Speed
        if 'cs' not in self.data:
            try:
                cs = numpy.sqrt(5.0*1e3*self.data.getData('t')*(1.38e-23)/(3.0*1.67e-27)/(1.0e6))
                self.data.append('cs', 'Sound Speed', r'$\mathrm{km/s}$', cs)
            except KeyError:
                raise KeyError('Could not find temperature \'t\'.  Cannot compute sound speed (cs) without Temperature (t)!')

        # --- Alfven speed
        if 'va' not in self.data:
            va = (self.data.getData('b') * 1.0e-10 /
                  numpy.sqrt(1.97e-24*4*numpy.pi*
                             self.data.getData('n')) )
            self.data.append('va', 'Alfven Speed', r'$\mathrm{km/s}$', va)
        
        # --- Magnetosonic mach number (dimensionless)
        if 'ms' not in self.data:
            ms = v / self.data.getData('cs')
            self.data.append('ms', 'Magnetosonic Mach Number', '', ms)

        # --- Alfvenic Mach Number (dimensionless)
        if 'ma' not in self.data:
            ma = v/va
            self.data.append('ma', 'Alfvenic Mach Number', '', ma)

        # --- Temperature (Kelvin)
        if 't' not in self.data:
            t = 1e-3*(self.data.getData('cs')**2)*1.0e6*1.67e-27/1.38e-23
            self.data.append('t', 'Temperature', r'$\mathrm{kK}$', t)

        # --- Hours since start
        if 'time_hr' not in self.data:
            hr = self.data.getData('time_min')/60.0
            self.data.append('time_hr', 'Time (hours since start)', r'$\mathrm{hour}$', hr)

        # --- datetime
        if 'time' not in self.data:
            time = []
            for minute in self.data.getData('time_min'):
                time.append( self.data.getData('meta')['Start date'] + datetime.timedelta(minutes=minute) )
            self.data.append('time', 'Date and Time', r'$\mathrm{Date/Time}$', time)
            
        # --- Compute & store day of year
        if 'day' not in self.data:
            doy = []
            for dt in self.data.getData('time'):
                tt = dt.timetuple()
                dayFraction = (tt.tm_hour+tt.tm_min/60.+tt.tm_sec/(60.*60.))/24.
                doy.append( float(tt.tm_yday) + dayFraction )
            self.data.append('time_doy', 'Day of Year', r'$\mathrm{Day}$', doy)
