# Custom
import kaipy.transform
from kaipy.solarWind.SolarWind import SolarWind
from kaipy.solarWind.OMNI import OMNI
from kaipy.kdefs import *

# 3rd party
import numpy
import netCDF4 as nc
from cdasws import CdasWs

# Standard
import datetime
import re
import os

class DSCOVR(OMNI):
    """
    OMNI Solar Wind file from CDAweb [http://cdaweb.gsfc.nasa.gov/].
    Data stored in GSE coordinates.
    """

    def __init__(self,t0,filename = None, doFilter = False, sigmaVal = 3.0):        
        SolarWind.__init__(self)

        self.filter = doFilter
        self.sigma = sigmaVal

        self.bad_data = [-999.900, 
                         99999.9, # V
                         9999.99, # B
                         999.990, # density
                         1.00000E+07, # Temperature
                         9999999.0, # Temperature
                         99999 # Activity indices 
                         ]

        self.__read(filename,t0)

    def __read(self, filename,t0):
        """
        Read the solar wind file & store results in self.data TimeSeries object.
        """
        (startDate, dates, data, datanames) = self.__readData(filename,t0)
        (dataArray, hasBeenInterpolated) = self._removeBadData(data,datanames)
        if self.filter:
            (dataArray, hasBeenInterpolated) = self._coarseFilter(dataArray, hasBeenInterpolated)
        self._storeDataDict(dates, dataArray, hasBeenInterpolated)
        self.__appendMetaData(startDate, filename)
        self._appendDerivedQuantities()

    def __readData(self, filename,t0):
        """
        return 2d array (of strings) containing data from file
        
        **TODO: read the fmt and figure out which column is which.  This would make things easier
        and more user friendly.  However, for now, we'll just assume the file is exactly 
        these quantities in this order
        """

        filedata = numpy.genfromtxt(filename)
        ntimes = numpy.shape(filedata)[0]
        
        #yrs = filedata[:,0]
        #doy = filedata[:,1]
        #hrs = filedata[:,2]
        #mns = filedata[:,3]
        n   = filedata[:,0] #n/cc
        vx  = filedata[:,1] #km/s
        vy  = filedata[:,2] #km/s
        vz  = filedata[:,3] #km/s
        cs  = filedata[:,4] #km/s
        bx  = filedata[:,5] #nT
        by  = filedata[:,6] #nT
        bz  = filedata[:,7] #nT
        
        dates = []
        rows  = []

        timeshift = 46 #minutes
        startTime = t0 + datetime.timedelta(minutes=timeshift)
        dst = self._readDst()

        #dst = [0 ,-6 ,-10,-10,-8 ,-7 ,-8 ,-10,-9  ,-8 ,-4 ,0  ,2  ,2   ,0  ,-1 ,-3 ,-4 ,-4 ,2  ,28 ,10 ,-14,-9 , 
        #       -2,-20,-27,-32,-36,-43,-73,-81,-107,-95,-81,-85,-97,-103,-87,-81,-78,-72,-63,-56,-50,-42,-34,-36,
        #       -41,-43,-43,-42,-37,-13, -12, -23,  -35, -38, -38, -36, -32, -28, -25, -38,  -35, -38, -39, -36, -34, -26, -17,  -4]

        print("Starting Time: ",startTime.isoformat())

        for i in range(ntimes):
            #currentTime = datetime.datetime(int(yrs[i]),1,1,hour=int(hrs[i]),minute=int(mns[i])) + datetime.timedelta(int(doy[i])-1)
            currentTime = t0 + datetime.timedelta(minutes=i)
            #calculating minutes from the start time
            #nMin = self.__deltaMinutes(currentTime,startTime)
            nMin = i

            if currentTime < startTime:
              data = [nMin,bx[0],by[0],bz[0],vx[0],vy[0],vz[0],n[0],cs[0],0.,0.,0.,dst[int(i/60.)]]
            else:
              ic = int(i-timeshift)
              data = [nMin,bx[ic],by[ic],bz[ic],vx[ic],vy[ic],vz[ic],n[ic],cs[ic],0.,0.,0.,dst[int(i/60.)]]

            dates.append( currentTime )
            rows.append( data )

        datanames = ['Time','Bx','By','Bz','Vx','Vy','Vz','n','Sound Speed','AE','AL','AU','Dst']

        return (t0, dates, rows, datanames)

    def _storeDataDict(self, dates, dataArray, hasBeenInterpolated):
        """
        Populate self.data TimeSeries object via the 2d dataArray read from file.
        """
        #self.__gse2gsm(dates, dataArray)

        self.data.append('time_min', 'Time (Minutes since start)', 'min', dataArray[:,0])

        # Magnetic field
        self.data.append('bx', 'Bx (gsm)', r'$\mathrm{nT}$', dataArray[:,1])
        self.data.append('isBxInterped', 'Is index i of By interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,0])
        
        self.data.append('by', 'By (gsm)', r'$\mathrm{nT}$', dataArray[:,2])
        self.data.append('isByInterped', 'Is index i of By interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,1])

        self.data.append('bz', 'Bz (gsm)', r'$\mathrm{nT}$', dataArray[:,3])
        self.data.append('isBzInterped', 'Is index i of Bz interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,2])

        # Velocity
        self.data.append('vx', 'Vx (gsm)', r'$\mathrm{km/s}$', dataArray[:,4])
        self.data.append('isVxInterped', 'Is index i of Vx interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,3])

        self.data.append('vy', 'Vy (gsm)', r'$\mathrm{km/s}$', dataArray[:,5])
        self.data.append('isVyInterped', 'Is index i of Vy interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,4])

        self.data.append('vz', 'Vz (gsm)', r'$\mathrm{km/s}$', dataArray[:,6])
        self.data.append('isVzInterped', 'Is index i of Vz interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,5])

        # Density
        self.data.append('n', 'Density', r'$\mathrm{1/cm^3}$', dataArray[:,7])
        self.data.append('isNInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,6])

        # Sound Speed (Thermal Speed)
        self.data.append('cs', 'Sound speed', r'$\mathrm{km/s}$', dataArray[:,8])
        self.data.append('isCsInterped', 'Is index i of Cs interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,7])

        ## Temperature
        #self.data.append('t', 'Temperature', r'$\mathrm{kK}$', dataArray[:,8]*1e-3)
        #self.data.append('isTInterped', 'Is index i of T interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,7])

        ## Activity Indices
        #self.data.append('ae', 'AE-Index', r'$\mathrm{nT}$', dataArray[:,9])
        #self.data.append('isAeInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,8])

        #self.data.append('al', 'AL-Index', r'$\mathrm{nT}$', dataArray[:,10])
        #self.data.append('isAlInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,9])

        #self.data.append('au', 'AU-Index', r'$\mathrm{nT}$', dataArray[:,11])
        #self.data.append('isAuInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,10])

        
        # Bowshock Location
        self.data.append('xBS', 'BowShockX (gsm)', r'$\mathrm{RE}$', dataArray[:,9])
        self.data.append('isxBSInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,8])

        self.data.append('yBS', 'BowShockY (gsm)', r'$\mathrm{RE}$', dataArray[:,10])
        self.data.append('isyBSInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,9])

        self.data.append('zBS', 'BowShockZ (gsm)', r'$\mathrm{RE}$', dataArray[:,11])
        self.data.append('iszBSInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,10])

        # DST
        self.data.append('symh', 'SYM/H', r'$\mathrm{nT}$', dataArray[:,12])
        self.data.append('isSymHInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,11])

    def __appendMetaData(self, date, filename):
        """
        Add standard metadata to the data dictionary.
        """
        metadata = {'Model': 'CUSTOM',
                    'Source': filename,
                    'Date processed': datetime.datetime.now(),
                    'Start date': date
                    }
        
        self.data.append(key='meta',
                         name='Metadata for OMNI Solar Wind file',
                         units='n/a',
                         data=metadata)


class WIND(SolarWind):
    """
    OMNI Solar Wind file from CDAweb [http://cdaweb.gsfc.nasa.gov/].
    Data stored in GSE coordinates.
    """

    def __init__(self, fSWE,fMFI,fOMNI,xloc,tOffset,t0,t1):        
        SolarWind.__init__(self)

        self.bad_data = [-999.900, 
                         99999.9, # V
                         9999.99, # B
                         999.990, # density
                         1.00000E+07, # Temperature
                         99999, # Activity indices 
                         9999000, # SWE del_time
                         -1e31 # SWE & MFI                         
                         ]
        self.good_quality = [4098,14338]
        #self.bad_times =[   '17-03-2015 17:27:07.488',
        #                    '18-03-2015 21:53:49.140'
        #                ]
        #self.bad_fmt = '%d-%m-%Y %H:%M:%S.f'        

        self.bad_datetime = [   datetime.datetime(2015,3,17,hour=17,minute=27,second=7,microsecond=488*1000),
                                datetime.datetime(2015,3,18,hour=21,minute=53,second=49,microsecond=140*1000)
                            ]

        #obtain 1 minute resolution observations from OMNI dataset
        print('Retrieving solar wind data from CDAWeb')
        self.__read(fSWE,fMFI,fOMNI,xloc,tOffset,t0,t1)

    def __read(self, fSWE,fMFI,fOMNI,xloc,tOffset,t0,t1):
        """
        Read the solar wind file & store results in self.data TimeSeries object.
        """
        (SWEstartDate, MFIstartDate, OMNIstartDate, SWEdata, MFIdata, OMNIdata, SWEqf) = self.__readData(fSWE,fMFI,fOMNI,tOffset,t0,xloc)


        (SWEdata) = self.__checkGoodData(SWEdata,SWEqf)
        (SWEdataArray, SWEhasBeenInterpolated)   = self.__removeBadData(SWEdata )
        (MFIdataArray, MFIhasBeenInterpolated)   = self.__removeBadData(MFIdata )
        (OMNIdataArray, OMNIhasBeenInterpolated) = self.__removeBadData(OMNIdata)

        (SWEdataArray, SWEhasBeenInterpolated)   = self.__coarseFilter(SWEdataArray , SWEhasBeenInterpolated )
        (MFIdataArray, MFIhasBeenInterpolated)   = self.__coarseFilter(MFIdataArray , MFIhasBeenInterpolated )
        (OMNIdataArray, OMNIhasBeenInterpolated) = self.__coarseFilter(OMNIdataArray, OMNIhasBeenInterpolated)

        #(SWEdataArray, SWEhasBeenInterpolated)   = self.__windowedFilter(SWEdataArray , SWEhasBeenInterpolated )
        #(MFIdataArray, MFIhasBeenInterpolated)   = self.__windowedFilter(MFIdataArray , MFIhasBeenInterpolated )
        #(OMNIdataArray, OMNIhasBeenInterpolated) = self.__windowedFilter(OMNIdataArray, OMNIhasBeenInterpolated)

        (dates, dataArray, hasBeenInterpolated)  = self.__joinData(SWEdataArray, SWEhasBeenInterpolated, 
                                                                    MFIdataArray, MFIhasBeenInterpolated, 
                                                                    OMNIdataArray, OMNIhasBeenInterpolated,
                                                                    t0,t1)
        self.__storeDataDict(dates, dataArray, hasBeenInterpolated)
        self.__appendMetaData(t0, SWEstartDate, fSWE)
        self._appendDerivedQuantities()

        
    def __readData(self, fhSWE,fhMFI,fhOMNI,tOffset,t0,xloc):
        """
        return 2d array (of strings) containing data from file
        """
        #print('__readData')
        #pulling variables from file
        tSWE = fhSWE.get('EPOCH')               #datetime
        vx   = fhSWE.get('VX_(GSE)')            #km/s
        vy   = fhSWE.get('VY_(GSE)')            #km/s
        vz   = fhSWE.get('VZ_(GSE)')            #km/s
        qfv  = fhSWE.get('QF_V')                #
        qfn  = fhSWE.get('QF_NP')               #
        n    = fhSWE.get('ION_NP')              ##/cc
        cs    = fhSWE.get('SW_VTH')              #km/s

        tMFI = fhMFI.get('EPOCH')               #datetime
        bx   = fhMFI.get('BX_(GSE)')            #nT
        by   = fhMFI.get('BY_(GSE)')            #nT
        bz   = fhMFI.get('BZ_(GSE)')            #nT

        tOMNI= fhOMNI.get('EPOCH_TIME')         #datetime
        ovx  = fhOMNI.get('VX_VELOCITY,_GSE')   #kms
        ae   = fhOMNI.get('1-M_AE')             #nT
        al   = fhOMNI.get('1-M_AL-INDEX')       #nT
        au   = fhOMNI.get('AU-INDEX')           #nT
        symh = fhOMNI.get('SYM/H_INDEX')        #nT
        xBow = fhOMNI.get('X_(BSN),_GSE')       #km

        tshift = ((0 - xloc) / vx[0])/60. # t = (x - x_0)/Vx where X = 0, x_0 = xloc, and Vx is Vx in first data block in km/s.
        print('tshift:',tshift,xloc,vx[0])        

        SWEdates = []
        SWErows = []
        SWEqf = []
        SWEstartTime = tSWE[0]
        #badtimes = []
        #for i in range(len(self.bad_times)):
        #    badtimes.append(datetime.datetime.strptime(self.bad_times[i],self.bad_fmt))

        for i in range(len(tSWE)):
            for itime in range(len(self.bad_datetime)):          
                if abs(self.__deltaMinutes(tSWE[i],self.bad_datetime[itime])) <= 3./60.:
                    qfv[i] = 0
                    qfn[i] = 0
                
            #calculating minutes from the start time
            nMin = self.__deltaMinutes(tSWE[i],t0)+tOffset+tshift

            data = [nMin,n[i],vx[i],vy[i],vz[i],cs[i]]

            qf = [qfv[i],qfn[i]]

            SWEdates.append( tSWE[i] )
            SWErows.append ( data    )
            SWEqf.append   ( qf      )

        MFIdates = []
        MFIrows = []
        MFIstartTime = tMFI[0]
        for i in range(len(tMFI)):
          
            #calculating minutes from the start time
            nMin = self.__deltaMinutes(tMFI[i],t0)+tOffset+tshift

            data = [nMin,bx[i],by[i],bz[i]]

            MFIdates.append( tMFI[i] )
            MFIrows.append( data )

        OMNIdates = []
        OMNIrows = []
        for i in range(len(tOMNI)):
          
            OMNIstartTime = tOMNI[0]
            #calculating minutes from the start time
            nMin = self.__deltaMinutes(tOMNI[i],t0)

            data = [nMin,ae[i],al[i],au[i],symh[i]]

            OMNIdates.append( tOMNI[i] )
            OMNIrows.append( data )

        return ( SWEstartTime, MFIstartTime, OMNIstartTime, SWErows, MFIrows, OMNIrows, SWEqf )

    def __checkGoodData(self, data, qf):
        """
        Check the quality flag and set to bad data if bad data
        """
        nvar = len(data[0])
        nqf = len(qf[0])
        ntime = len(data)
        #print(numpy.shape(data),nvar,nqf,ntime)
        for itime in range(ntime):
            for iq in range(nqf):
                if qf[itime][iq] not in self.good_quality:
                    for ivar in range(1,nvar):
                        data[itime][ivar] = self.bad_data[-1]
        return ( data )

    def __removeBadData(self, data):
        """
        Linearly interpolate over bad data (defined by self.bad_data
        list) for each variable in dataStrs.
        
        data: 2d list.  Each row is a list containing:
          [nMinutes, Bx, By, Bz, Vx, Vy, Vz, rho, temp, ae, al, au, symh]

        Returns:
          data: interpolated floating-point numpy array
          hasBeenInterpolated: 2d array that identifies if bad values were removed/interpolated.

        NOTE: This is remarkably similar to __coarseFilter!
          Refactoring to keep it DRY wouldn't be a bad idea. . .
        """
        #assert( len(data[0]) == 13 )
        nvar = len(data[0])
        hasBeenInterpolated = numpy.empty((len(data), nvar-1))
        hasBeenInterpolated.fill(False)

        for varIdx in range(1,nvar):

            lastValidIndex = -1
            for curIndex,row in enumerate(data):
                if row[varIdx] in self.bad_data:
                    # This item has bad data.
                    hasBeenInterpolated[curIndex, varIdx-1] = True
                    if (lastValidIndex == -1) & (curIndex == len(data)-1):
                        # Data does not have at least one valid element!
                        # Setting all values to 0 so that file can still be made
                        print("No good elements, setting all values to 0 for variable ID: ", varIdx)
                        data[curIndex][varIdx] = 0.
                        #raise Exception("First & Last datapoint(s) in OMNI "+
                        #                  "solar wind file are invalid.  Not sure "+
                        #                  "how to interpolate across bad data.")
                    elif (curIndex == len(data)-1):
                        # Clamp last bad data to previous known good data.
                        data[curIndex][varIdx] = data[lastValidIndex][varIdx]
                    else:
                        # Note the bad data & skip this element for now.
                        # We will linearly interpolate between valid data
                        continue

                # At this point, curIndex has good data.
                if (lastValidIndex+1) == curIndex:
                    # Set current element containing good data.
                    data[curIndex][varIdx] = float( row[varIdx] )
                else:
                    # If first index is invalid, clamp to first good value.
                    if lastValidIndex == -1:
                        lastValidIndex = 0
                        data[lastValidIndex][varIdx] = data[curIndex][varIdx]

                    # Linearly interpolate over bad data.
                    interpolated = numpy.interp(range(lastValidIndex, curIndex), # x-coords of interpolated values
                                                [lastValidIndex, curIndex],  # x-coords of data.
                                                [float(data[lastValidIndex][varIdx]), float(data[curIndex][varIdx])]) # y-coords of data.
                    # Store the results.
                    for j,val in enumerate(interpolated):
                        data[lastValidIndex+j][varIdx] = val
                lastValidIndex = curIndex

        return (numpy.array(data, numpy.float), hasBeenInterpolated)

    def __coarseFilter(self, dataArray, hasBeenInterpolated):
        """
         Use coarse noise filtering to remove values outside 3
         deviations from mean of all values in the plotted time
         interval.

         Parameters:

           dataArray: 2d numpy array.  Each row is a list
             containing [nMinutes, Bx, By, Bz, Vx, Vy, Vz, rho, temp, ae, al, au, symh]

           hasBeenInterpolated: 2d boolean list.  Each row is a list
             of boolean values denoting whether dataArray[:,1:9] was
             derived/interpolated from the raw data (ie. bad points
             removed).

         Output:
           dataArray:  same structure as input array with bad elements removed
           hasBeenInterpolated: same as input array with interpolated values stored.

        NOTE: This is remarkably similar to __removeBadData!
          Refactoring to keep it DRY wouldn't be a bad idea. . .
        """
        
        nvar = len(dataArray[0])

        stds = []
        means = []
        for varIdx in range(1,nvar):
            stds.append( dataArray[:,varIdx].std() )
            means.append( dataArray[:,varIdx].mean() )
            
            # Linearly interpolate over data that exceeds 3 standard
            # deviations from the mean
            lastValidIndex = -1
            for curIndex,row in enumerate(dataArray):
                # Are we outside 3 sigma from mean?
                if abs(means[varIdx-1] - row[varIdx]) > 3*stds[varIdx-1]:
                    hasBeenInterpolated[curIndex, varIdx-1] = True
                    if (curIndex == len(dataArray)-1):
                        # Clamp last bad data to previous known good data.
                        dataArray[curIndex][varIdx] = dataArray[lastValidIndex][varIdx]
                    else:
                        # Note the bad data & skip this element for now.
                        # We will linearly interpolate between valid data
                        continue

                if (lastValidIndex+1) != curIndex:
                    # If first index is invalid, clamp to first good value.
                    if lastValidIndex == -1:
                        lastValidIndex = 0
                        dataArray[lastValidIndex][varIdx] = dataArray[curIndex][varIdx]

                    # Linearly interpolate over bad data.
                    interpolated = numpy.interp(range(lastValidIndex, curIndex), # x-coords of interpolated values
                                                [lastValidIndex, curIndex],  # x-coords of data.
                                                [float(dataArray[lastValidIndex][varIdx]), float(dataArray[curIndex][varIdx])]) # y-coords of data.
                    # Store the results.
                    for j,val in enumerate(interpolated):
                        dataArray[lastValidIndex+j][varIdx] = val
                lastValidIndex = curIndex

        return (dataArray, hasBeenInterpolated)


    def __windowedFilter(self, dataArray, hasBeenInterpolated):
        """
         Use coarse noise filtering to remove values outside 3
         deviations from mean of all values in the plotted time
         interval.

         Parameters:

           dataArray: 2d numpy array.  Each row is a list
             containing [nMinutes, Bx, By, Bz, Vx, Vy, Vz, rho, temp, ae, al, au, symh]

           hasBeenInterpolated: 2d boolean list.  Each row is a list
             of boolean values denoting whether dataArray[:,1:9] was
             derived/interpolated from the raw data (ie. bad points
             removed).

         Output:
           dataArray:  same structure as input array with bad elements removed
           hasBeenInterpolated: same as input array with interpolated values stored.

        NOTE: This is remarkably similar to __removeBadData!
          Refactoring to keep it DRY wouldn't be a bad idea. . .
        """
        
        nvar = len(dataArray[0])
        windowsize = 60*4

        stds = []
        means = []
        for varIdx in range(1,nvar):
            stds.append( dataArray[:,varIdx].std() )
            means.append( dataArray[:,varIdx].mean() )
            
            # Linearly interpolate over data that exceeds 3 standard
            # deviations from the mean
            lastValidIndex = -1
            for curIndex,row in enumerate(dataArray):
                # Are we outside 3 sigma from mean?
                if abs(means[varIdx-1] - row[varIdx]) > 3*stds[varIdx-1]:
                    hasBeenInterpolated[curIndex, varIdx-1] = True
                    if (curIndex == len(dataArray)-1):
                        # Clamp last bad data to previous known good data.
                        dataArray[curIndex][varIdx] = dataArray[lastValidIndex][varIdx]
                    else:
                        # Note the bad data & skip this element for now.
                        # We will linearly interpolate between valid data
                        continue

                if (lastValidIndex+1) != curIndex:
                    # If first index is invalid, clamp to first good value.
                    if lastValidIndex == -1:
                        lastValidIndex = 0
                        dataArray[lastValidIndex][varIdx] = dataArray[curIndex][varIdx]

                    # Linearly interpolate over bad data.
                    interpolated = numpy.interp(range(lastValidIndex, curIndex), # x-coords of interpolated values
                                                [lastValidIndex, curIndex],  # x-coords of data.
                                                [float(dataArray[lastValidIndex][varIdx]), float(dataArray[curIndex][varIdx])]) # y-coords of data.
                    # Store the results.
                    for j,val in enumerate(interpolated):
                        dataArray[lastValidIndex+j][varIdx] = val
                lastValidIndex = curIndex

        return (dataArray, hasBeenInterpolated)

    def __joinData(self, SWEdataArray, SWEhasBeenInterpolated, MFIdataArray, MFIhasBeenInterpolated, OMNIdataArray, OMNIhasBeenInterpolated,t0,t1):
        #print('joinData')
        ntime = self.__deltaMinutes(t1,t0)
        nMin = range(int(ntime))
        n  = numpy.interp(nMin,SWEdataArray[:,0],SWEdataArray[:,1])
        vx = numpy.interp(nMin,SWEdataArray[:,0],SWEdataArray[:,2])
        vy = numpy.interp(nMin,SWEdataArray[:,0],SWEdataArray[:,3])
        vz = numpy.interp(nMin,SWEdataArray[:,0],SWEdataArray[:,4])
        cs = numpy.interp(nMin,SWEdataArray[:,0],SWEdataArray[:,5])
        bx = numpy.interp(nMin,MFIdataArray[:,0],MFIdataArray[:,1])
        by = numpy.interp(nMin,MFIdataArray[:,0],MFIdataArray[:,2])
        bz = numpy.interp(nMin,MFIdataArray[:,0],MFIdataArray[:,3])
        ae = numpy.interp(nMin,OMNIdataArray[:,0],OMNIdataArray[:,1])
        al = numpy.interp(nMin,OMNIdataArray[:,0],OMNIdataArray[:,2])
        au = numpy.interp(nMin,OMNIdataArray[:,0],OMNIdataArray[:,3])
        symh = numpy.interp(nMin,OMNIdataArray[:,0],OMNIdataArray[:,4])
        nI  = numpy.interp(nMin,SWEdataArray[:,0],SWEhasBeenInterpolated[:,0])
        vxI = numpy.interp(nMin,SWEdataArray[:,0],SWEhasBeenInterpolated[:,1])
        vyI = numpy.interp(nMin,SWEdataArray[:,0],SWEhasBeenInterpolated[:,2])
        vzI = numpy.interp(nMin,SWEdataArray[:,0],SWEhasBeenInterpolated[:,3])
        csI = numpy.interp(nMin,SWEdataArray[:,0],SWEhasBeenInterpolated[:,4])
        bxI = numpy.interp(nMin,MFIdataArray[:,0],MFIhasBeenInterpolated[:,0])
        byI = numpy.interp(nMin,MFIdataArray[:,0],MFIhasBeenInterpolated[:,1])
        bzI = numpy.interp(nMin,MFIdataArray[:,0],MFIhasBeenInterpolated[:,2])
        aeI = numpy.interp(nMin,OMNIdataArray[:,0],OMNIhasBeenInterpolated[:,0])
        alI = numpy.interp(nMin,OMNIdataArray[:,0],OMNIhasBeenInterpolated[:,1])
        auI = numpy.interp(nMin,OMNIdataArray[:,0],OMNIhasBeenInterpolated[:,2])
        symhI = numpy.interp(nMin,OMNIdataArray[:,0],OMNIhasBeenInterpolated[:,3])
        
        dates = []
        dataArray = []
        interped = []
        hasBeenInterpolated = []
        for i in nMin:
            dates.append(t0+datetime.timedelta(minutes=i))

            arr = [nMin[i],bx[i],by[i],bz[i],vx[i],vy[i],vz[i],n[i],cs[i],ae[i],al[i],au[i],symh[i]]
            dataArray.append(arr)
            arr = [bxI[i],byI[i],bzI[i],vxI[i],vyI[i],vzI[i],nI[i],csI[i],aeI[i],alI[i],auI[i],symhI[i]]
            hasBeenInterpolated.append(arr)
    
        return (dates, numpy.array(dataArray,numpy.float), numpy.array(hasBeenInterpolated))

    def __storeDataDict(self, dates, dataArray, hasBeenInterpolated):
        """
        Populate self.data TimeSeries object via the 2d dataArray read from file.
        """
        #print('__storeDataDict')
        self.__gse2gsm(dates, dataArray)

        #print(numpy.shape(dataArray))
        #print(numpy.shape(dataArray[:,0]))
        self.data.append('time_min', 'Time (Minutes since start)', 'min', dataArray[:,0])

        # Magnetic field
        self.data.append('bx', 'Bx (gsm)', r'$\mathrm{nT}$', dataArray[:,1])
        self.data.append('isBxInterped', 'Is index i of By interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,0])
        
        self.data.append('by', 'By (gsm)', r'$\mathrm{nT}$', dataArray[:,2])
        self.data.append('isByInterped', 'Is index i of By interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,1])

        self.data.append('bz', 'Bz (gsm)', r'$\mathrm{nT}$', dataArray[:,3])
        self.data.append('isBzInterped', 'Is index i of Bz interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,2])

        # Velocity
        self.data.append('vx', 'Vx (gsm)', r'$\mathrm{km/s}$', dataArray[:,4])
        self.data.append('isVxInterped', 'Is index i of Vx interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,3])

        self.data.append('vy', 'Vy (gsm)', r'$\mathrm{km/s}$', dataArray[:,5])
        self.data.append('isVyInterped', 'Is index i of Vy interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,4])

        self.data.append('vz', 'Vz (gsm)', r'$\mathrm{km/s}$', dataArray[:,6])
        self.data.append('isVzInterped', 'Is index i of Vz interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,5])

        # Density
        self.data.append('n', 'Density', r'$\mathrm{1/cm^3}$', dataArray[:,7])
        self.data.append('isNInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,6])

        ## Temperature
        #self.data.append('t', 'Temperature', r'$\mathrm{kK}$', dataArray[:,8]*1e-3)
        #self.data.append('isTInterped', 'Is index i of T interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,7])

        # Sound Speed (Thermal Speed)
        self.data.append('cs', 'Sound speed', r'$\mathrm{km/s}$', dataArray[:,8])
        self.data.append('isCsInterped', 'Is index i of Cs interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,7])

        # Activity Indices
        self.data.append('ae', 'AE-Index', r'$\mathrm{nT}$', dataArray[:,9])
        self.data.append('isAeInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,8])

        self.data.append('al', 'AL-Index', r'$\mathrm{nT}$', dataArray[:,10])
        self.data.append('isAlInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,9])

        self.data.append('au', 'AU-Index', r'$\mathrm{nT}$', dataArray[:,11])
        self.data.append('isAuInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,10])

        self.data.append('symh', 'SYM/H', r'$\mathrm{nT}$', dataArray[:,12])
        self.data.append('isSymHInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,11])
        
    def __appendMetaData(self, date, dateshift, filename):
        """
        Add standard metadata to the data dictionary.
        """
        metadata = {'Model': 'WIND',
                    'Source': filename,
                    'Date processed': datetime.datetime.now(),
                    'Start date': date,
                    }
        
        self.data.append(key='meta',
                         name='Metadata for WIND Solar Wind file',
                         units='n/a',
                         data=metadata)

    
    def __deltaMinutes(self, t1, startDate):
        """
        Returns: Number of minutes elapsed between t1 and startDate.
        """
        diff = t1 - startDate

        return (diff.days*24.0*60.0 + diff.seconds/60.0)

    def __gse2gsm(self, dates, dataArray):
        """
        Transform magnetic field B and velocity V from GSE to GSM
        coordinates.  Store results by overwriting dataArray contents.
        """
        for i,data in enumerate(dataArray):
            d = dates[i]

            # Update magnetic field
            b_gsm = kaipy.transform.GSEtoGSM(data[1], data[2], data[3], d)        
            data[1] = b_gsm[0]
            data[2] = b_gsm[1]
            data[3] = b_gsm[2]

            # Update Velocity
            v_gsm = kaipy.transform.GSEtoGSM(data[4], data[5], data[6], d)
            data[4] = v_gsm[0]
            data[5] = v_gsm[1]
            data[6] = v_gsm[2]

class DSCOVRNC(OMNI):
    """
    OMNI Solar Wind file from CDAweb [http://cdaweb.gsfc.nasa.gov/].
    Data stored in GSE coordinates.
    """

    def __init__(self,t0,t1,doFilter = False, sigmaVal = 3.0):        
        SolarWind.__init__(self)

        self.filter = doFilter
        self.sigma = sigmaVal

        self.bad_data = [-999.900, 
                         99999.9, # V
                         9999.99, # B
                         999.990, # density
                         1.00000E+07, # Temperature
                         9999999.0, # Temperature
                         99999, # Activity indices 
                         -99999,
                         1e+20
                         ]

        self.__read(t0,t1)

    def __read(self, t0,t1):
        """
        Read the solar wind file & store results in self.data TimeSeries object.
        """
        (startDate, dates, data) = self.__readData(t0,t1)
        (dataArray, hasBeenInterpolated) = self._removeBadData(data)
        if self.filter:
            (dataArray, hasBeenInterpolated) = self._coarseFilter(dataArray, hasBeenInterpolated)
        self._storeDataDict(dates, dataArray, hasBeenInterpolated)
        self.__appendMetaData(startDate)
        self._appendDerivedQuantities()

    def __readData(self, t0,t1):
        """
        return 2d array (of strings) containing data from file
        
        **TODO: read the fmt and figure out which column is which.  This would make things easier
        and more user friendly.  However, for now, we'll just assume the file is exactly 
        these quantities in this order
        """

        filelist = os.listdir()
        pop = []
        f1m = []
        m1m = []
        fmt = '%Y%m%d'
        jud0 = datetime.datetime(1970,1,1,0,0,0,0)
        
        for f in filelist:
            if f[0:2] == 'oe':
                ctime = datetime.datetime.strptime(f[15:23],fmt)
                if ctime >= t0 and ctime <=t1:
                    if 'pop' in f:
                        pop.append(f)
                    if 'f1m' in f:
                        f1m.append(f)
                    if 'm1m' in f:
                        m1m.append(f)
                        
        pop = np.sort(pop)
        f1m = np.sort(f1m)
        m1m = np.sort(m1m)
        
        if len(pop) != len(f1m) or len(f1m) != len(m1m) or len(pop) != len(m1m):
            raise Exception('file list not the same')
        if len(pop) == 0 or len(f1m) == 0 or len(m1m) == 0:
            raise Exception('missing files for this daterange')
        
        mtime = []
        ftime = []
        ptime = []
        n = []
        vx = []
        vy = []
        vz = []
        temp = []
        bx = []
        by = []
        bz = []
        satx = []
        for i in range(len(pop)):
            pfn = pop[i]
            ffn = f1m[i]
            mfn = m1m[i]
            pds = nc.Dataset(pfn) #time, sat_x_gse
            fds = nc.Dataset(ffn) #time,proton_density, proton_vx_gse, proton_vy_gse, proton_vz_gse, proton_temperature
            mds = nc.Dataset(mfn) #time, bx_gse, by_gse, bz_gse
            for k in range(len(mds['time'])):
                mtime.append(jud0 + datetime.timedelta(milliseconds=mds['time'][:][k]))
                bx.append(mds['bx_gse'][:][k])
                by.append(mds['by_gse'][:][k])
                bz.append(mds['bz_gse'][:][k])
            for k in range(len(fds['time'])):
                ftime.append(jud0 + datetime.timedelta(milliseconds=fds['time'][:][k]))
                '''
                if fds['overall_quality'][:][k] == 0:
                    n.append(fds['proton_density'][:][k])
                    vx.append(fds['proton_vx_gse'][:][k])
                    vy.append(fds['proton_vy_gse'][:][k])
                    vz.append(fds['proton_vz_gse'][:][k])
                    temp.append(fds['proton_temperature'][:][k])
                else:
                    n.append(numpy.nan)
                    vx.append(numpy.nan)
                    vy.append(numpy.nan)
                    vz.append(numpy.nan)
                    temp.append(numpy.nan)
                '''
                n.append(fds['proton_density'][:][k])
                vx.append(fds['proton_vx_gse'][:][k])
                vy.append(fds['proton_vy_gse'][:][k])
                vz.append(fds['proton_vz_gse'][:][k])
                temp.append(fds['proton_temperature'][:][k])
            for k in range(len(pds['time'])):
                ptime.append(jud0 + datetime.timedelta(milliseconds=pds['time'][:][k]))
                satx.append(pds['sat_x_gse'][:][k])
        
        dates = []
        rows  = []

        #timeshift = 46 #minutes
        # simple projectile motion
        # t = (x - x0)/v
        timeshift = np.int(np.round((np.mean(satx)*-1)/(np.nanmean(vx))/60.0))
        startTime = t0 + datetime.timedelta(minutes=timeshift)
        #endTime = t0 + datetime.timedelta(minutes=timeshift+ntimes)
        endTime = t1
        dsttime,dst = self._readDst(t0,t1)
        ntimes = t1 - t0
        ntimes = int(ntimes.total_seconds()/60.0)

        print("Starting Time: ",startTime.isoformat())
        print("We are using a constant timeshift of: ", timeshift ," minutes")
        #itp = 0 #ptime
        itf = 0 #ftime
        itm = 0 #mtime
        itd = 0 #dsttime

        for i in range(ntimes):
            #currentTime = datetime.datetime(int(yrs[i]),1,1,hour=int(hrs[i]),minute=int(mns[i])) + datetime.timedelta(int(doy[i])-1)
            currentTime = t0 + datetime.timedelta(minutes=i)
            #calculating minutes from the start time
            #nMin = self.__deltaMinutes(currentTime,startTime)
            while(mtime[itm] + datetime.timedelta(minutes=timeshift) < currentTime):
                itm = itm+1
            while(ftime[itf] + datetime.timedelta(minutes=timeshift) < currentTime):
                itf = itf+1
            while(dsttime[itd] < currentTime):
                itd = itd+1
            nMin = i
            
            data = [nMin,bx[itm],by[itm],bz[itm],vx[itf],vy[itf],vz[itf],n[itf],temp[itf],0,0,0,dst[itd],0,0,0]

            #if currentTime < startTime:
            #  data = [nMin,bx[0],by[0],bz[0],vx[0],vy[0],vz[0],n[0],temp[0],0.,0.,0.,dst[int(i/60.)]]
            #else:
            #  ic = int(i-timeshift)
            #  data = [nMin,bx[ic],by[ic],bz[ic],vx[ic],vy[ic],vz[ic],n[ic],temp[ic],0.,0.,0.,dst[int(i/60.)]]

            dates.append( currentTime )
            rows.append( data )

        return (t0, dates, rows)

    def _readDst(self,startTime,endTime):
        dstfile = open("dst.dat",'r')
        text = dstfile.readlines()
        for i,j in enumerate(text):
            if j[0] == '2':
                iskip = i
                break
        dstfile.close()

        dat = np.genfromtxt("dst.dat",skip_header=iskip, autostrip=True,dtype=None)
        dsttime = []
        dst = []
        fmt='%Y-%m-%dT%H:%M:%S.000'
        for i in dat:
            timestr = i[0].decode()+"T"+i[1].decode()
            currenttime = datetime.datetime.strptime(timestr,fmt)
            if currenttime >= startTime and currenttime <= endTime:
                dsttime.append(currenttime)
                dst.append(i[3])
                
        return (dsttime, dst)
            
    def __appendMetaData(self, date):
        """
        Add standard metadata to the data dictionary.
        """
        metadata = {'Model': 'CUSTOM',
                    'Source': 'NOAA DSCOVR NC',
                    'Date processed': datetime.datetime.now(),
                    'Start date': date
                    }
        
        self.data.append(key='meta',
                         name='Metadata for OMNI Solar Wind file',
                         units='n/a',
                         data=metadata)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
