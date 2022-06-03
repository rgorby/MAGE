# Custom
import kaipy.transform
from kaipy.solarWind.SolarWind import SolarWind

# 3rd party
import numpy
import cdasws as cdas

# Standard
import datetime
import re

#'WI_K0_SWE'
#'Np','V_GSM','QF_V','QF_Np','THERMAL_SPD'

#'WI_H0_MFI'
#'BGSM'

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


class WINDF(SolarWind):
    """
    OMNI Solar Wind file from CDAweb [http://cdaweb.gsfc.nasa.gov/].
    Data stored in GSE coordinates.
    """

    def __init__(self, fWIND):        
        SolarWind.__init__(self)

        self.bad_data = [-999.900, 
                         99999.9, # V
                         9999.99, # B
                         999.990, # density
                         1.00000E+07, # Temperature
                         9999999.,    # Temperature
                         99999, # Activity indices 
                         9999000, # SWE del_time
                         -1e31 # SWE & MFI                         
                         ]
        self.good_quality = [4098,14338]
        #self.bad_times =[   '17-03-2015 17:27:07.488',
        #                    '18-03-2015 21:53:49.140'
        #                ]
        #self.bad_fmt = '%d-%m-%Y %H:%M:%S.f'        
        #
        #self.bad_datetime = [   datetime.datetime(2015,3,17,hour=17,minute=27,second=7,microsecond=488*1000),
        #                        datetime.datetime(2015,3,18,hour=21,minute=53,second=49,microsecond=140*1000)
        #                    ]
        #
        #obtain 1 minute resolution observations from OMNI dataset
        print('Retrieving solar wind data from CDAWeb')
        self.__read(fWIND)

    def __read(self, fWIND):
        """
        Read the solar wind file & store results in self.data TimeSeries object.
        """
        (startDate, dates, data) = self.__readData(fWIND)
        (dataArray, hasBeenInterpolated) = self.__removeBadData(data)
        #(dataArray, hasBeenInterpolated) = self.__coarseFilter(dataArray, hasBeenInterpolated)
        self.__storeDataDict(dates, dataArray, hasBeenInterpolated)
        self.__appendMetaData(startDate, fWIND)
        self._appendDerivedQuantities()

        
    def __readData(self, fWIND):
        """
        return 2d array (of strings) containing data from file
        
        **TODO: read the fmt and figure out which column is which.  This would make things easier
        and more user friendly.  However, for now, we'll just assume the file is exactly 
        these quantities in this order
        """

        filedata = numpy.genfromtxt(fWIND)
        ntimes = numpy.shape(filedata)[0]
        
        yrs = filedata[:,0]
        doy = filedata[:,1]
        hrs = filedata[:,2]
        mns = filedata[:,3]
        bx  = filedata[:,4] #nT
        by  = filedata[:,5] #nT
        bz  = filedata[:,6] #nT
        vx  = filedata[:,7] #km/s
        vy  = filedata[:,8] #km/s
        vz  = filedata[:,9] #km/s
        n   = filedata[:,10] #n/cc
        T   = filedata[:,11] #K
        xBow   = filedata[:,12] #RE
        yBow   = filedata[:,13] #RE
        zBow   = filedata[:,14] #RE
        
        t0 = datetime.datetime(int(yrs[0]),1,1,hour=int(hrs[0]),minute=int(mns[0])) + datetime.timedelta(int(doy[0])-1)
        t1 = datetime.datetime(int(yrs[-1]),1,1,hour=int(hrs[-1]),minute=int(mns[-1])) + datetime.timedelta(int(doy[-1])-1)
        
        fOMNI = cdas.get_data(
               'sp_phys',
               'OMNI_HRO_1MIN',
               t0,
               t1,#+datetime.timedelta(minutes=tBuffer),
               ['BX_GSE,BY_GSE,BZ_GSE,Vx,Vy,Vz,proton_density,T,AE_INDEX,AL_INDEX,AU_INDEX,SYM_H,BSN_x,BSN_y,BSN_z']
            )
            
        tOMNI= fOMNI.get('EPOCH_TIME')         #datetime
        ae   = fOMNI.get('1-M_AE')             #nT
        al   = fOMNI.get('1-M_AL-INDEX')       #nT
        au   = fOMNI.get('AU-INDEX')           #nT
        symh = fOMNI.get('SYM/H_INDEX')        #nT
        ovx  = fOMNI.get('VX_VELOCITY,_GSE')   #kms
        ovy  = fOMNI.get('VY_VELOCITY,_GSE')   #kms
        ovz  = fOMNI.get('VZ_VELOCITY,_GSE')   #kms
        oxBow = fOMNI.get('X_(BSN),_GSE')       #km
        oyBow = fOMNI.get('Y_(BSN),_GSE')       #km
        ozBow = fOMNI.get('Z_(BSN),_GSE')       #km
        
        dates = []
        rows  = []
        
        for i in range(ntimes):
            currentTime = datetime.datetime(int(yrs[i]),1,1,hour=int(hrs[i]),minute=int(mns[i])) + datetime.timedelta(int(doy[i])-1)
            startTime = t0
            #calculating minutes from the start time
            nMin = self.__deltaMinutes(currentTime,startTime)

            data = [nMin,bx[i],by[i],bz[i],vx[i],vy[i],vz[i],n[i],T[i],ae[i],al[i],au[i],symh[i],xBow[i],yBow[i],zBow[i]]

            dates.append( currentTime )
            rows.append( data )

        return (startTime, dates, rows)

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

        # Temperature
        self.data.append('t', 'Temperature', r'$\mathrm{kK}$', dataArray[:,8]*1e-3)
        self.data.append('isTInterped', 'Is index i of T interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,7])

        ## Sound Speed (Thermal Speed)
        #self.data.append('cs', 'Sound speed', r'$\mathrm{km/s}$', dataArray[:,8])
        #self.data.append('isCsInterped', 'Is index i of Cs interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,7])

        # Activity Indices
        self.data.append('ae', 'AE-Index', r'$\mathrm{nT}$', dataArray[:,9])
        self.data.append('isAeInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,8])

        self.data.append('al', 'AL-Index', r'$\mathrm{nT}$', dataArray[:,10])
        self.data.append('isAlInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,9])

        self.data.append('au', 'AU-Index', r'$\mathrm{nT}$', dataArray[:,11])
        self.data.append('isAuInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,10])

        self.data.append('symh', 'SYM/H', r'$\mathrm{nT}$', dataArray[:,12])
        self.data.append('isSymHInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,11])

        # Bowshock Location
        self.data.append('xBS', 'BowShockX (gsm)', r'$\mathrm{RE}$', dataArray[:,13])
        self.data.append('isxBSInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,12])

        self.data.append('yBS', 'BowShockY (gsm)', r'$\mathrm{RE}$', dataArray[:,14])
        self.data.append('isyBSInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,13])

        self.data.append('zBS', 'BowShockZ (gsm)', r'$\mathrm{RE}$', dataArray[:,15])
        self.data.append('iszBSInterped', 'Is index i of N interpolated from bad data?', r'$\mathrm{boolean}$', hasBeenInterpolated[:,14])
        
        
    def __appendMetaData(self, date, filename):
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
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()
