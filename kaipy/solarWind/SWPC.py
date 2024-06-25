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

class DSCOVRNC(OMNI):
    """
    OMNI Solar Wind file from CDAweb [http://cdaweb.gsfc.nasa.gov/].
    Data stored in GSE coordinates.
    """

    def __init__(self,t0,t1,doFilter = False, sigmaVal = 3.0,dodcx=False,dcxfile=None):
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
        self.__read(t0,t1,dodcx=dodcx,dcxfile=dcxfile)

    def __read(self, t0,t1,dodcx=False,dcxfile=None):
        """
        Read the solar wind file & store results in self.data TimeSeries object.
        """
        (startDate, dates, data) = self.__readData(t0,t1,dodcx=dodcx,dcxfile=dcxfile)
        (dataArray, hasBeenInterpolated) = self._removeBadData(data)
        if self.filter:
            (dataArray, hasBeenInterpolated) = self._coarseFilter(dataArray, hasBeenInterpolated)
        self._storeDataDict(dates, dataArray, hasBeenInterpolated)
        self.__appendMetaData(startDate)
        self._appendDerivedQuantities()

    def __readData(self, t0,t1,dodcx=False,dcxfile=None):
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
        fmt1 = '%Y%m%d%H%M%S'
        fmt2 = '%Y%m%d%H%M%S'
        jud0 = datetime.datetime(1970,1,1,0,0,0,0)
        
        for f in filelist:
            if f[0:2] == 'oe':
                ctime = datetime.datetime.strptime(f[15:29],fmt1)
                etime = datetime.datetime.strptime(f[31:45],fmt2)
                if (ctime >= t0 and ctime <=t1) or (t0 <= ctime and ctime <= t1) or (t0 <= etime and etime <= t1):
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
        timeshift = int(np.round((np.mean(satx)*-1)/(np.nanmean(vx))/60.0))
        startTime = t0 + datetime.timedelta(minutes=timeshift)
        #endTime = t0 + datetime.timedelta(minutes=timeshift+ntimes)
        endTime = t1
        if dodcx or dcxfile != None:
            dsttime,dst = self._readDcx(t0,t1,dcxfile)
        else:
            #dsttime,dst = self._readDst(t0,t1)
            dsttime,dst = self._getDst(t0,t1)
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
            
    def _readDcx(self,startTime,endTime,dcxstr = None):
        if dcxstr == None:
            dcxstr = "dcx.txt"
        dat = np.genfromtxt(dcxstr,autostrip=True,dtype=None)
        dsttime = []
        dst = []
        for i in range(len(dat)):
            currenttime = datetime.datetime(dat[i][1],dat[i][2],dat[i][3],dat[i][4])
            print(i,currenttime,startTime,endTime)
            if currenttime >= startTime and currenttime <= endTime:
                dsttime.append(currenttime)
                dst.append(dat[i][5])
                print(currenttime,dat[i][5])
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


class ACESWPC(DSCOVRNC):
    """
    OMNI Solar Wind file from CDAweb [http://cdaweb.gsfc.nasa.gov/].
    Data stored in GSE coordinates.
    """

    def __init__(self,t0,t1,doFilter = False, sigmaVal = 3.0,dodcx=False,dcxfile=None):
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
        self.__read(t0,t1,dodcx=dodcx,dcxfile=dcxfile)

    def __read(self, t0,t1,dodcx=False,dcxfile=None):
        """
        Read the solar wind file & store results in self.data TimeSeries object.
        """
        (startDate, dates, data) = self.__readData(t0,t1,dodcx=dodcx,dcxfile=dcxfile)
        (dataArray, hasBeenInterpolated) = self._removeBadData(data)
        if self.filter:
            (dataArray, hasBeenInterpolated) = self._coarseFilter(dataArray, hasBeenInterpolated)
        self._storeDataDict(dates, dataArray, hasBeenInterpolated)
        self.__appendMetaData(startDate)
        self._appendDerivedQuantities()

    def __readData(self, t0,t1,dodcx=False,dcxfile=None):
        """
        return 2d array (of strings) containing data from file

        **TODO: read the fmt and figure out which column is which.  This would make things easier
        and more user friendly.  However, for now, we'll just assume the file is exactly
        these quantities in this order
        """

        self.__downloadACE(t0,t1)

        filelist = os.listdir()
        pop = []
        f1m = []
        m1m = []
        fmt1 = '%Y%m%d%H%M%S'
        fmt2 = '%Y%m%d%H%M%S'
        jud0 = datetime.datetime(1970,1,1,0,0,0,0)

        for f in filelist:
            if f[0:2] == 'oe':
                ctime = datetime.datetime.strptime(f[15:29],fmt1)
                etime = datetime.datetime.strptime(f[31:45],fmt2)
                if (ctime >= t0 and ctime <=t1) or (t0 <= ctime and ctime <= t1) or (t0 <= etime and etime <= t1):
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
        timeshift = int(np.round((np.mean(satx)*-1)/(np.nanmean(vx))/60.0))
        startTime = t0 + datetime.timedelta(minutes=timeshift)
        #endTime = t0 + datetime.timedelta(minutes=timeshift+ntimes)
        endTime = t1
        if dodcx or dcxfile != None:
            dsttime,dst = self._readDcx(t0,t1,dcxfile)
        else:
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

    def __downloadACE(self, t0, t1):

        swpcdir = "https://sohoftp.nascom.nasa.gov/sdb/goes/ace/daily/"
        magpost = "_ace_mag_1m.txt"
        swepost = "_ace_swepam_1m.txt"
        dday = (t1-t0).days
        current_time = t0
        for i in range(dday):
            current_time = t0 + datetime.timedelta(days=i)
            datestr = current_time.strftime('%Y%m%d')
            magstr = swpcdir + datestr + magpost
            swestr = swpcdir + datestr + swepost
            print("Downloading ACE data: ",datestr + magpost)
            print("Downloading ACE data: ",datestr + swepost)
            urllib.request.urlretrieve(magstr,datestr + magpost)
            urllib.request.urlretrieve(swestr,datestr + swepost)
            self.__readACE(datestr,magpost,swepost)

    def __readACE(self,datestr,magpost,swepost):
        rdfile = open(datestr+magpost,'r')
        text = rdfile.readlines()
        for i,j in enumerate(text):
            if j[0] == '2' or j[0] == '1':
                magskip = i
                break
        rdfile.close()
        rdfile = open(datestr+swepost,'r')
        text = rdfile.readlines()
        for i,j in enumerate(text):
            if j[0] == '2' or j[0] == '1':
                sweskip = i
                break
        rdfile.close()

        dat = np.genfromtxt(datestr + magpost,skip_header=magskip, autostrip=True,dtype=None)
        magtime = []
        bx = []
        by = []
        bz = []
        for i in dat:
            currenttime = datetime.datetime(i[0],i[1],i[2],i[3],i[5])
            print(currenttime)
            if currenttime >= t0 and currenttime <= t1:
                magtime.append(currenttime)
                bx.append(i[7])
                bx.append(i[8])
                bx.append(i[9])
        dat = np.genfromtxt(datestr + swepost,skip_header=sweskip, autostrip=True,dtype=None)
        swetime = []
        n = []
        v = []
        t = []
        for i in dat:
            currenttime = datetime.datetime(i[0],i[1],i[2],i[3],i[5])
            print(currenttime)
            if currenttime >= t0 and currenttime <= t1:
                magtime.append(currenttime)
                bx.append(i[7])
                bx.append(i[8])
                bx.append(i[9])

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
