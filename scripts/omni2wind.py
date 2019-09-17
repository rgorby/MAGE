#!/usr/bin/env python

#Converts OMNI output data to Gamera solar wind file to be used as boundary conditions

#Reads from ASCII file
#Time(min) Density (AMU/cm^-3) Vx(km/s) Vy(km/s) Vz(km/s) Cs(km/s) Bx(nT) By(nT) Bz(nT) B(nT) tilt(rad)

#Writes to HDF5 Gamera wind file
#t,D,V,P,B = [s],[#/cm3],[m/s],[nPa],[nT]

#Utilizes ai.cdas and geopack, make sure to install modules before running. For more info go to https://bitbucket.org/aplkaiju/kaiju/wiki/Gamerasphere

Mp = 1.67e-27 #Proton mass [kg]
gamma = 5/3.0

import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import kaipy.solarWind
from  kaipy.solarWind import swBCplots
from  kaipy.solarWind.OMNI import OMNI
import datetime
from ai import cdas

def bxFit(sw, fileType, filename):
    def bxFitPlot(bxFit_array):
        kaipy.solarWind.swBCplots.BasicPlot(sw.data, 'time_doy', 'bx', color='k')
        plt.plot(sw.data.getData('time_doy'), bxFit_array, 'g')
        plt.title('Bx Fit Coefficients ('+fileType+'):\n$Bx_{fit}(0)$=%f      $By_{coef}$=%f      $Bz_{coef}$=%f' % (coef[0], coef[1], coef[2]) )
        plt.legend(('$Bx$','$Bx_{fit}$'))

    coef = sw.bxFit()

    print('Bx Fit Coefficients are ', coef)
    by = sw.data.getData('by')
    bz = sw.data.getData('bz')
    bxFit = coef[0] + coef[1] * by + coef[2] * bz

    # Save plot
    bxFitPlot(bxFit)
    bxPlotFilename = os.path.basename(filename) + '_bxFit.png'
    print('Saving "%s"' % bxPlotFilename)
    plt.savefig(bxPlotFilename)

if __name__ == "__main__":
        fOut = "bcwind.h5"
        mod = "LFM"
        t0="1/1/2010, 00:00:00"
        t1="1/1/2010, 00:59:59"
        Ts = 0.0
        obs="OMNI"
        MainS = """ This script does several things:
                      1. Fetch OMNI data from CDAWeb between the specified times 
                      2. Generate standard plots of solar wind data
                      3. Write output in a model file format.
                         - "LFM" format will:
                             a. Generate coefficients for Bx Fit
                             b. Save a bcwind.h5 file
                         - "TIEGCM" format will:
                             a. Compute 15-minute boxcar average lagged by 5 minutes
                             b. Sub-sample at 5-minutes
                             c. Write NetCDF IMF data file
        """

        parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
        parser.add_argument('-t0',type=str,metavar="TStart",default=t0,help="Start time in a string 'MM/DD/YYY, HH:MM:SS' (default: %(default)s)")
        parser.add_argument('-t1',type=str,metavar="TStop",default=t1,help="End time in a string'MM/DD/YYY, HH:MM:SS' (default: %(default)s)")
        parser.add_argument('-obs',type=str,metavar="OMNI",default=obs,help="Select spacecraft to obtain observations from (default: %(default)s)")
        parser.add_argument('-o',type=str,metavar="wind.h5",default=fOut,help="Output Gamera wind file (default: %(default)s)")
        parser.add_argument('-m',type=str,metavar="LFM",default=mod,help="Format to write.  Options are LFM or TIEGCM (default: %(default)s)")
        parser.add_argument('-TsG',type=float,metavar="TStart",default=Ts,help="Gamera start time [min] (default: %(default)s)")
        parser.add_argument('-TsL',type=float,metavar="TStart",default=Ts,help="LFM start time [min] (default: %(default)s)")
        parser.add_argument('-nobx', action='store_true',default=False,help="Zero out Bx (default: %(default)s)")
        #Finalize parsing
        args = parser.parse_args()

        fOut = args.o
        mod = args.m
        TsG = args.TsG
        TsL = args.TsL
        nobx = args.nobx

        t0 = args.t0
        t1 = args.t1

        fmt='%m/%d/%Y, %H:%M:%S'

        if (obs == 'OMNI'):
            fileType = 'OMNI'
            filename = 'OMNI_HRO_1MIN.txt'
            
            #obtain 1 minute resolution observations from OMNI dataset
            fIn = cdas.get_data(
               'sp_phys',
               'OMNI_HRO_1MIN',
               datetime.datetime.strptime(t0,fmt),
               datetime.datetime.strptime(t1,fmt),
               ['BX_GSE,BY_GSE,BZ_GSE,Vx,Vy,Vz,proton_density,T']
            )

        else:
            raise Exception('Error:  Not able to obtain dataset from spacecraft. Please select another mission.')

        # Read the solar wind data into 'sw' object and interpolate over the bad data.
        sw = eval('kaipy.solarWind.'+fileType+'.'+fileType)(fIn)

        # Do output format-specific tasks:
        if (mod == 'TIEGCM'):
            # Write TIEGCM IMF solar wind file
            #FIXME: need to update when want to include, example code in pyLTR.SolarWind.Writer.TIEGCM
            raise Exception('Error:  Cannot currently produce TIEGCM output.')
        elif (mod == 'LFM'):
            # Bx Fit 
            bxFit( sw, fileType, filename)
        
            # Interpolate to one minute:
            time_1minute = range(int(sw.data.getData('time_min').min()),
                                 int(sw.data.getData('time_min').max()) )
            n    = np.interp(time_1minute, sw.data.getData('time_min'), sw.data.getData('n'))
            vx   = np.interp(time_1minute, sw.data.getData('time_min'), sw.data.getData('vx'))
            vy   = np.interp(time_1minute, sw.data.getData('time_min'), sw.data.getData('vy'))
            vz   = np.interp(time_1minute, sw.data.getData('time_min'), sw.data.getData('vz'))
            cs   = np.interp(time_1minute, sw.data.getData('time_min'), sw.data.getData('cs'))
            bx   = np.interp(time_1minute, sw.data.getData('time_min'), sw.data.getData('bx'))
            by   = np.interp(time_1minute, sw.data.getData('time_min'), sw.data.getData('by'))
            bz   = np.interp(time_1minute, sw.data.getData('time_min'), sw.data.getData('bz'))
            b    = np.interp(time_1minute, sw.data.getData('time_min'), sw.data.getData('b'))

            #initalize matrix to hold solar wind data
            lfmD = np.zeros((n.shape[0],11))

            date = sw.data.getData('meta')['Start date']
            for i,time in enumerate(time_1minute):
                # Convert relevant quantities to SM Coordinates
                v_sm = sw._gsm2sm(date+datetime.timedelta(minutes=time), vx[i],vy[i],vz[i])
                b_sm = sw._gsm2sm(date+datetime.timedelta(minutes=time), bx[i],by[i],bz[i])
                tilt = sw._getTiltAngle(date+datetime.timedelta(minutes=time))

                lfmD[i] = [time,n[i],v_sm[0],v_sm[1],v_sm[2],cs[i],b_sm[0],b_sm[1],b_sm[2],b[i],tilt]

            # Save a plot of the solar wind data.
            kaipy.solarWind.swBCplots.MultiPlot(sw.data, 'time_doy', ['n', 'vx','vy','vz','t','bx','by','bz'])
            plt.title('Solar Wind data for\n %s' % filename)
            swPlotFilename = os.path.basename(filename) + '.png'
            print('Saving "%s"' % swPlotFilename)
            plt.savefig(swPlotFilename)

            print("Converting to Gamera solar wind file")
            Nt,Nv = lfmD.shape
            print("\tFound %d variables and %d lines"%(Nv,Nt))
            if (nobx):
                print("\tNot using Bx fields")
            #Create holders for Gamera data
            T  = np.zeros(Nt)
            D  = np.zeros(Nt)
            P  = np.zeros(Nt)
            Vx = np.zeros(Nt)
            Vy = np.zeros(Nt)
            Vz = np.zeros(Nt)
            Bx = np.zeros(Nt)
            By = np.zeros(Nt)
            Bz = np.zeros(Nt)
            ThT= np.zeros(Nt) #Tilt

            #Convert LFM time to seconds and reset to start at 0
            print("\tOffsetting from LFM start (%5.2f min) to Gamera start (%5.2f min)"%(TsL,TsG))
            T0 = lfmD[:,0].min()
            T = (lfmD[:,0]-TsL+TsG)*60
            
            #Calculating time in UT
            UT = []
            [UT.append(np.string_(date+datetime.timedelta(seconds=i)).strip()) for i in T]

            #Density, magnetic field, and tilt don't require scaling
            D   = lfmD[:,1]
            ThT = lfmD[:,10]
            if (nobx):
                Bx[:] = 0.0
            else:
                Bx  = lfmD[:,6]
            By  = lfmD[:,7]
            Bz  = lfmD[:,8]

            #Velocity
            vScl = 1.0e+3 #km/s->m/s
            Vx  = vScl*lfmD[:,2]
            Vy  = vScl*lfmD[:,3]
            Vz  = vScl*lfmD[:,4]

            #Now get pressure (nPa) from sound speed (km/s)
            #Convert density: AMU/cm3 -> kg/m3
            Cs = lfmD[:,5]
            Dmks = (D*Mp)*( (1.0e+2)**3.0 ) #kg/m3
            Cmks = Cs*1.0e+3 #m/s
            P = (1.0e+9)*Dmks*Cmks*Cmks/gamma #nPa

            print("Writing Gamera solar wind to %s"%(fOut))
            with h5py.File(fOut,'w') as hf:
                hf.create_dataset("T" ,data=T)
                hf.create_dataset("UT",data=UT)
                hf.create_dataset("D" ,data=D)
                hf.create_dataset("P" ,data=P)
                hf.create_dataset("Vx",data=Vx)
                hf.create_dataset("Vy",data=Vy)
                hf.create_dataset("Vz",data=Vz)
                hf.create_dataset("Bx",data=Bx)
                hf.create_dataset("By",data=By)
                hf.create_dataset("Bz",data=Bz)
                hf.create_dataset("tilt",data=ThT)
                
        else:
            raise Exception('Error:  Misunderstood output file format.')
