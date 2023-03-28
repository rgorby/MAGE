# Functions to calculate Supermag indices from simulated outputs
import os
import numpy as np
import h5py
import kaipy.kaiH5 as kh5
from astropy.time import Time
from scipy.spatial import qhull
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import warnings
import math
import datetime
import json

#### NEED TO POINT TO SUPERMAG API SCRIPT
#### /glade/p/hao/msphere/gamshare/supermag/supermag_api.py 
import kaipy.supermag_api as smapi
####

# this warning is very annoying
warnings.filterwarnings("ignore", category=DeprecationWarning) 

###############################################################################
###############################################################################

def Time2Float(x):
    """Converts datetime to float, so that interpolation/smoothing can be performed"""
    if (type(x) == np.ndarray) or (type(x) == list):
        emptyarray = []
        for i in x:
            z = (i - datetime.datetime(1970, 1, 1, 0)).total_seconds()
            emptyarray.append(z)
        emptyarray = np.array([emptyarray])
        return emptyarray[0]
    else:
        return (x - datetime.datetime(1970, 1, 1, 0)).total_seconds()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def Float2Time(x):
    """Converts array back to datetime so that it can be plotted with time on the axis"""
    if (type(x) == np.ndarray) or (type(x) == list):
        emptyarray = []
        for i in x:
            z = datetime.datetime.utcfromtimestamp(i)
            emptyarray.append(z)
        emptyarray = np.array([emptyarray])
        return emptyarray[0]
    else:
        return datetime.datetime.utcfromtimestamp(x)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def interp_tri(xy):
    """
    do first part of griddata calculation so it can be done only once for the
    main grid rather than being called every time we need to do the
    interpolation for each path (this speeds the code up significantly)
    use interp_grid for the second part
    from stackoverflow.com/questions/20915502/
    """
    tri = qhull.Delaunay(xy)
    return tri

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def interp_grid(values, tri, uv, d=2):
    """
    grid interpolation - done for each connection
    This is esentially what griddata does but with the first part of the
    triangulation separated out (in interp_tri) so it can be performed just
    once for the E-field grid. Then this part is for finding the interpolated
    values for all the paths along each of the connections
    Doing it this way gave ~7x speed improvement for a uniform field
    Inputs
    ------
    values : np.array
       values for each point on the grid - i.e. the E-field values
    tri : scipy.spatial.qhull.Delaunay
       from interp_grid
    uv : np.array
       new grid positions for interpolation - i.e. pathsteps along the
       connection
    d : scalar, default = 2
       ?? dimensions? would be 3 for a 3d grid
    """
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
    return np.einsum('nj,nj->n', np.take(values, vertices),
                     np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True))))

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def FetchSMData(user, start, numofdays, savefolder, badfrac=0.1, nanflags=True, doDB=True):
    """Retrieve all available SuperMagnet data for a specified period
    If data has not already been downloaded, fetches data from Supermag
    
    Parameters
    -----------
    user = username for downloading SuperMag Data
    start = start day (datetime obj)
    numofdays = number of days from start to download
    savefolder = folder where downloaded data will be saved as json. This function
        looks here first for saved data before downloading.
    badfrac = tolerable fraction of data that is 99999.0. Sites with more bad data
        than this fraction will be ignored
    nanflags = will set 99999.0 values to nans if True (True by default)    
    doDB = Whether to pull the pre-baselined values from supermag

    Returns
    -----------
    Dictionary which has the following data as keys:
    {td, sitenames, glon, glat, mlon, mlat, mcolat, 
                BNm, BEm, BZm, BNg, BEg, BZg, MLT, DECL, SZA}
    """

    if (doDB):
        smFlags = "all,delta=start,baseline=all"
    else:
        smFlags = "all"
        
    # Look at all saved .jsons
    filenames = [x for x in sorted(os.listdir(savefolder)) if '.json' in x]
    startstr = start.strftime('%Y-%m-%d')

    exists = False
    for filename in filenames:
        if startstr in filename:
            # first split grads numberOfDays.json and second returns the numberOfDays
            daystring = filename.split('_')[-1].split('.')[0] 

            if int(daystring) >= numofdays:
                print("Supermag data already exists locally")
                print(filename)
 
                exists = True
                break
            continue

    if exists == False:
        print("Supermag data not local, fetching:")

        STATUS, master, badindex = [], [], []

        #ZZZ
        status, stations = smapi.SuperMAGGetInventory(user, startstr, extent = 86400*numofdays)
        for iii in stations:
            print("Fetching: ", iii)
            #ZZZ
            
            status, A = smapi.SuperMAGGetData(user, startstr, extent=86400*numofdays, 
                                           flagstring=smFlags, station = iii, FORMAT = 'list')
            
            if status:
                quickvals = np.array([x['N']['nez'] for x in A])

                # get rid of data if too many bad values
                if np.sum(quickvals>999990.0) >= badfrac*len(quickvals):
                    badindex.append(False)
                    print(iii, "BAD")
                else:
                    badindex.append(True)

                STATUS.append(status)
                master.append(A)
            else:
                STATUS.append(status)
                badindex.append(False)
                master.append(['BAD'])

        badindex = np.array(badindex)
        master, stations = np.array(master)[badindex], np.array(stations)[badindex]

        # Make the Supermag data a dict for saving later
        output = {}
        for i in master:
            output[i[0]['iaga']] = list(i)

        filename = "SM_DATA_" + startstr + '_%1d.json' % (numofdays)
        with open(os.path.join(savefolder,filename), mode='w') as f:
            #print(savefolder + filename)
            json.dump(output, f)
        f.close()

    # Now read in the data
    with open(os.path.join(savefolder,filename),mode='r') as r:
        #print(savefolder + filename)
        rr = json.load(r)

    sitenames = np.array(list(rr.keys()))

    timedate = np.array([Float2Time(x['tval']) for x in rr[sitenames[0]]])
    glon = np.array([rr[x][0]['glon'] for x in sitenames])
    glon[glon>180] -= 360 
    glat = np.array([rr[x][0]['glat'] for x in sitenames])
    mlon = np.array([rr[x][0]['mlon'] for x in sitenames])
    mlat = np.array([rr[x][0]['mlat'] for x in sitenames])
    mcolat = np.array([rr[x][0]['mcolat'] for x in sitenames])

    BNm = np.zeros((len(timedate), len(sitenames)))

    BEm, BZm, BNg, BEg, BZg = np.copy(BNm), np.copy(BNm), np.copy(BNm), np.copy(BNm), np.copy(BNm)
    MLT, DECL, SZA = np.copy(BNm), np.copy(BNm), np.copy(BNm)

    for j, v2 in enumerate(sitenames):
        for i, v1 in enumerate(rr[v2]):
            BNm[i][j] = v1['N']['nez']
            BEm[i][j] = v1['E']['nez']
            BZm[i][j] = v1['Z']['nez']
            
            BNg[i][j] = v1['N']['geo']
            BEg[i][j] = v1['E']['geo']
            BZg[i][j] = v1['Z']['geo']

            MLT[i][j] = v1['mlt']
            DECL[i][j] = v1['decl']
            SZA[i][j] = v1['sza']

    if (nanflags == True):
        BNm[BNm==999999.0] = np.nan
        BEm[BEm==999999.0] = np.nan
        BZm[BZm==999999.0] = np.nan

        BNg[BNg==999999.0] = np.nan
        BEg[BEg==999999.0] = np.nan
        BZg[BZg==999999.0] = np.nan

    # only use points after start of sim data
    i = (timedate>=start)
    output = {'td':timedate[i], 'sitenames':sitenames, 'glon':glon, 'glat':glat,
             'mlon':mlon, 'mlat':mlat, 'mcolat':mcolat, 'BNm':BNm[i], 'BEm':BEm[i], 'BZm':BZm[i],
             'BNg':BNg[i], 'BEg':BEg[i], 'BZg':BZg[i], 'mlt':MLT[i], 'decl':DECL[i], 'sza':SZA[i]}     
    return output

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def MJD2Str(m0):
    """Returns timestrings and datetime objects for simulation float"""
    dtObj = Time(m0,format='mjd').datetime
    tStr = dtObj.strftime("%H:%M:%S") + " " +  dtObj.strftime("%m/%d/%Y")
    return tStr, dtObj

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def ReadSimData(filename, quiet = True):
    """Read in all of the needed ground mag data from a .h5 file
    
    Parameters
    -----------
    filename = address of .h5 delta-b object
    
    Output
    -----------
    dictionary of:
    td = time for each step
    glon, glat = geo coordinates
    mlt = magnetic local time
    smlon,smlat = sm coordinates
    dBt, dBp, dBr = Btheta (*-1), phi, radius
    dBn = magnetic North component
    """
    f = h5py.File(filename, 'r')

    # Assume that the XYZ0 shell is the R=1.0 shell
#    X0, X1 = np.array(f['Xcc'])[0], np.array(f['Xcc'])[1]
#    Y0, Y1 = np.array(f['Ycc'])[0], np.array(f['Ycc'])[1]
#    Z0, Z1 = np.array(f['Zcc'])[0], np.array(f['Zcc'])[1]

    # get the timesteps from the simulation
    nSteps,sIds = kh5.cntSteps(filename)
    MJDs = kh5.getTs(filename,sIds,"MJD",-np.inf)
    TIMETHING = np.array([MJD2Str(x) for x in MJDs])
    tStr, td = TIMETHING[:,0], TIMETHING[:,1]

    # Get lat lon
    thetacc = np.rad2deg(np.array(f['Thetacc'])[0].ravel())
    lat = 90 - thetacc
    phicc   = np.rad2deg(np.array(f['Phicc'])[0].ravel())
    lon = np.copy(phicc)
    lon[lon>180] -= 360


    first_step_name = "Step#%s" % sIds[0]
    data = f[first_step_name]
    dBt, dBp = np.array(data['dBt'])[0], np.array(data['dBp'])[0]
    BH = np.sqrt(dBt**2 + dBp**2)
    maglat = np.array(data['smlat'])[0].flatten()

    # Create the arrays to populate with data
    masterdBn = np.zeros((nSteps, len(dBt.flatten())))
    masterdBt, masterdBp, masterdBr = np.copy(masterdBn), np.copy(masterdBn), np.copy(masterdBn)
    mastersmlon, masterMLT = np.copy(masterdBn), np.copy(masterdBn)

    # make big array of all of the data
    for i in sIds:
        data = f['Step#%01d' % i]
        dBt = np.array(data['dBt'])[0].flatten()
        dBp = np.array(data['dBp'])[0].flatten()
        dBr = np.array(data['dBr'])[0].flatten()
        dBn = np.array(data['dBn'])[0].flatten()        # change back for STPAT!

        masterdBt[i] = dBt*-1 # Bx (mag coords) points South, opposite of SM Z
        masterdBp[i] = dBp
        masterdBr[i] = dBr
        masterdBn[i] = dBn

        # Get MLT
        smlon = np.array(data['smlon'])[0].flatten()

        smlon[smlon>180] -= 360 # have smlon between -180 and 180
        mastersmlon[i] = smlon
        UT = td[i].hour + (td[i].minute/60.)
        MLT = ((np.deg2rad(smlon)/np.pi) * 12) + 12
        masterMLT[i] = MLT

        if quiet == False:
            print("Reading in: ", i)

    masterdBh = np.sqrt(masterdBt**2 + masterdBp**2)
    output = {'td':td, 'glon':lon, 'glat':lat, 'mlat':maglat, 'mlt':masterMLT, 'smlon':mastersmlon,
              'dBt':masterdBt, 'dBp':masterdBp, 'dBr':masterdBr, 'dBn':masterdBn}
    
    return output

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def FetchSMIndices(user, start, numofdays, wanted = 'ALL'):
    """Retrieve SME, SML, SMU indices from SuperMag

    Parameters
    -----------
    user = username for downloading SuperMag Data
    start = start day (datetime obj)
    numofdays = number of days from start to download
    wanted = list of wanted attrs (e.g., ['SME', 'SML']
        downloads all by default

    Returns
    -----------
    output = dictionary of wanted values as arrays + 'td' array
    """
    #ZZZ
    status, vals = smapi.SuperMAGGetIndices(user, start, 86400*numofdays, 'all', FORMAT='list')

    if (wanted == 'ALL'):
        wanted = list(vals[0].keys())[1:]

    output = {x:[] for x in wanted}
    output['td'] = []

    for step in vals:
        output['td'].append(Float2Time(step['tval']))

        for j in wanted:
            output[j].append(step[j])

    ind = np.array(output['td']) >= start

    for i in output.keys():
        output[i] = np.array(output[i])[ind]

    return output

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def CalculateSMRBins(mlt, B, mlat):
    """Calculate SMR and bins for a given mlt, B, and mlat
    Does the latitudinal scaling

    Parameters
    -----------
    mlt = array of magnetic local time
    B = array of B
    mlat = array of magnetic latitudes

    Returns
    -----------
    SMR = calculated SMR (SYMH equivalent) as array
    SMR00, SMR06, SMR12, SMR18 = 6-hour SMR bins as arrays
    """
    latscaling = np.array([np.cos(np.deg2rad(mlat))])
    Bnorth = np.copy(B)/latscaling

    ind = np.abs(mlat) <= 50
    BNorthSMR = Bnorth.T[ind].T
    mltSMR    = np.copy(mlt).T[ind].T

    ind1 = (mltSMR >= 21) + (mltSMR <  3)      
    ind2 = (mltSMR >= 3)  * (mltSMR < 9)
    ind3 = (mltSMR >= 9)  * (mltSMR < 15)
    ind4 = (mltSMR >= 15) * (mltSMR < 21)

    B1, B2, B3, B4, SMR = [], [], [], [], []
    for i in range(BNorthSMR.shape[0]):

        bin1 = np.nanmean(BNorthSMR[i][ind1[i]])
        bin2 = np.nanmean(BNorthSMR[i][ind2[i]])
        bin3 = np.nanmean(BNorthSMR[i][ind3[i]])
        bin4 = np.nanmean(BNorthSMR[i][ind4[i]])

        B1.append(bin1)
        B2.append(bin2)
        B3.append(bin3)
        B4.append(bin4)
        
        SMR.append(bin1 + bin2 + bin3 + bin4)

    SMR, B1, B2, B3, B4 = np.array(SMR)/4., np.array(B1), np.array(B2), np.array(B3), np.array(B4)

    return SMR, B1, B2, B3, B4

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def InterpolateSimData(SIM, SM):
    """Interpolates simulated data to SuperMag station coordinates, and calculates 
    different indices (SME/U/L/R).

    This function first interpolates sim data so that it is on the same timestamps
    as the SM data, then interpolates to the same geographic positions. Then it calculates 
    SME/U/L/R for the SM interpolated data, + all simulated data (superSME/U/L/R).
    
    Parameters
    -----------
    SIM = dictionary of simulated data (from ReadSimData())
    SM  = dictionary of supermag data (from FetchSMData())
    
    Output
    -----------
    dictionary of 
    td = timestamps
    glon, glat, mlon, mlat = Supermag coordinates
    dBn, dBt, dBp = interpolated magnetic components
    mlt = interpolated mlt

    dBnsmall, mltsmall = simulated data for all points, but limited to Supermag timestamps

    SME/U/L/R = indices at interpolated points
    superSME/U/L/R = indices using all data
    
    SMRbins = list of 4 bins using interpolated data(SMR00/06/12/18)
    superSMRbins = as before with all data
    """

    # find limits of overlapping data
    starttimei = np.max([SIM['td'][0], SM['td'][0]])
    endtimei = np.min([SIM['td'][-1], SM['td'][-1]])

    index2 = (SM['td'] >= starttimei) * (SM['td'] <= endtimei)   # for real SM data
    smalltd = SM['td'][index2]

    # convert to float for initial time-series interpolation
    SIMtdf = Time2Float(SIM['td'])
    SMtdf = Time2Float(smalltd)

    # interpolate the simulated values so that the times match real timestamps (i.e., on each minute exactly)
    newdBn = np.zeros((len(SMtdf), SIM['dBn'].shape[1]))
    newmlt, newdBt, newdBp = np.copy(newdBn), np.copy(newdBn), np.copy(newdBn)

    for i,v in enumerate(SIM['dBn'].T):
        newdBn[:,i] = np.interp(SMtdf, SIMtdf, v)
        newmlt[:,i] = np.interp(SMtdf, SIMtdf, SIM['mlt'].T[i])
        newdBt[:,i] = np.interp(SMtdf, SIMtdf, SIM['dBt'].T[i])
        newdBp[:,i] = np.interp(SMtdf, SIMtdf, SIM['dBp'].T[i])

    # now interpolate to the SM coordinates
    sim_points = np.vstack((SIM['glon'], SIM['glat'])).T
    wanted_points = np.vstack((SM['glon'], SM['glat'])).T
    interptri = interp_tri(sim_points)

    # make empty arrays: len(smalltd) steps * len(SM['glon']) sites
    interp_dBt = np.zeros((len(smalltd), len(SM['glon'])))
    interp_dBn, interp_dBp = np.copy(interp_dBt), np.copy(interp_dBt)
    #interp_dBp, interp_dBr, interp_dBn = np.copy(interp_dBt), np.copy(interp_dBt), np.copy(interp_dBt)
    interp_smlon, interp_mlt = np.copy(interp_dBt), np.copy(interp_dBt)

    for i in range(len(smalltd)):
        interp_dBn[i] = interp_grid(newdBn[i], interptri, wanted_points)
        interp_mlt[i] = interp_grid(newmlt[i], interptri, wanted_points)
        interp_dBt[i] = interp_grid(newdBt[i], interptri, wanted_points)
        interp_dBp[i] = interp_grid(newdBp[i], interptri, wanted_points)

        if (i%10) == 0:
            pass
            #print(i)

    ##### Now to calculate indices #####
    # calculate SME, SML, SMU
    ind = (SM['mlat'] >= 40) * (SM['mlat'] <= 80) # equivalent to SME/SMU/SML
    Northonly = interp_dBn.T[ind].T
    SMU_calc = np.max(Northonly, axis = 1)
    SML_calc = np.min(Northonly, axis = 1)
    SME_calc = SMU_calc - SML_calc

    # calculate SME indices with ALL points between 40 and 80 (supersupermag)
    ind = (SIM['mlat'] >= 40) * (SIM['mlat'] <= 80)
    Northonly = newdBn.T[ind].T
    SMU_calc2 = np.max(Northonly, axis = 1)
    SML_calc2 = np.min(Northonly, axis = 1)
    SME_calc2 = SMU_calc2 - SML_calc2

    # calculate SMR and SMR00, SMR06 etc.
    # interpolated data first
    SMRi, SMR00i, SMR06i, SMR12i, SMR18i = CalculateSMRBins(interp_mlt, interp_dBn, SM['mlat'])
    SMR, SMR00, SMR06, SMR12, SMR18 = CalculateSMRBins(newmlt, newdBn, SIM['glat'])

    output = {'td':smalltd, 'glon':SM['glon'], 'glat':SM['glat'], 'mlon':SM['mlon'], 'mlat':SM['mlat'],
              'dBn':interp_dBn, 'dBt':interp_dBt, 'dBp':interp_dBp, 'mlt':interp_mlt, 
              'dBnsmall':newdBn, 'mltsmall':newmlt, 'SMU':SMU_calc, 
              'SME':SME_calc, 'SML':SML_calc, 'superSME':SME_calc2, 'superSML':SML_calc2, 'superSMU':SMU_calc2,
              'SMR':SMRi, 'SMRbins':[SMR00i, SMR06i, SMR12i, SMR18i], 'superSMR':SMR, 
              'superSMRbins':[SMR00, SMR06, SMR12, SMR18]}

    return output

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
 
def SMContourPlotPrep(mlt, td, bx):
    """Bins max and min Bn to make hourly SMU/SML contour plots 
    
    Parameters
    -----------
    mlt = magnetic local time array
    td  = array of datetime
    bx = array of Bx (magnetic north component)
    
    Output:
    -----------
    hourbins = array of 1.5*hour (24.5, 23.5...)
    SMUbig = array of binned SMU values
    SMLbig = array of binned SML values
    """

    SMUbig = np.zeros((24, len(td)))    # 24 3-hour bins, for each min
    SMLbig = np.zeros((24, len(td)))    # 24 3-hour bins, for each min

    hourbins = np.arange(24, 0, -1) + 0.5

    for i, hour in enumerate(hourbins):
        #print(i)
        diff1 = np.abs(mlt - hour) <= 1.5
        diff2 = np.abs(mlt + 24 - hour) <= 1.5
        diff3 = np.abs(mlt - 24 - hour) <= 1.5

        ind = diff1 + diff2 + diff3
        bx2 = np.copy(bx)
        bx2[np.invert(ind)] = np.nan

        SMUbig[i] = np.nanmax(bx2, axis = 1)
        SMLbig[i] = np.nanmin(bx2, axis = 1)

    return hourbins, SMUbig, SMLbig

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def MakeContourPlots(SM, SMinterp, maxx = 'default', fignumber = 1):
    """Makes hourly contour plots for SMU/SML
    
    Parameters
    -----------
    SM, SMinterp = SM and interpolated dictionaries
    maxx = absolute max value for colorbar. If left as 'default', will adapt to data
    fignumber = number of figure

    Output:
    -----------
    Tasty contour plot
    """

    # Real data
    hourbins, SMUbigreal, SMLbigreal = SMContourPlotPrep(SM['mlt'], SM['td'], SM['BNm'])

    # Reduced Sim data 
    hourbins, SMUbigsimR, SMLbigsimR = SMContourPlotPrep(SMinterp['mltsmall'], SMinterp['td'], SMinterp['dBnsmall'])

    # Sim data at SM sites
    hourbins, SMUbigsim, SMLbigsim = SMContourPlotPrep(SMinterp['mlt'], SMinterp['td'], SMinterp['dBn'])

    if maxx == 'default':
        max1 = np.max([np.max(SMUbigreal), np.max(SMLbigreal), np.max(SMUbigsimR), np.max(SMLbigsimR), np.max(SMUbigsimR), np.max(SMLbigsimR)])
        min1 = np.min([np.min(SMUbigreal), np.min(SMLbigreal), np.min(SMUbigsimR), np.min(SMLbigsimR), np.min(SMUbigsimR), np.min(SMLbigsimR)])

        maxx = np.max([max1, np.abs(min1)])

    xxx, yyy, qqq = 0.5, 0.9, 16
    cmapp = plt.cm.RdBu
      
    #######################
    fig = plt.figure(fignumber)
    plt.clf()

    ax1 = plt.subplot(3, 2, 1)
    plt.pcolormesh(SM['td'], hourbins, SMUbigreal, vmin = -1*maxx, vmax = maxx, cmap = cmapp, shading = 'auto')
    plt.xlim([SMinterp['td'][0], SMinterp['td'][-1]])
    plt.grid(True)
    plt.text(xxx, yyy, "Real SMU", horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes, fontsize = qqq)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    ax2 = plt.subplot(3, 2, 2)
    plt.pcolormesh(SM['td'], hourbins, SMLbigreal, vmin = -1*maxx, vmax = maxx, cmap = cmapp, shading = 'auto')
    plt.xlim([SMinterp['td'][0], SMinterp['td'][-1]])
    plt.grid(True)
    plt.text(xxx, yyy, "Real SML", horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes, fontsize = qqq)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    #######################

    ax3 = plt.subplot(3, 2, 3)
    plt.pcolormesh(SMinterp['td'], hourbins, SMUbigsim, vmin = -1*maxx, vmax = maxx, cmap = cmapp, shading = 'auto')
    plt.xlim([SMinterp['td'][0], SMinterp['td'][-1]])
    plt.grid(True)
    plt.text(xxx, yyy, "Interpolated SMU", horizontalalignment='center', verticalalignment='center', transform=ax3.transAxes, fontsize = qqq)
    plt.ylabel("MLT", fontsize = qqq)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    ax4 = plt.subplot(3, 2, 4)
    plt.pcolormesh(SMinterp['td'], hourbins, SMLbigsim, vmin = -1*maxx, vmax = maxx, cmap = cmapp, shading = 'auto')
    plt.xlim([SMinterp['td'][0], SMinterp['td'][-1]])
    plt.grid(True)
    plt.text(xxx, yyy, "Interpolated SML", horizontalalignment='center', verticalalignment='center', transform=ax4.transAxes, fontsize = qqq)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    #######################
       
    ax5 = plt.subplot(3, 2, 5)
    plt.pcolormesh(SMinterp['td'], hourbins, SMUbigsimR, vmin = -1*maxx, vmax = maxx, cmap = cmapp, shading = 'auto')
    plt.xlim([SMinterp['td'][0], SMinterp['td'][-1]])
    plt.grid(True)
    plt.text(xxx, yyy, "Super-SMU", horizontalalignment='center', verticalalignment='center', transform=ax5.transAxes, fontsize = qqq)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    ax6 = plt.subplot(3, 2, 6)
    im = plt.pcolormesh(SMinterp['td'], hourbins, SMLbigsimR, vmin = -1*maxx, vmax = maxx, cmap = cmapp, shading = 'auto')
    plt.xlim([SMinterp['td'][0], SMinterp['td'][-1]])
    plt.grid(True)
    plt.text(xxx, yyy, "Super-SML", horizontalalignment='center', verticalalignment='center', transform=ax6.transAxes, fontsize = qqq)
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))

    plt.tight_layout()

    fig.subplots_adjust(right=0.85)
    cbar_ax = fig.add_axes([0.875, 0.05, 0.03, 0.90])
    cb = fig.colorbar(im, cax = cbar_ax)
    cb.set_label("Mag Perturbation (nT)", fontsize = qqq)
    plt.show()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def MakeIndicesPlot(SMI, SMinterp, fignumber = 1):
    fig = plt.figure(fignumber)
    plt.clf()

    ax1 = plt.subplot(311)
    plt.plot(SMinterp['td'], SMinterp['SME'], 'r', label = 'interpolated')
    plt.plot(SMinterp['td'], SMinterp['SMU'], 'r')
    plt.plot(SMinterp['td'], SMinterp['SML'], 'r')

    plt.plot(SMinterp['td'], SMinterp['superSME'], 'g', label = 'super')
    plt.plot(SMinterp['td'], SMinterp['superSMU'], 'g')
    plt.plot(SMinterp['td'], SMinterp['superSML'], 'g')

    plt.plot(SMI['td'], SMI['SME'], 'b', label = 'real')
    plt.plot(SMI['td'], SMI['SMU'], 'b')
    plt.plot(SMI['td'], SMI['SML'], 'b')
    plt.grid(True)
    plt.legend()
    plt.xlim([SMinterp['td'][0], SMinterp['td'][-1]])
    plt.ylabel("SME/U/L", fontsize = 20)

    ax2 = plt.subplot(312)
    plt.plot(SMinterp['td'], SMinterp['SMR'], 'r', label = 'interpolated')
    plt.plot(SMinterp['td'], SMinterp['superSMR'], 'g', label = 'super')
    plt.plot(SMI['td'], SMI['smr'], 'b', label = 'real')
    plt.grid(True)
    plt.legend()
    plt.xlim([SMinterp['td'][0], SMinterp['td'][-1]])
    plt.ylabel("SMR", fontsize = 20)

    ax3 = plt.subplot(313)
    plt.plot(SMinterp['td'], SMinterp['SMRbins'][0], 'r', label = 'interpolated')
    plt.plot(SMinterp['td'], SMinterp['SMRbins'][1], 'r')
    plt.plot(SMinterp['td'], SMinterp['SMRbins'][2], 'r')
    plt.plot(SMinterp['td'], SMinterp['SMRbins'][3], 'r')

    plt.plot(SMinterp['td'], SMinterp['superSMRbins'][0], 'g', label = 'super')
    plt.plot(SMinterp['td'], SMinterp['superSMRbins'][1], 'g')
    plt.plot(SMinterp['td'], SMinterp['superSMRbins'][2], 'g')
    plt.plot(SMinterp['td'], SMinterp['superSMRbins'][3], 'g')

    plt.plot(SMI['td'], SMI['smr00'], 'b', label = 'real')
    plt.plot(SMI['td'], SMI['smr06'], 'b')
    plt.plot(SMI['td'], SMI['smr12'], 'b')
    plt.plot(SMI['td'], SMI['smr18'], 'b')

    plt.xlim([SMinterp['td'][0], SMinterp['td'][-1]])
    plt.legend()
    plt.grid(True)
    plt.ylabel("SMR 6-hour bins", fontsize = 20)
    plt.show()

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def Z_Tensor_1D(resistivities, thicknesses, frequencies):
    """Calculate 1D Z-Tensor for given ground resistivity profile.
    Parameters
    -----------
    resistivities = array or list of resistivity values in Ohm.m
    thicknesses = array or list of thicknesses in m.
        **len(resistivities) must be len(thicknesses) + 1**
    frequencies = array or list of frequencies to get response of
    
    Returns
    -----------
    Z = complex array of Z tensor values
    
    Taken from:
    http://www.digitalearthlab.com/tutorial/tutorial-1d-mt-forward/"""
    
    if len(resistivities) != len(thicknesses) + 1:
        print("Length of inputs incorrect!")
        return 
    
    mu = 4*np.pi*1E-7; #Magnetic Permeability (H/m)
    n = len(resistivities);
    master_Z, master_absZ, master_phase = [], [], []

    for frequency in frequencies:   
        w =  2*np.pi*frequency;       
        impedances = list(range(n));
        #compute basement impedance
        impedances[n-1] = np.sqrt(w*mu*resistivities[n-1]*1j);
       
        for j in range(n-2,-1,-1):
            resistivity = resistivities[j];
            thickness = thicknesses[j];
      
            # 3. Compute apparent resistivity from top layer impedance
            #Step 2. Iterate from bottom layer to top(not the basement) 
            # Step 2.1 Calculate the intrinsic impedance of current layer
            dj = np.sqrt((w * mu * (1.0/resistivity))*1j);
            wj = dj * resistivity;
            # Step 2.2 Calculate Exponential factor from intrinsic impedance
            ej = np.exp(-2*thickness*dj);                     
        
            # Step 2.3 Calculate reflection coeficient using current layer
            #          intrinsic impedance and the below layer impedance
            belowImpedance = impedances[j + 1];
            rj = (wj - belowImpedance)/(wj + belowImpedance);
            re = rj*ej; 
            Zj = wj * ((1 - re)/(1 + re));
            impedances[j] = Zj;    
    
        # Step 3. Compute apparent resistivity from top layer impedance
        Z = impedances[0];
        phase = math.atan2(Z.imag, Z.real)
        master_Z.append(Z)
        master_absZ.append(abs(Z))
        master_phase.append(phase)
        #master_res.append((absZ * absZ)/(mu * w))
    return np.array(master_Z)

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def E_Field_1D(bx, by, resistivities, thicknesses, timestep = 60., Z = None, calc_Z = True, pad = True, padnum = 150):
    """Calculate horizontal E-field components given Bx, By, resistivities and thicknesses.
    
    Parameters
    -----------
    bx, by = array of Bx, By timeseries in nT
    resistivities = array or list of resistivity values in Ohm.m
    thicknesses = array or list of thicknesses in m.
        **len(resistivities) must be len(thicknesses) + 1**
    timestep = time between samples (default is 60. for minute sampling)
    
    Z = complex Z-tensor array. If not supplied, Z will be calculated from input
        resistivities and thicknesses
    
    Returns
    -----------
    ext, eyt = arrays of electric field components in mV/km"""
    
    if pad == False:
        new_bx = bx
        new_by = by
    else:
        new_bx = np.concatenate((bx[:padnum], bx, bx[-padnum:][::-1]))
        new_by = np.concatenate((by[:padnum], by, by[-padnum:][::-1]))
    
    mu0 = 4*np.pi * 1e-7
    freq = np.fft.fftfreq(new_bx.size, d = timestep)
    freq[0] = 1e-100

    if calc_Z == True:  # if you need to calculate Z
        Z = Z_Tensor_1D(resistivities, thicknesses, freq)
        
    bx_fft = np.fft.fft(new_bx)
    by_fft = np.fft.fft(new_by)

    exw = Z * by_fft/mu0; 
    eyw = -1 * Z * bx_fft/mu0

    ext = 1e-3 * np.fft.ifft(exw).real
    eyt = 1e-3 * np.fft.ifft(eyw).real

    if pad == False:
        return ext, eyt
    else:
        return ext[padnum:-padnum], eyt[padnum:-padnum]
        
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

def EField1DCalculation(BX, BY, TD):
    """Calculates Ex and Ey using a resistive 1-D ground model

    This function assumes that it is given 60 second data
    
    **WARNING** and time-series with nans will return all nans for that location
    I'm not sure what the best way to handle these are, probably a linear interpolation
    before calculating E-fields
    
    Also, be wary of FFT edge values at start and end of E-field calculation

    Parameters
    -----------
    BX, BY = arrays of Bx(north) and By(east) from SM functions (e.g., SM['dBn'])
    TD = timedate array

    Returns
    -----------
    EX, EY = arrays of Ex, Ey in mV/km, same shape as input Bx, By
    """
    Qres = np.array([500., 150., 20., 300., 100., 10., 1.])
    Qthick = 1000. * np.array([4., 6., 5., 65., 300., 200.])

    freq = np.fft.fftfreq(len(TD), d = 60.)
    freq[0] = 1e-100
    ZZ = Z_Tensor_1D(Qres, Qthick, freq)

    EX, EY = np.zeros((BX.shape)), np.zeros((BY.shape))

    for i, v in enumerate(EX.T):
        bx, by = BX[:,i], BY[:,i]
        
        bx[np.isnan(bx)] = 0
        bx[np.isnan(bx)] = 0
        
        ex, ey = E_Field_1D(bx, by, Qres, Qthick, timestep = 60., Z = ZZ, calc_Z = True, pad = True, padnum = 150)

        EX[:,i] = ex
        EY[:,i] = ey

        print("E-field calc: ", i)

    return EX, EY
    
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
"""
# EXAMPLE RUN
user = 'blakese' # username for Supermag downloads
savefolder = '/glade/u/home/sblake/indices/FINAL/DATA/' # folder where the Supermag jsons are stored/fetched

# Working Example
# Read in all needed data from .h5 simulation output
filename = '/glade/u/home/skareem/Work/gemmie/stpatex/stpatdb.deltab.h5'
#filename = '/glade/u/home/skareem/Work/dass/Data/newdassh/dassdb.deltab.h5'
print("Reading SIM data")
SIM = ReadSimData(filename)
start = SIM['td'][0] # start time of simulation data
numofdays = 3        # going to  

print("Fetching SM indices")
SMI  = FetchSMIndices(user, start, numofdays)

print("Fetching SM data")
SM = FetchSMData(user, start, numofdays, savefolder, badfrac = 0.1)

print("Interpolating SIM data") # interpolates and calculates SM indices
SMinterp = InterpolateSimData(SIM, SM)

print("Making Indices Plot")
MakeIndicesPlot(SMI, SMinterp, fignumber = 1)

print("Making Contour Plot")
MakeContourPlots(SM, SMinterp, maxx = 1000, fignumber = 2)

print("Calculating E-Field for SM data")
BX, BY, TD = SM['BNm'], SM['BEm'], SM['td']
EX, EY = EField1DCalculation(BX, BY, TD)
SM['Ex'] = EX
SM['Ey'] = EY

"""



