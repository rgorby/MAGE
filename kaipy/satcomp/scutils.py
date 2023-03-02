
# Import standard modules.
import datetime
import os
import subprocess
from xml.dom import minidom

# Import supplemental modules.
from astropy.coordinates import SkyCoord
import astropy.units as u
from cdasws import CdasWs
import h5py
import numpy as np
from spacepy.coordinates import Coords
import spacepy.datamodel as dm
from spacepy.time import Ticktock
from sunpy.coordinates import frames

# Import project-specific modules.
import kaipy.kaijson as kj
import kaipy.kaiTools as kaiTools
import kaipy.kdefs


# Module constants.

# A very small value.
TINY = 1.0e-8

# Compute the path to the directory containing this module.
package_directory = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))

# Compute the path to the magnetospheric spacecraft metadata file.
MAGNETOSPHERIC_SPACECRAFT_METADATA_PATH = os.path.join(
    PACKAGE_DIRECTORY, "sc_cdasws_strs.json"
)

# Radius of Sun in kilometers.
R_SUN_KILOMETERS = u.Quantity(1*u.Rsun, u.km).value

# Internal name to use to refer to the ephemeris positions.
SC_DATA_EPHEMERIS_NAME = "Ephemeris"

# Internal name to use to refer to the dates/times for the ephemeris positions.
SC_DATA_EPHEMERIS_EPOCH_NAME = "Ephemeris_Epoch"

# HTTP status codes.
HTTP_STATUS_OK = 200
HTTP_STATUS_NOT_FOUND = 404

# Minimum and maximum radius of heliospheric results grid, in units of
# solar radii.
R_MIN_HELIO = 21.5
R_MAX_HELIO = 220.0
scstrs_fname = os.path.join(package_directory, 'sc_cdasws_strs.json')

#======
#General
#======

def trilinterp(xbnd, ybnd, zbnd, valbnd, x, y, z):
    """3D linear interpolation
        xbnd,ybnd,zbnd: 2-element arrays each of the bounding dimensions
        valbnd: 2x2x2 array of variable to be interpolated
        x, y, z: point inside bounds
    """

    xd = (x - xbnd[0])/(xbnd[1]-xbnd[0])
    yd = (y - ybnd[0])/(ybnd[1]-ybnd[0])
    zd = (z - zbnd[0])/(zbnd[1]-zbnd[0])
    v00 = valbnd[0,0,0]*(1-xd) + valbnd[1,0,0]*xd
    v01 = valbnd[0,0,1]*(1-xd) + valbnd[1,0,1]*xd
    v10 = valbnd[0,1,0]*(1-xd) + valbnd[1,1,0]*xd
    v11 = valbnd[0,1,1]*(1-xd) + valbnd[1,1,1]*xd
    v0 = v00*(1-yd) + v10*yd
    v1 = v01*(1-yd) + v11*yd
    v = v0*(1-zd) + v1*zd

    return v

def varMap_1D(og, ng, var):
    """Map variable from one grid to another
    og: old grid
    ng: new grid
    var: variable to re-map
    """
    varnew = np.zeros((len(ng)))

    for e in range(len(ng)):
        if ng[e] < og[0] or ng[e] > og[-1]:
            continue

        idx = 0
        while og[idx+1] < ng[e]: idx += 1

        glow = og[idx]
        ghigh = og[idx+1]
        d = (ng[e] - glow)/(ghigh-glow)
        

        varnew[e] = var[idx]*(1-d) + var[idx+1]*d
    return varnew

def getWeights_ConsArea(og, og_lower, og_upper, ng, ng_lower, ng_upper):
    """Calculate overlap (weights) to map values on one grid to another,
        where total width in grid dimension are conserved
        (i.e. properly map RCM eetas to uniform grid)
        og, og_lower, og_upper: old grid and lower/upper bounds of each grid point
        ng, ng_lower, ng_upper: new grid and lower/upper bounds of each grid point
    """
    """Example:
        og: || |  |  |   |   |     |     |
        ng: |  |  |  |  |  |  |  |  |  |  |
        For each cell center on ng, calc which og cells overlap and the fraction of overlap
    """

    Nog = len(og)
    Nng = len(ng)
    weightMap = [[] for e in range(Nng) ]  # Ne x (nx2)
    for iNG in range(Nng):
        ng_l = ng_lower[iNG]
        ng_u = ng_upper[iNG]
        ng_w = ng_u - ng_l  # cell width
        frac_arr = []
        for k in range(Nog):
            #Do these two cells overlap
            if ng_l <= og_upper[k] and ng_u >= og_lower[k]:
                #Get overlap bounds and width
                ovl_lower = og_lower[k] if og_lower[k] > ng_l else ng_l
                ovl_upper = og_upper[k] if og_upper[k] < ng_u else ng_u
                ovl_width = ovl_upper - ovl_lower
                #og_width = og_upper[k] - og_lower[k]
                frac_arr.append([k, ovl_width/ng_w])
        weightMap[iNG] = frac_arr
    return weightMap

def computeErrors(obs,pred):
    MAE = 1./len(obs) *np.sum(np.abs(obs-pred))
    MSE = 1./len(obs) * np.sum((obs-pred)**2)
    RMSE = np.sqrt(MSE)
    MAPE = 1./len(obs) * np.sum(np.abs(obs-pred)/
        np.where(abs(obs) < TINY,TINY,abs(obs)))
    RSE = (np.sum((obs-pred)**2)/
        np.where(np.sum((obs-np.mean(obs))**2)<TINY,TINY,np.sum((obs-np.mean(obs))**2)))
    PE = 1-RSE
    return MAE,MSE,RMSE,MAPE,RSE,PE

#======
#Cdaweb-related
#======


def getScIds(spacecraft_data_file:str=scstrs_fname, doPrint:bool=False):
    """Fetch spacecraft descriptions from the database file.

    Load the spacecraft descriptions (a file) containing information needed
    to get spacecraft data from CDAWeb.

    Parameters
    ----------
    spacecraft_data_file: str, default "sc_cdasws_strs.json"
        Name of file containing spacecraft descriptions.
    doPrint: bool, default False
        If True, print a summary of the spacecraft descriptions.

    Returns
    -------
    scdict: dict
        Dictionary of spacecraft descriptions.
    """
    # Read the spacecraft database.
    scdict = kj.load(spacecraft_data_file)

    # Print a summary, if requested.
    if doPrint:
        print("Retrievable spacecraft data:")
        for sc in scdict:
            print('	 ' + sc)
            for v in scdict[sc]:
                print('	   ' + v)

    # Return the dictionary of spacraft descriptions.
    return scdict


def getCdasDsetInterval(dsName):
    cdas = CdasWs()

    data = cdas.get_datasets(idPattern=dsName)
    if len(data) == 0:
        return None, None
    tInt = data[0]['TimeInterval']
    return tInt['Start'], tInt['End']

def pullVar(cdaObsId:str, cdaDataId, t0:str, t1:str, deltaT:float=60,
            epochStr:str="Epoch", doVerbose:bool=False):
    """Pull spacecraft data from CDAWeb.

    Pull spacecraft data from CDAWeb.

    Parameters
    ----------
    cdaObsId: str
        Dataset name.
    cdaDataId: str or list of str
        Desired variable(s) from dataset.
    t0: str
        Data start time, formatted as '%Y-%m-%dT%H:%M:%S.%f'.
    t1: str
        Data end time, formatted as '%Y-%m-%dT%H:%M:%S.%f'.
    deltaT: float, default 60
        Time cadence (seconds), used when interpolating through time with
        no data.
    epochStr: str, default "Epoch"
        Name of time variable in dataset.
    doVerbose: bool, default False
        Helpful for debugging/diagnostics.

    Returns
    -------
    status: dict
        Status information returned for the query.
    data: spacepy.pycdf.CDFCopy
        Object containing data returned by the query, None if no results.
    """

    # Specify how CDAWeb should bin the data.
    binData = {
        'interval': deltaT, 
        'interpolateMissingValues': True,
        'sigmaMultipler': 4
    }

    # Create the CDAWeb query object.
    cdas = CdasWs()

    # Perform the query.
    status, data = cdas.get_data(cdaObsId, cdaDataId, t0, t1, binData=binData)

    # Process the query status.
    if status["http"]["status_code"] in (204, 404):
        # 204 = No Content
        # 404 = Not Found
        if doVerbose:
            print("No data found.")
    elif status['http']['status_code'] != 200 or data is None:
        # Handle the case where CdasWs just doesn't work if you give it variables in arg 2
        # If given empty var list instead, it'll return the full day on day in t0, and that's it
        # So, call for as many days as we need data for and build one big data object
        if doVerbose: print("Bad pull, trying to build day-by-day")

        if '.' in t0:
            t0dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%S.%fZ")
            t1dt = datetime.datetime.strptime(t1, "%Y-%m-%dT%H:%M:%S.%fZ")
        else:
            t0dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%SZ")
            t1dt = datetime.datetime.strptime(t1, "%Y-%m-%dT%H:%M:%SZ")
        numDays = t1dt.day-t0dt.day + 1 #Number of days we want data from
        if doVerbose: print("numDays: " + str(numDays))

        tstamp_arr = []
        tstamp_deltas = []
        for i in range(numDays):
            tstamp_arr.append((t0dt + datetime.timedelta(days=i)).strftime("%Y-%m-%dT%H:%M:%SZ"))
            tstamp_deltas.append((t0dt + datetime.timedelta(days=i+1)).strftime("%Y-%m-%dT%H:%M:%SZ"))
        if doVerbose: print("Tstamp_arr: " + str(tstamp_arr))
        #Get first day
        status, data = cdas.get_data(cdaObsId, [], tstamp_arr[0], tstamp_deltas[0], binData=binData)
        if doVerbose: print("Pulling " + t0)
        
        if status['http']['status_code'] != 200:
            # If it still fails, its some other problem and we'll die
            if doVerbose: print("Still bad pull. Dying.")
            return status,data
        if data is None:
            if doVerbose: print("Cdas responded with 200 but returned no data")
            return status,data
        if epochStr not in data.keys():
            if doVerbose: print(epochStr + " not in dataset, can't build day-by-day")
            data = None
            return status,data
        if isinstance(cdaDataId,str):
            if doVerbose: print(cdaDataId + " not in dataset, can't build day-by-day")
            # Mimic the cdasws return code for case when id isn't provided
            if not (cdaDataId in data.keys()):
                status = {'http': {'status_code': 404}}
                data = None
                return status,data
        if isinstance(cdaDataId,list):
            for item in cdaDataId:
                print(item + " not in dataset, can't build day-by-day")
                if not (item in data.keys()):
                    # Mimic the cdasws return code for case when id isn't provided
                    status = {'http': {'status_code': 404}}
                    data = None
                    return status,data

        #Figure out which axes are the epoch axis in each dataset so we can concatenate along it
        dk = list(data.keys())
        nTime = len(data[epochStr])
        cataxis = np.array([-1 for i in range(len(dk))])
        for k in range(len(dk)):
            shape = np.array(data[dk[k]].shape)
            for i in range(len(shape)):
                if shape[i] == nTime:
                    cataxis[k] = i
                    continue

        #Then append rest of data accordingly
        for i in range(1,numDays):
            if doVerbose: print("Pulling " + str(tstamp_arr[i]))
            status, newdata = cdas.get_data(cdaObsId, [], tstamp_arr[i], tstamp_deltas[i], binData=binData)
            for k in range(len(dk)):
                if cataxis[k] == -1:
                    continue
                else:
                    key = dk[k]
                    data[key] = np.concatenate((data[key], newdata[key]), axis=cataxis[k])
    else:
        if doVerbose:
            print("Got data in one pull.")

    # Return the query status and results.
    return status, data

def addVar(mydata,scDic,varname,t0,t1,deltaT,epochStr='Epoch'):
    #print(scDic,varname,idname,dataname,scDic[idname])
    if scDic[varname]['Id'] is not None:
        status,data = pullVar(scDic[varname]['Id'],scDic[varname]['Data'],t0,t1,deltaT,epochStr=epochStr)
        #print(status)
        if status['http']['status_code'] == 200 and data is not None:
            mydata[varname] = dm.dmarray(data[scDic[varname]['Data']],
                                         attrs=data[scDic[varname]['Data']].attrs)
            #mydata.tree(attrs=True)
    else:
        #Mimic the cdasws return code for case when id isn't provided
        status = {'http': {'status_code': 404}}
    return status


def getSatData(scDic:dict, t0:str, t1:str, deltaT:float):
    """Fetch spacecraft data in the specified time range.

    Fetch spacecraft data for the specified time range, at the specified cadence.

    Parameters
    ----------
    scDic: dict
        Spacecraft descriptive information.
    t0: str
        Start time for data, in format "%Y-%m-%dT%H:%M:%SZ".
    t1: str
        Stop time for data, in format "%Y-%m-%dT%H:%M:%SZ".
    deltaT: float
        Cadence for requested spacecraft data (seconds).

    Returns
    -------
    status: dict
        Query status returned by CDAWeb.
    mydata: spacepy.datamodel.SpaceData
        All of the spacecraft for the specified time range and cadence.
    """
    # Fetch the empheris data. If not found, abort this query.
    status, data = pullVar(
        scDic['Ephem']['Id'], scDic['Ephem']['Data'], t0, t1, deltaT
    )
    if status['http']['status_code'] != 200 or data is None:
        print('Unable to get data for ', scDic['Ephem']['Id'])
        return status, data

    # Create a new SpaceData object to hold the results of the query.
    mydata = dm.SpaceData(attrs={'Satellite':data.attrs['Source_name']})

    # Determine which of the returned data should be used for the time of the
    # ephemeris positions.
    if 'Epoch_bin' in data.keys():
        mytime = data['Epoch_bin']
        epochStr = 'Epoch_bin'
    elif 'Epoch' in data.keys():
        mytime = data['Epoch']
        epochStr = 'Epoch'
    elif ([key for key in data.keys() if key.endswith('_state_epoch')]):
        epochStr = [key for key in data.keys() if key.endswith('_state_epoch')][0]
        #mytime = data[[key for key in data.keys()
        #if key.endswith('_state_epoch')][0]]
        mytime = data[epochStr]
    else:
        print('Unable to determine time type')
        status = {'http':{'status_code':404}}
        return status, data

    # Extract the times assigned to the ephemeris positions.
    mydata['Epoch_bin'] = dm.dmarray(mytime, attrs=mytime.attrs)

    # Extract the ephemeris positions.
    mydata['Ephemeris'] = dm.dmarray(
        data[scDic['Ephem']['Data']],
        attrs=data[scDic['Ephem']['Data']].attrs
    )

    # Now fetch the data measured by the spacecraft.
    keys = ['MagneticField', 'Velocity', 'Density', 'Pressure', "Speed", "Temperature"]
    for key in keys:
        if key in scDic:
            status1 = addVar(
                mydata, scDic, key, t0, t1, deltaT, epochStr=epochStr
            )

    #Add any metavar since they might be needed for unit/label determination
    search_key = 'metavar'
    res = [key for key,val in data.items() if search_key in key]
    for name in res:
        try:
            len(mydata[name])
        except:
            mydata[name] = dm.dmarray([data[name]],attrs=data[name].attrs)
        else:
            mydata[name] = dm.dmarray(data[name],attrs=data[name].attrs)

    return status,mydata


#======
#Shared data derivations
#======

def xyz_to_L(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = np.arctan2(z, np.sqrt(x**2 + y**2))
    #Convert sc location to L shell, assuming perfect dipole
    lat = lat*np.pi/180.0  # deg to rad
    return r/np.cos(lat)**2

def getJScl(Bmag,Beq,en=2.0):
    #Given sin^n(alpha) dep. on intensity calculate fraction based on accessible Alpha

    Na = 360
    A = np.linspace(0,0.5*np.pi,Na)
    da = A[1]-A[0]
    Ia = np.sin(A)**en
    Ic = np.zeros(Ia.shape)
    Nt = len(Bmag)
    I0 = Ia.sum()

    It = np.zeros(Nt)
    for n in range(Nt):
        if (Bmag[n]<TINY):
            It[n] = 0.0
        else:
            Ac = np.arcsin(np.sqrt(Beq[n]/Bmag[n]))
            Ic[:] = Ia[:]
            Icut = (A>Ac)
            Ic[Icut] = 0.0
            It[n] = Ic.sum()/I0
    return It

def genSCXML(fdir,ftag,
    scid="sctrack_A",h5traj="sctrack_A.h5",numSegments=1):

    (fname,isMPI,Ri,Rj,Rk) = kaiTools.getRunInfo(fdir,ftag)
    root = minidom.Document()
    xml = root.createElement('Kaiju')
    root.appendChild(xml)
    chimpChild = root.createElement('Chimp')
    scChild = root.createElement("sim")
    scChild.setAttribute("runid",scid)
    chimpChild.appendChild(scChild)
    fieldsChild = root.createElement("fields")
    fieldsChild.setAttribute("doMHD","T")
    fieldsChild.setAttribute("grType","LFM")
    fieldsChild.setAttribute("ebfile",ftag)
    if isMPI:
        fieldsChild.setAttribute("isMPI","T")
    chimpChild.appendChild(fieldsChild)
    if isMPI:
        parallelChild = root.createElement("parallel")
        parallelChild.setAttribute("Ri","%d"%Ri)
        parallelChild.setAttribute("Rj","%d"%Rj)
        parallelChild.setAttribute("Rk","%d"%Rk)
        chimpChild.appendChild(parallelChild)
    unitsChild = root.createElement("units")
    unitsChild.setAttribute("uid","EARTH")
    chimpChild.appendChild(unitsChild)
    trajChild = root.createElement("trajectory")
    trajChild.setAttribute("H5Traj",h5traj)
    trajChild.setAttribute("doSmooth","F")
    chimpChild.appendChild(trajChild)
    # if doTrc:
    #     outChild = root.createElement('output')
    #     outChild.setAttribute('doTrc', "T")
    #     chimpChild.appendChild(outChild)
    if numSegments > 1:
        parInTimeChild = root.createElement("parintime")
        parInTimeChild.setAttribute("NumB","%d"%numSegments)
        chimpChild.appendChild(parInTimeChild)
    xml.appendChild(chimpChild)
    return root


#======
#SCTrack
#======

def convertGameraVec(x,y,z,ut,fromSys,fromType,toSys,toType):
    invec = Coords(np.column_stack((x,y,z)),fromSys,fromType, use_irbem=False)
    invec.ticks = Ticktock(ut)
    outvec = invec.convert(toSys,toType)
    return outvec

def createInputFiles(data,scDic,scId,mjd0,sec0,fdir,ftag,numSegments):
    Re = 6380.0
    toRe = 1.0
    if 'UNITS' in data['Ephemeris'].attrs:
        if "km" in data['Ephemeris'].attrs['UNITS']:
            toRe = 1.0/Re
    elif 'UNIT_PTR' in data['Ephemeris'].attrs:
        if data[data['Ephemeris'].attrs['UNIT_PTR']][0]:
            toRe = 1.0/Re
    if 'SM' == scDic['Ephem']['CoordSys']:
        smpos = Coords(data['Ephemeris'][:,0:3]*toRe,'SM','car', use_irbem=False)
        smpos.ticks = Ticktock(data['Epoch_bin'])
    elif 'GSM' == scDic['Ephem']['CoordSys']:
        scpos = Coords(data['Ephemeris'][:,0:3]*toRe,'GSM','car', use_irbem=False)
        scpos.ticks = Ticktock(data['Epoch_bin'])
        smpos = scpos.convert('GSE','car')
        scpos = Coords(data['Ephemeris'][:,0:3]*toRe,'GSM','car', use_irbem=False)
        scpos.ticks = Ticktock(data['Epoch_bin'])
        smpos = scpos.convert('SM','car')
    elif 'GSE'== scDic['Ephem']['CoordSys']:
        scpos = Coords(data['Ephemeris'][:,0:3]*toRe,'GSE','car', use_irbem=False)
        scpos.ticks = Ticktock(data['Epoch_bin'])
        smpos = scpos.convert('SM','car')
    else:
        print('Coordinate system transformation failed')
        return
    elapsedSecs = (smpos.ticks.getMJD()-mjd0)*86400.0+sec0
    scTrackName = os.path.join(fdir,scId+".sc.h5")
    with h5py.File(scTrackName,'w') as hf:
        hf.create_dataset("T" ,data=elapsedSecs)
        hf.create_dataset("X" ,data=smpos.x)
        hf.create_dataset("Y" ,data=smpos.y)
        hf.create_dataset("Z" ,data=smpos.z)
    chimpxml = genSCXML(fdir,ftag,
        scid=scId,h5traj=os.path.basename(scTrackName),numSegments=numSegments)
    xmlFileName = os.path.join(fdir,scId+'.xml')
    with open(xmlFileName,"w") as f:
        f.write(chimpxml.toprettyxml())

    return (scTrackName,xmlFileName)


def addGAMERA(data,scDic,h5name):
    h5file = h5py.File(h5name, 'r')
    ut = kaiTools.MJD2UT(h5file['MJDs'][:])

    bx = h5file['Bx']
    by = h5file['By']
    bz = h5file['Bz']

    if not 'MagneticField' in scDic:
        toCoordSys = 'GSM'
    else:
        toCoordSys = scDic['MagneticField']['CoordSys']
    lfmb_out = convertGameraVec(bx[:],by[:],bz[:],ut,
        'SM','car',toCoordSys,'car')
    data['GAMERA_MagneticField'] = dm.dmarray(lfmb_out.data,
        attrs={'UNITS':bx.attrs['Units'],
        'CATDESC':'Magnetic Field, cartesian'+toCoordSys,
        'FIELDNAM':"Magnetic field",'AXISLABEL':'B'})
    vx = h5file['Vx']
    vy = h5file['Vy']
    vz = h5file['Vz']
    if not 'Velocity' in scDic:
        toCoordSys = 'GSM'
    else:
        toCoordSys = scDic['Velocity']['CoordSys']
    lfmv_out = convertGameraVec(vx[:],vy[:],vz[:],ut,
        'SM','car',toCoordSys,'car')
    data['GAMERA_Velocity'] = dm.dmarray(lfmv_out.data,
        attrs={'UNITS':vx.attrs['Units'],
        'CATDESC':'Velocity, cartesian'+toCoordSys,
               'FIELDNAM':"Velocity",'AXISLABEL':'V'})
    speed = h5file["Vx"]
    data['GAMERA_Speed'] = dm.dmarray(-speed[:],
        attrs={'UNITS':speed.attrs['Units'],
        'CATDESC':'Speed','FIELDNAM':"Speed",'AXISLABEL':'vr'})
    den = h5file['D']
    data['GAMERA_Density'] = dm.dmarray(den[:],
        attrs={'UNITS':den.attrs['Units'],
        'CATDESC':'Density','FIELDNAM':"Density",'AXISLABEL':'n'})
    pres = h5file['P']
    data['GAMERA_Pressure'] = dm.dmarray(pres[:],
        attrs={'UNITS':pres.attrs['Units'],
        'CATDESC':'Pressure','FIELDNAM':"Pressure",'AXISLABEL':'P'})
    inDom = h5file['inDom']
    data['GAMERA_inDom'] = dm.dmarray(inDom[:],
        attrs={'UNITS':inDom.attrs['Units'],
        'CATDESC':'In GAMERA Domain','FIELDNAM':"InDom",
        'AXISLABEL':'In Domain'})
    return


def matchUnits(data):
    vars = ['Density','Pressure','Temperature','Velocity','MagneticField']
    for var in vars:
        try:
            data[var]
        except:
            print(var,'not in data')
        else:
            if (data[var].attrs['UNITS'] == data['GAMERA_'+var].attrs['UNITS'].decode()):
                print(var,'units match')
            else:
                if 'Density' == var:
                    if (data[var].attrs['UNITS'] == 'cm^-3' or data[var].attrs['UNITS'] == '/cc'):
                        data[var].attrs['UNITS'] = data['GAMERA_'+var].attrs['UNITS']
                        print(var,'units match')
                    else:
                        print('WARNING ',var,'units do not match')
                if 'Velocity' == var:
                    if (data[var].attrs['UNITS'] == 'km/sec'):
                        data[var].attrs['UNITS'] = data['GAMERA_'+var].attrs['UNITS']
                        print(var,'units match')
                    else:
                        print('WARNING ',var,'units do not match')
                if 'MagneticField' == var:
                    if (data[var].attrs['UNITS'] == '0.1nT'):
                        print('Magnetic Field converted from 0.1nT to nT')
                        data[var]=data[var]/10.0
                        data[var].attrs['UNITS'] = 'nT'
                    else:
                        print('WARNING ',var,'units do not match')
                if 'Pressure' == var:
                    print('WARNING ',var,'units do not match')
                if 'Temperature' == var:
                    print('WARNING ',var,'units do not match')

    return

def extractGAMERA(data,scDic,scId,mjd0,sec0,fdir,ftag,cmd,numSegments,keep):

    (scTrackName,xmlFileName) = createInputFiles(data,scDic,scId,
        mjd0,sec0,fdir,ftag,numSegments)

    if 1 == numSegments:
        sctrack = subprocess.run([cmd, xmlFileName], cwd=fdir, capture_output=True,
                            # stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True)
        with open(os.path.join(fdir, "sctrack.out"), "w") as f:
            f.write(sctrack.stdout)

        #print(sctrack)
        h5name = os.path.join(fdir, scId + '.sc.h5')

    else:
        process = []
        for seg in range(1,numSegments+1):
            process.append(subprocess.Popen([cmd, xmlFileName,str(seg)],
                            cwd=fdir,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                            text=True))
        for proc in process:
            proc.communicate()
        h5name = mergeFiles(scId,fdir,numSegments)


    addGAMERA(data,scDic,h5name)

    if not keep:
        subprocess.run(['rm',h5name])
        subprocess.run(['rm',xmlFileName])
        subprocess.run(['rm',scTrackName])
        if numSegments > 1:
            h5parts = os.path.join(fdir,scId+'.*.sc.h5')
            subprocess.run(['rm',h5parts])
    return

def copy_attributes(in_object, out_object):
    '''Copy attributes between 2 HDF5 objects.'''
    for key, value in list(in_object.attrs.items()):
        out_object.attrs[key] = value


def createMergeFile(fIn,fOut):
    iH5 = h5py.File(fIn,'r')
    oH5 = h5py.File(fOut,'w')
    copy_attributes(iH5,oH5)
    for Q in iH5.keys():
        oH5.create_dataset(Q,data=iH5[Q],maxshape=(None,))
        copy_attributes(iH5[Q],oH5[Q])
    iH5.close()
    return oH5

def addFileToMerge(mergeH5,nextH5):
    nS = nextH5.attrs['nS']
    nE = nextH5.attrs['nE']
    for varname in mergeH5.keys():
        dset = mergeH5[varname]
        dset.resize(dset.shape[0]+nextH5[varname].shape[0],axis=0)
        dset[nS-1:nE]=nextH5[varname][:]
    return

def mergeFiles(scId,fdir,numSegments):
    seg = 1
    inH5Name = os.path.join(fdir,scId+'.%04d'%seg+'.sc.h5')
    mergeH5Name = os.path.join(fdir,scId+'.sc.h5')
    mergeH5 = createMergeFile(inH5Name,mergeH5Name)
    #print(inH5Name,mergeH5Name)
    for seg in range(2,numSegments+1):
        nextH5Name = os.path.join(fdir,scId+'.%04d'%seg+'.sc.h5')
        nextH5 = h5py.File(nextH5Name,'r')
        addFileToMerge(mergeH5,nextH5)

    return mergeH5Name

def genSatCompPbsScript(scId,fdir,cmd,account='P28100045'):
    headerString = """#!/bin/tcsh
#PBS -A %s
#PBS -N %s
#PBS -j oe
#PBS -q casper
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1
"""
    moduleString = """module purge
module load git/2.22.0 intel/18.0.5 hdf5/1.10.5 impi/2018.4.274
module load ncarenv/1.3 ncarcompilers/0.5.0 python/3.7.9 cmake/3.14.4
module load ffmpeg/4.1.3 paraview/5.8.1 mkl/2018.0.5
ncar_pylib /glade/p/hao/msphere/gamshare/casper_satcomp_pylib
module list
"""
    commandString = """cd %s
setenv JNUM ${PBS_ARRAY_INDEX}
date
echo 'Running analysis'
%s %s $JNUM
date
"""
    xmlFileName = os.path.join(fdir,scId+'.xml')
    pbsFileName = os.path.join(fdir,scId+'.pbs')
    pbsFile = open(pbsFileName,'w')
    pbsFile.write(headerString%(account,scId))
    pbsFile.write(moduleString)
    pbsFile.write(commandString%(fdir,cmd,xmlFileName))
    pbsFile.close()

    return pbsFileName

def genSatCompLockScript(scId,fdir,account='P28100045'):
    headerString = """#!/bin/tcsh
#PBS -A %s
#PBS -N %s
#PBS -j oe
#PBS -q casper
#PBS -l walltime=0:15:00
#PBS -l select=1:ncpus=1
"""
    commandString = """cd %s
touch %s
"""
    pbsFileName = os.path.join(fdir,scId+'.done.pbs')
    pbsFile = open(pbsFileName,'w')
    pbsFile.write(headerString%(account,scId))
    pbsFile.write(commandString%(fdir,scId+'.lock'))
    pbsFile.close()

    return pbsFileName

def errorReport(errorName,scId,data):

    keysToCompute = []
    keys = data.keys()

    print('Writing Error to ',errorName)
    f = open(errorName,'w')
    if 'Density' in keys:
        keysToCompute.append('Density')
    if 'Pressue' in keys:
        keysToCompute.append('Pressue')
    if 'Temperature' in keys:
        keysToCompute.append('Temperature')
    if 'MagneticField' in keys:
        keysToCompute.append('MagneticField')
    if 'Velocity' in keys:
        keysToCompute.append('Velocity')

    for key in keysToCompute:
        #print('Plotting',key)
        if 'MagneticField' == key or 'Velocity' == key:
            for vecComp in range(3):
                maskedData = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
                    data[key][:,vecComp])
                maskedGamera = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
                    data['GAMERA_'+key][:,vecComp])
                MAE,MSE,RMSE,MAPE,RSE,PE = computeErrors(maskedData,maskedGamera)
                f.write(f'Errors for: {key},{vecComp}\n')
                f.write(f'MAE: {MAE}\n')
                f.write(f'MSE: {MSE}\n')
                f.write(f'RMSE: {RMSE}\n')
                f.write(f'MAPE: {MAPE}\n')
                f.write(f'RSE: {RSE}\n')
                f.write(f'PE: {PE}\n')
        else:
                maskedData = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
                    data[key][:])
                maskedGamera = np.ma.masked_where(data['GAMERA_inDom'][:]==0.0,
                    data['GAMERA_'+key][:])
                MAE,MSE,RMSE,MAPE,RSE,PE = computeErrors(maskedData,maskedGamera)
                f.write(f'Errors for: {key}\n')
                f.write(f'MAE: {MAE}\n')
                f.write(f'MSE: {MSE}\n')
                f.write(f'RMSE: {RMSE}\n')
                f.write(f'MAPE: {MAPE}\n')
                f.write(f'RSE: {RSE}\n')
                f.write(f'PE: {PE}\n')
    f.close()
    return


#-----------------------------------------------------------------------------

# HELIO ONLY BELOW THIS POINT.

# ADD DOCSTRINGS FOR RAISED EXCEPTIONS.

# Methods specific for comparing gamhelio output to data from heliospheric
# spacecraft.


def read_MJDc(path):
    """Read the MJD of the center of the WSA input map.

    Read the MJD of the center of the WSA input map. This value was originally
    read from the FITS file used to provide the WSA initial conditions for a
    gamhelio run. It is now stored as a global attribute in all of the HDF5-
    format files generated by gamhelio. This MJD value should be used as
    the date for constructing the gamhelio coordinate frame, needed to convert
    values from CDAWeb to gamhelio coordinates. 

    Parameters
    ----------
    path : str
        Path to HDF5 results file from gamhelio.

    Returns
    -------
    mjdc : float
        Value of MJDc global attribute.
    """
    with h5py.File(path, "r") as hf:
        mjdc = hf.attrs["MJDc"]
    return mjdc


def helio_pullVar(
        cdaObsId, cdaDataId,
        t0, t1, deltaT=60,
        epochStr="Epoch", doVerbose=False
):
    """Pull spacecraft data from CDAWeb.

    Pull spacecraft data from CDAWeb.

    Parameters
    ----------
    cdaObsId: str
        Dataset name.
    cdaDataId: str or list of str
        Desired variable(s) from dataset.
    t0: str
        Data start time, formatted as '%Y-%m-%dT%H:%M:%S.%f'.
    t1: str
        Data end time, formatted as '%Y-%m-%dT%H:%M:%S.%f'.
    deltaT: float, default 60
        Time cadence (seconds), used when interpolating through time with
        no data.
    epochStr: str, default "Epoch"
        Name of time variable in dataset.
    doVerbose: bool, default False
        Helpful for debugging/diagnostics.

    Returns
    -------
    status: dict
        Status information returned for the query.
    data: spacepy.pycdf.CDFCopy
        Object containing data returned by the query, None if no results.
    """

    # Specify how CDAWeb should bin the data.
    binData = {
        "interval": deltaT, 
        "interpolateMissingValues": True,
        "sigmaMultipler": 4
    }

    # Create the CDAWeb query object.
    cdas = CdasWs()

    # Perform the query.
    status, data = cdas.get_data(cdaObsId, cdaDataId, t0, t1, binData=binData)

    # Process the query status.
    if status["http"]["status_code"] in (204, HTTP_STATUS_NOT_FOUND):
        # 204 = No Content
        # HTTP_STATUS_NOT_FOUND = Not Found
        if doVerbose:
            print("No data found.")

    # Return the query status and results.
    return status, data


def helio_addVar(my_data, scDic, varname, t0, t1, deltaT, epochStr="Epoch"):
    cdaweb_dataset_name = scDic[varname]["Id"]
    cdaweb_variable_name = scDic[varname]["Data"]
    if cdaweb_variable_name is not None:
        status,data = helio_pullVar(
            cdaweb_variable_name, cdaweb_dataset_name,
            t0, t1, deltaT, epochStr=epochStr)
        if status["http"]["status_code"] == HTTP_STATUS_OK and data is not None:
            my_data[varname] = dm.dmarray(data[scDic[varname]["Data"]],
                                         attrs=data[scDic[varname]["Data"]].attrs)
    else:
        #Mimic the cdasws return code for case when id isn't provided
        status = {"http": {"status_code": HTTP_STATUS_NOT_FOUND}}
    return status


def get_helio_cdaweb_data(
    sc_id, sc_metadata, start_time, end_time, cdaweb_data_interval,
    verbose=False, debug=False
):
    """Fetch heliosphere spacecraft data in the specified time range.

    Fetch heliosphere spacecraft data for the specified time range, at the
    specified cadence. This function copies the raw CDAWeb data from the
    query result object into a local data object. The raw data must then
    be ingested to convert it to a form that gamhelio understands.

    Parameters
    ----------
    sc_id : str
        ID string for spacecraft.
    sc_metadata : dict
        Spacecraft descriptive information for the heliosphere spacecraft YML
        file.
    start_time : str
        Start time for data, in format "%Y-%m-%dT%H:%M:%SZ".
    end_time : str
        Stop time for data, in format "%Y-%m-%dT%H:%M:%SZ".
    cdaweb_data_interval : float
        Cadence for requested spacecraft data (seconds).
    verbose : bool, default False
        Set to True for printing verbose progress messages.
    debug : bool, default False
        Set to True for printing debugging messages.

    Returns
    -------
    my_data : spacepy.datamodel.SpaceData
        All of the spacecraft for the specified time range and cadence.

    Raises
    ------

    """
    # Fetch the ephemeris variable(s) from CDAWeb.
    cdaweb_dataset_name = sc_metadata["Ephem"]["Id"]
    cdaweb_variable_name = sc_metadata["Ephem"]["Data"]
    print("Fetching spacecraft %s, dataset %s, variable(s) %s from CDAWeb." %
          (sc_id, cdaweb_dataset_name, cdaweb_variable_name))
    cdaweb_query_status, cdaweb_query_results = helio_pullVar(
        cdaweb_dataset_name, cdaweb_variable_name,
        start_time, end_time, cdaweb_data_interval
    )
    if debug:
        print("cdaweb_query_status = %s" % cdaweb_query_status)
        print("cdaweb_query_results = %s" % cdaweb_query_results)

    # Return if no data found (which can happen if the requested time range
    # if before the spacecraft launch date).
    if cdaweb_query_results is None:
        return None

    # Create a new SpaceData object to hold the ingested results of the
    # query.
    # NOTE: This object is treated as a Python dictionary for adding and
    # removing data. However, it also has the "attrs" attribute which
    # is a dictionary holding attributes used by SpacePy.
    sc_data = dm.SpaceData(
        attrs={
            "Satellite": sc_id
        }
    )
    if debug:
        print("data = %s" % sc_data)
        print("data.attrs = %s" % sc_data.attrs)

    # Add any CDAWeb metavar since they might be needed for unit/label
    # determination. These metavars are added as new dm.dmarray objects,
    # and so must be copied as lists, even if they are scalars.
    # THIS IS MESSY. REWRITE IT.
    search_key = "metavar"
    for key in cdaweb_query_results:
        if not search_key in key:
            continue
        if verbose:
            print("Copying CDAWeb metavariable %s." % key)
        try:
            len(sc_data[key])
        except:
            sc_data[key] = dm.dmarray(
                [cdaweb_query_results[key]],
                attrs=cdaweb_query_results[key].attrs
            )
        else:
            sc_data[key] = dm.dmarray(
                cdaweb_query_results[key],
                attrs=cdaweb_query_results[key].attrs
            )

    # Extract the times assigned to the ephemeris times.
    sc_data[sc_metadata["Ephem"]["Epoch_Name"]] = dm.dmarray(
        cdaweb_query_results[sc_metadata["Ephem"]["Epoch_Name"]],
        attrs=cdaweb_query_results[sc_metadata["Ephem"]["Epoch_Name"]].attrs
    )

    # Extract the variable(s) which define the ephemeris positions.
    if isinstance(cdaweb_variable_name, list):
        for variable_name in cdaweb_variable_name:
            sc_data[variable_name] = dm.dmarray(
                cdaweb_query_results[variable_name],
                attrs=cdaweb_query_results[variable_name].attrs
            )
    else:
        sc_data[cdaweb_variable_name] = dm.dmarray(
            cdaweb_query_results[cdaweb_variable_name],
            attrs=cdaweb_query_results[cdaweb_variable_name].attrs
        )

    # Now fetch the physical data measured by the spacecraft.
    # Each spacecraft entry in the spacecraft metadata database can have
    # entries for any number of physical variables. For this comparison,
    # we are only interested in the following keys:
    # 1. Speed = radial velocity of the solar wind (km/s)
    # 2. MagneticField = components of the magnetic field (nT)
    # 3. Density = number density of solar wind (#/cc)
    # 4. Temperature = temperature of solar wind (K)
    # Note that in this case, the addVar() function directly adds each
    # new variable directly to the my_data object.
    variable_names = ["Speed", "MagneticField", "Density", "Temperature"]
    for variable_name in variable_names:
        cdaweb_dataset_name = sc_metadata[variable_name]["Id"]
        cdaweb_variable_name = sc_metadata[variable_name]["Data"]
        if verbose:
            print("Fetching spacecraft %s, dataset %s, variable %s from "
                  "CDAWeb." % (sc_id, cdaweb_dataset_name,
                               cdaweb_variable_name))
        cdaweb_query_status, cdaweb_query_results = helio_pullVar(
            cdaweb_dataset_name, cdaweb_variable_name,
            start_time, end_time, cdaweb_data_interval
        )
        if (cdaweb_query_status["http"]["status_code"] != HTTP_STATUS_OK or
            cdaweb_query_results is None) :
            print("No data found for spacecraft %s, dataset %s, variable %s!"
                  % (sc_id, cdaweb_dataset_name, cdaweb_variable_name))
            continue

        # Extract the CDAWeb variable(s) which define this local variable.
        # Then mask out the measured values that were interpolated by CDAWeb.
        # Such points are identified by entries where the VARNAME_NBIN
        # value is 0.
        if isinstance(cdaweb_variable_name, list):
            for variable_name in cdaweb_variable_name:
                sc_data[variable_name] = dm.dmarray(
                    cdaweb_query_results[variable_name],
                    attrs=cdaweb_query_results[variable_name].attrs
                )
                nbin_varname = variable_name.upper() + "_NBIN"
                w = np.where(cdaweb_query_results[nbin_varname] == 0)
                sc_data[variable_name][w] = np.nan
        else:
            sc_data[cdaweb_variable_name] = dm.dmarray(
                cdaweb_query_results[cdaweb_variable_name],
                attrs=cdaweb_query_results[cdaweb_variable_name].attrs
            )
            nbin_varname = cdaweb_variable_name.upper() + "_NBIN"
            w = np.where(cdaweb_query_results[nbin_varname] == 0)
            sc_data[cdaweb_variable_name][w] = np.nan


    # Return the accumulated data from CDAWeb.
    return sc_data


def ingest_cdaweb_ephemeris(sc_data, sc_metadata, MJDc, verbose=False, debug=False):
    """Convert CDAWeb spacecraft ephemeris to gamhelio format.

    Convert CDAWeb spacecraft ephemeris to gamhelio format. This process
    converts the CDAWeb coordinates (which may be in any supported coordinate
    system) to the Cartesian GH(MJDc) frame used by gamhelio. This function
    also makes a copy of the times associated with the ephemeris, and saves
    them with a standard name.

    Once the variables are ingested, they are referred to using the names of
    the local variables listed in the spacecraft metadata.

    Parameters
    ----------
    sc_data : dm.SpaceData
        SpaceData object for all data as originally returned from CDAWeb.
    sc_metadata : dict
        Spacecraft descriptive information for the heliosphere spacecraft YML
        file.
    MJDc : float
        Value of MJDc attribute from gamhelio result files.
    verbose : bool, default False
        Set to True for printing verbose progress messages.
    debug : bool, default False
        Set to True for printing debugging messages.

    Returns
    -------
    None
    """
    # Fetch the coordinate system for the ephemeris.
    cdaweb_coordinate_system = sc_metadata["Ephem"]['CoordSys']

    # Ingest the trajectory by converting it to the gamhelio  frame.
    if cdaweb_coordinate_system == "GSE":

        # GSE(t) coordinates from CDAWeb are provided as a single 2-D array
        # of Cartesian (x, y, z) values, of shape (n, 3), called "XYZ_GSE".
        # These values are in units of kilometers. The ephemeris time is
        # in a variable called "Epoch_bin".

        # Convert the GSE(t) ephemeris locations to the gamhelio frame
        # (GH(MJDc) so that sctrack.x can interpolate gamhelio model data to
        # the spacecraft ephemeris locations.

        # Convert the GSE(t) Cartesian coordinates from kilometers to R_sun.
        x = sc_data["XYZ_GSE"][:, 0]/R_SUN_KILOMETERS
        y = sc_data["XYZ_GSE"][:, 1]/R_SUN_KILOMETERS
        z = sc_data["XYZ_GSE"][:, 2]/R_SUN_KILOMETERS
        t = sc_data["Epoch_bin"]

        # Create astropy.coordinates.SkyCoord objects for each GSE(t)
        # position and time.
        c = SkyCoord(
            x*u.Rsun, y*u.Rsun, z*u.Rsun,
            frame=frames.GeocentricSolarEcliptic, obstime=t,
            representation_type="cartesian"
        )

        # Convert the MJDc from a float MJD to a UTC datetime.
        obstime = kaiTools.MJD2UT(MJDc)

        # Create the gamhelio frame (GH(MJDc)) for this data, which is a
        # modified Heliographic Stonyhurst (HGS) frame in SunPy. The frame
        # is defined at the date of MJDc from the WSA solar wind file used
        # for the initial conditions of the simulation. Additionally, the
        # GH(MJDc) frame has the x- and y-axes reversed relative to the HGS
        # frame.
        gh_frame = frames.HeliographicStonyhurst(obstime=obstime)

        # Convert the GSE(t) Cartesian positions to HGS(MJDc) coordinates.
        # As a SkyCoord object, the converted coordinates are available
        # in a variety of coordinate systems.
        c = c.transform_to(gh_frame)

        # Save the HGS(MJDc) Cartesian coordinates as GH(MJDc) coordinates.
        # These variable names (X, Y, Z) are the same as those used in
        # the gamhelio output files.
        sc_data["Ephemeris_time"] = dm.dmarray(t)
        sc_data["X"] = dm.dmarray(-c.cartesian.x)
        sc_data["Y"] = dm.dmarray(-c.cartesian.y)
        sc_data["Z"] = dm.dmarray(c.cartesian.z)

    elif cdaweb_coordinate_system == "HGI":

        # Fetch the value of 1 AU in kilometers.
        AU_km = u.Quantity(1*u.astrophys.AU, u.km).value

        # Fetch the value of 1 Rsun in kilometers.
        Rsun_km = u.Quantity(1*u.Rsun, u.km).value

        # Compute the conversion factor from AU to Rsun.
        Rsun_per_AU = AU_km/Rsun_km

        # Fetch time, HGI latitude, longitude, and radius (convert to Rsun).
        # Note that these vaariables have different names for different
        # spacecraft.
        if "HGI_LAT" in sc_metadata["Ephem"]["Data"]:
            t = sc_data["Epoch_bin"]
            lat = sc_data["HGI_LAT"]
            lon = sc_data["HGI_LON"]
            # Convert radius to Rsun.
            rad = sc_data["RAD_AU"]*Rsun_per_AU
        elif "heliographicLatitude" in sc_metadata["Ephem"]["Data"]:
            t = sc_data["Epoch_bin"]
            lat = sc_data["heliographicLatitude"]
            lon = sc_data["heliographicLongitude"]
            # Convert radius to Rsun.
            rad = sc_data["radialDistance"]*Rsun_per_AU
        else:
            raise TypeError("Unexpected HGI variable names: %s" %
                            sc_metadata["Ephem"]["Data"])

        # Create SkyCoord objects for each HGI(t) position and time.
        # Note HGI = Heliographic Inertial = Heliocentric Inertial
        c = SkyCoord(
            lon*u.deg, lat*u.deg, rad*u.Rsun,
            frame=frames.HeliocentricInertial, obstime=t,
            representation_type="spherical"
        )

        # Create the HGS(t0) coordinate frame.
        mjdc_frame = frames.HeliographicStonyhurst(obstime=kaiTools.MJD2UT(MJDc))

        # Convert the HGI(t) spherical (lon, lat, radius) positions to HGS(MJDc)).
        c = c.transform_to(mjdc_frame)

        # Save the HGS(MJDc) Cartesian coordinates as GH(MJDc) coordinates.
        # These variable names (X, Y, Z) are the same as those used in
        # the gamhelio output files.
        sc_data["Ephemeris_time"] = dm.dmarray(t)
        sc_data["X"] = dm.dmarray(-c.cartesian.x)
        sc_data["Y"] = dm.dmarray(-c.cartesian.y)
        sc_data["Z"] = dm.dmarray(c.cartesian.z)

    else:
        raise KeyError("Unknown ephemeris coordinate system: %s!" %
                       cdaweb_coordinate_system)

    # At this point, sc_data contains all of the ephemeris times, and the
    # Cartesian ephemeris positions as (X, Y, Z) in the GH(MJDc) frame.


def ingest_cdaweb_speed(sc_data, sc_metadata, MJDc, verbose=False, debug=False):
    """Convert CDAWeb speed to gamhelio format.

    Convert CDAWeb speed to gamhelio format. If needed, compute
    the radial component of the speed in the GH(MJDc) frame.

    Parameters
    ----------
    sc_data : dm.SpaceData
        SpaceData object for all data as originally returned from CDAWeb.
    sc_metadata : dict
        Spacecraft descriptive information for the heliosphere spacecraft YML
        file.
    MJDc : float
        Value of MJDc attribute from gamhelio result files.
    verbose : bool, default False
        Set to True for printing verbose progress messages.
    debug : bool, default False
        Set to True for printing debugging messages.

    Returns
    -------
    None
    """
    # If this variable was not found, skip ingest.
    if not sc_metadata["Speed"]["Data"] in sc_data:
        return

    # Fetch the dataset and name of the variable.
    cdaweb_dataset_name = sc_metadata["Speed"]["Id"]
    cdaweb_variable_name = sc_metadata["Speed"]["Data"]

    # Add the radial component of the *spacecraft-measured* speed as
    # a new variable.
    if cdaweb_variable_name == "Vp":
        # The proton velocity is the desired radial velocity.
        sc_data["Speed"] = dm.dmarray(
            sc_data["Vp"],
            attrs = {
                "UNITS": "km/s",
                "CATDESC": "Radial speed",
                "FIELDNAM": "Radial speed",
                "AXISLABEL": "Vr"
            }
        )
    elif cdaweb_variable_name == "plasmaSpeed":
        sc_data["Speed"] = dm.dmarray(
            sc_data["plasmaSpeed"],
            attrs = {
                "UNITS": "km/s",
                "CATDESC": "Radial speed",
                "FIELDNAM": "Radial speed",
                "AXISLABEL": "Vr"
            }
        )
    elif cdaweb_variable_name == "VR":
        sc_data["Speed"] = dm.dmarray(
            sc_data["VR"],
            attrs = {
                "UNITS": "km/s",
                "CATDESC": "Radial speed",
                "FIELDNAM": "Radial speed",
                "AXISLABEL": "Vr"
            }
        )
    elif cdaweb_variable_name == "V_RTN":
        # This is a set of Vr, Vt, Vn components.
        sc_data["Speed"] = dm.dmarray(
            sc_data["V_RTN"][:, 0],
            attrs = {
                "UNITS": "km/s",
                "CATDESC": "Radial speed",
                "FIELDNAM": "Radial speed",
                "AXISLABEL": "Vr"
            }
        )
    else:
        raise TypeError("Unexpected variable: dataset %s, "
                        "variable %s!" %
                        (cdaweb_dataset_name, cdaweb_variable_name))


def ingest_cdaweb_magnetic_field(sc_data, sc_metadata, MJDc, verbose=False, debug=False):
    """Convert CDAWeb magnetic field to gamhelio format.

    Convert CDAWeb magnetic field to gamhelio format. If needed, compute
    the radial component of the magnetic field in the GH(MJDc) frame.

    Parameters
    ----------
    sc_data : dm.SpaceData
        SpaceData object for all data as originally returned from CDAWeb.
    sc_metadata : dict
        Spacecraft descriptive information for the heliosphere spacecraft YML
        file.
    MJDc : float
        Value of MJDc attribute from gamhelio result files.
    verbose : bool, default False
        Set to True for printing verbose progress messages.
    debug : bool, default False
        Set to True for printing debugging messages.

    Returns
    -------
    None
    """
    # If this variable was not found, skip ingest.
    if isinstance(sc_metadata["MagneticField"]["Data"], list):
        cdaweb_variable_name = sc_metadata["MagneticField"]["Data"][0]
    else:
        cdaweb_variable_name = sc_metadata["MagneticField"]["Data"]
    if not cdaweb_variable_name in sc_data:
        return

    # Fetch the dataset and name of the variable.
    cdaweb_dataset_name = sc_metadata["MagneticField"]["Id"]
    cdaweb_variable_name = sc_metadata["MagneticField"]["Data"]

    # Add the radial component of the *spacecraft-measured* magnetic field as
    # a new variable.
    magnetic_field_coordinate_frame = sc_metadata["MagneticField"]["CoordSys"]
    if magnetic_field_coordinate_frame == "GSE":
        if cdaweb_variable_name == "BGSEc":
            # The magnetic fiels was returned as a 2-D array of shape (n, 3)
            # containing (Bx, By, Bz) in units of nT.
            # The radial component of B in GSE is just the negative
            # x-component, since the GSE +x axis points in the opposite
            # direction as the GH(MHDc) +x axis.
            sc_data["Br"] = dm.dmarray(
                -sc_data["BGSEc"][:, 0],
                attrs = {
                    "UNITS": "nT",
                    "CATDESC": "Radial magnetic field",
                    "FIELDNAM": "Radial magnetic field",
                    "AXISLABEL": "Br"
                }
            )
    elif magnetic_field_coordinate_frame == "RTN":
        # Just use the radial component.
        if "Br" in sc_data:
            sc_data["Br"] = dm.dmarray(
                sc_data["BR"][:],
                attrs = {
                    "UNITS": "nT",
                    "CATDESC": "Radial magnetic field",
                    "FIELDNAM": "Radial magnetic field",
                    "AXISLABEL": "Br"
                }
            )
        elif "BR" in sc_data:
            sc_data["Br"] = dm.dmarray(
                sc_data["BR"][:],
                attrs = {
                    "UNITS": "nT",
                    "CATDESC": "Radial magnetic field",
                    "FIELDNAM": "Radial magnetic field",
                    "AXISLABEL": "Br"
                }
            )
        elif "B_RTN" in sc_data:
            sc_data["Br"] = dm.dmarray(
                sc_data["B_RTN"][:, 0],
                attrs = {
                    "UNITS": "nT",
                    "CATDESC": "Radial magnetic field",
                    "FIELDNAM": "Radial magnetic field",
                    "AXISLABEL": "Br"
                }
            )
        else:
            raise TypeError("Unexpected variable: dataset %s, "
                "variable %s!" %
                (cdaweb_dataset_name, cdaweb_variable_name))

    else:
        raise TypeError("Unexpected magnetic coordinate frame: "
                        "dataset %s, variable %s!" %
                        (cdaweb_dataset_name, cdaweb_variable_name))


def ingest_cdaweb_density(sc_data, sc_metadata, MJDc, verbose=False, debug=False):
    """Convert CDAWeb density to gamhelio format.

    Convert CDAWeb density to gamhelio format.

    Parameters
    ----------
    sc_data : dm.SpaceData
        SpaceData object for all data as originally returned from CDAWeb.
    sc_metadata : dict
        Spacecraft descriptive information for the heliosphere spacecraft YML
        file.
    MJDc : float
        Value of MJDc attribute from gamhelio result files.
    verbose : bool, default False
        Set to True for printing verbose progress messages.
    debug : bool, default False
        Set to True for printing debugging messages.

    Returns
    -------
    None
    """
    # Fetch the dataset and name of the variable.
    cdaweb_dataset_name = sc_metadata["Density"]["Id"]
    cdaweb_variable_name = sc_metadata["Density"]["Data"]

    # Add the *spacecraft-measured* number density as a new variable.
    if cdaweb_variable_name == "Np":
        # The proton density is the desired density.
        sc_data["Density"] = dm.dmarray(
            sc_data["Np"],
            attrs = {
                "UNITS": "#/cc",
                "CATDESC": "Number density",
                "FIELDNAM": "Number density",
                "AXISLABEL": "N"
            }
        )
    elif cdaweb_variable_name == "plasmaDensity":
        # The proton density is the desired density.
        sc_data["Density"] = dm.dmarray(
            sc_data["plasmaDensity"],
            attrs = {
                "UNITS": "#/cc",
                "CATDESC": "Number density",
                "FIELDNAM": "Number density",
                "AXISLABEL": "N"
            }
        )
    elif cdaweb_variable_name == "protonDensity":
        # The proton density is the desired density.
        sc_data["Density"] = dm.dmarray(
            sc_data["protonDensity"],
            attrs = {
                "UNITS": "#/cc",
                "CATDESC": "Number density",
                "FIELDNAM": "Number density",
                "AXISLABEL": "N"
            }
        )
    elif cdaweb_variable_name == "N":
        # The proton density is the desired density.
        sc_data["Density"] = dm.dmarray(
            sc_data["N"],
            attrs = {
                "UNITS": "#/cc",
                "CATDESC": "Number density",
                "FIELDNAM": "Number density",
                "AXISLABEL": "N"
            }
        )
    else:
        raise TypeError("Unexpected variable: dataset %s, "
                        "variable %s!" %
                        (cdaweb_dataset_name, cdaweb_variable_name))


def ingest_cdaweb_temperature(sc_data, sc_metadata, MJDc, verbose=False, debug=False):
    """Convert CDAWeb temperature to gamhelio format.

    Convert CDAWeb temperature to gamhelio format.

    Parameters
    ----------
    sc_data : dm.SpaceData
        SpaceData object for all data as originally returned from CDAWeb.
    sc_metadata : dict
        Spacecraft descriptive information for the heliosphere spacecraft YML
        file.
    MJDc : float
        Value of MJDc attribute from gamhelio result files.
    verbose : bool, default False
        Set to True for printing verbose progress messages.
    debug : bool, default False
        Set to True for printing debugging messages.

    Returns
    -------
    None
    """
    # Fetch the dataset and name of the variable.
    cdaweb_dataset_name = sc_metadata["Temperature"]["Id"]
    cdaweb_variable_name = sc_metadata["Temperature"]["Data"]

    # Add the *spacecraft-measured* number density as a new variable.
    if cdaweb_variable_name == "Tpr":
        # The proton density is the desired density.
        sc_data["Temperature"] = dm.dmarray(
            sc_data["Tpr"],
            attrs = {
                "UNITS": "K",
                "CATDESC": "Temperature",
                "FIELDNAM": "Temperature",
                "AXISLABEL": "T"
            }
        )
    elif cdaweb_variable_name == "plasmaTemp":
        sc_data["Temperature"] = dm.dmarray(
            sc_data["plasmaTemp"],
            attrs = {
                "UNITS": "K",
                "CATDESC": "Temperature",
                "FIELDNAM": "Temperature",
                "AXISLABEL": "T"
            }
        )
    elif cdaweb_variable_name == "protonTemp":
        sc_data["Temperature"] = dm.dmarray(
            sc_data["protonTemp"],
            attrs = {
                "UNITS": "K",
                "CATDESC": "Temperature",
                "FIELDNAM": "Temperature",
                "AXISLABEL": "T"
            }
        )
    elif cdaweb_variable_name == "T":
        sc_data["Temperature"] = dm.dmarray(
            sc_data["T"],
            attrs = {
                "UNITS": "K",
                "CATDESC": "Temperature",
                "FIELDNAM": "Temperature",
                "AXISLABEL": "T"
            }
        )
    else:
        raise TypeError("Unexpected variable: dataset %s, "
                        "variable %s!" %
                        (cdaweb_dataset_name, cdaweb_variable_name))


def ingest_helio_cdaweb_data(sc_id, sc_data, sc_metadata, MJDc, verbose=False, debug=False):
    """Ingest CDAWeb data to gamhelio format.

    Convert CDAWeb data to gamhelio format.

    Parameters
    ----------
    sc_data : dm.SpaceData
        SpaceData object for all data as originally returned from CDAWeb.
    sc_metadata : dict
        Spacecraft descriptive information for the heliosphere spacecraft YML
        file.
    MJDc : float
        Value of MJDc attribute from gamhelio result files.
    verbose : bool, default False
        Set to True for printing verbose progress messages.
    debug : bool, default False
        Set to True for printing debugging messages.

    Returns
    -------
    None
    """
    # Ingest the ephemeris.
    ingest_cdaweb_ephemeris(
        sc_data, sc_metadata, MJDc,
        verbose=verbose, debug=debug
    )

    # Ingest the speed measurements, if available.
    if sc_metadata["Speed"]["Data"] in sc_data:
        ingest_cdaweb_speed(
            sc_data, sc_metadata, MJDc,
            verbose=verbose, debug=debug
        )

    # Ingest the magnetic field measurements.
    # if sc_metadata["MagneticField"]["Data"] in sc_data:
    ingest_cdaweb_magnetic_field(
        sc_data, sc_metadata, MJDc,
        verbose=verbose, debug=debug
    )

    # Ingest the density measurements.
    if sc_metadata["Density"]["Data"] in sc_data:
        ingest_cdaweb_density(
            sc_data, sc_metadata, MJDc,
            verbose=verbose, debug=debug
        )

    # Ingest the temperature measurements.
    if sc_metadata["Temperature"]["Data"] in sc_data:
        ingest_cdaweb_temperature(
            sc_data, sc_metadata, MJDc,
            verbose=verbose, debug=debug
        )


def create_sctrack_helio_trajectory_file(
    sc_data, sc_metadata, sc_id, mjd0, sec0, gamhelio_results_directory,
    ftag, numSegments, mjdc
):
    """Create the input trajectory file needed by the sctrack.x interpolator.

    Create the input trajectory file needed by the sctrack.x interpolator.

    Parameters
    ----------
    XXX

    Returns
    -------
    XXX
    """
    # Fetch the times and coordinates of the trajectory. They should already
    # be in the Cartesian GH(MJDc) frame.
    t = sc_data["Ephemeris_time"]
    x = sc_data["X"]
    y = sc_data["Y"]
    z = sc_data["Z"]

    # Compute the elapsed times in seconds since the start of the gamhelio
    # simulation for each ephemeris position. sctrack.x needs this value in
    # order to perform the interpolation of gamhelio model output to the
    # ephemeris points.
    elapsed = [(tt - t[0]).total_seconds() for tt in t]

    # Create the HDF5 file containing the spacecraft data transformed to the
    # GH(MJDc) frame and elapsed time.
    trajectory_path = os.path.join(
        gamhelio_results_directory, sc_id + ".sc.h5"
    )
    with h5py.File(trajectory_path, "w") as hf:
        hf.create_dataset("T", data=elapsed)
        hf.create_dataset("X", data=x)
        hf.create_dataset("Y", data=y)
        hf.create_dataset("Z", data=z)

    # Return the path to the trajectory file.
    return trajectory_path


def create_sctrack_helio_xml(
    output_directory, run_id, spacecraft_id,
    trajectory_file, num_parallel_segments=1
):
    """Generate heliosphere XML input file for sctrack.x.

    Generate heliosphere XML input file for sctrack.x. This file provides
    the information required by sctrack.x to interpolate gamhelio model
    results to the times and locations of the spacecraft data.

    Parameters
    ----------
    output_directory: str
        Path to directory containing gamhelio model output.
    run_id: str
        Identifying tag used in gamhelio result filenames.
    spacecraft_id: str
        Name string (no spaces allowed) for spacecraft used in trajectory.
    trajectory_file: str
        Name of file containing spacecraft trajectory.
    num_parallel_segments: int, default 1
        Number of threads for sctrack.x to use.

    Returns
    -------
    root: minidom.Document
        Root document object for XML.
    """

    # Determine if the gamhelio results were generated by an MPI run, and
    # the organization of the MPI ranks.
    (file_name, isMPI, Ri, Rj, Rk) = kaiTools.getRunInfo(output_directory, run_id)

    # Create the XML Document.
    xml_doc = minidom.Document()

    # Create the top-level <Kaiju> element.
    kaiju_el = xml_doc.createElement("Kaiju")
    xml_doc.appendChild(kaiju_el)

    # Create the <Chimp> element.
    chimp_el = xml_doc.createElement("Chimp")
    kaiju_el.appendChild(chimp_el)

    # Create the <sim> element.
    sim_el = xml_doc.createElement("sim")
    sim_el.setAttribute("runid", spacecraft_id)
    chimp_el.appendChild(sim_el)

    # Create the <fields> element.
    fields_el = xml_doc.createElement("fields")
    fields_el.setAttribute("doMHD", "T")
    fields_el.setAttribute("grType", "SPH")
    fields_el.setAttribute("ebfile", run_id)
    if isMPI:
        fields_el.setAttribute("isMPI", "T")
    chimp_el.appendChild(fields_el)

    # If the results were created by an MPI run, create the <parallel>
    # element.
    if isMPI:
        parallel_el = xml_doc.createElement("parallel")
        parallel_el.setAttribute("Ri", "%d" % Ri)
        parallel_el.setAttribute("Rj", "%d" % Rj)
        parallel_el.setAttribute("Rk", "%d" % Rk)
        chimp_el.appendChild(parallel_el)

    # Create the <units> element.
    units_el = xml_doc.createElement("units")
    units_el.setAttribute("uid", "HELIO")
    chimp_el.appendChild(units_el)

    # Create the <trajectory> element.
    trajectory_el = xml_doc.createElement("trajectory")
    trajectory_el.setAttribute("H5Traj", trajectory_file)
    trajectory_el.setAttribute("doSmooth", "F")
    chimp_el.appendChild(trajectory_el)

    # Create the <domain> element.
    domain_el = xml_doc.createElement("domain")
    domain_el.setAttribute("dtype", "SPH")
    domain_el.setAttribute("rmin", "%s" % R_MIN_HELIO)
    domain_el.setAttribute("rmax", "%s" % R_MAX_HELIO)
    chimp_el.appendChild(domain_el)

    # If data interpolation will be run in parallel, create the <parintime>
    # element.
    if num_parallel_segments > 1:
        parintime_el = xml_doc.createElement("parintime")
        parintime_el.setAttribute("NumB", "%d" % num_parallel_segments)
        chimp_el.appendChild(parintime_el)

    # Return the entire XML document object.
    return xml_doc


def create_sctrack_helio_input_files(
    data, scDic, scId, mjd0, sec0, gamhelio_results_directory,
    ftag, numSegments, mjdc
):
    """Create the input files needed by the sctrack.x interpolator.

    Create the input files needed by the sctrack.x interpolator.

    Parameters
    ----------
    XXX

    Returns
    -------
    XXX
    """
    # Create the HDF5 file containing the trajectory data.
    trajectory_path = create_sctrack_helio_trajectory_file(
        data, scDic, scId, mjd0, sec0, gamhelio_results_directory,
        ftag, numSegments, mjdc
    )

    # Extract the trajectory file name from the path.
    trajectory_file = os.path.split(trajectory_path)[-1]

    # Create the XML describing the required interpolations.
    chimp_xml = create_sctrack_helio_xml(
        gamhelio_results_directory, ftag, scId,
        trajectory_file, num_parallel_segments=0
    )

    # Write the XML to a file.
    xml_path = os.path.join(gamhelio_results_directory, scId + ".xml")
    with open(xml_path, "w") as f:
        f.write(chimp_xml.toprettyxml())

    # Return the paths to the HDF5 and XML files.
    return (trajectory_path, xml_path)


def ingest_interpolated_radial_velocity(h5file, sc_data):
    """Ingest the interpolated radial velocity and transform as needed.

    Ingest the interpolated radial velocity and transform as needed.

    Parameters
    ----------
    XXX

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Fetch the Cartesian components of the spacecraft position in the
    # GH(MJDc) frame used by the gamhelio model. Each is shape (n,).
    # Then create a shape (n, 3) array for further use.
    # Then compute the radius.
    X = h5file["X"]
    Y = h5file["Y"]
    Z = h5file["Z"]
    R = np.vstack([X, Y, Z]).T
    radius = np.sqrt(X[:]**2 + Y[:]**2 + Z[:]**2)

    # Fetch the Cartesian components of the interpolated model solar wind
    # velocity in the GH(MJDc) frame.
    Vx = h5file["Vx"]
    Vy = h5file["Vy"]
    Vz = h5file["Vz"]

    # Count the number of interpolated velocities.
    n = len(Vx)

    # Create an array of the interpolated velocity as a shape (n, 3) array.
    V = np.vstack([Vx, Vy, Vz]).T

    # Compute the radial component of the interpolated velocity.
    Vr = np.empty(n)
    for i in range(n):
        Vr[i] = V[i].dot(R[i])/radius[i]

    # Add the interpolated model solar wind radial velocity as a new variable.
    sc_data["GAMHELIO_Speed"] = dm.dmarray(
        Vr,
        attrs = {
            "UNITS": "km/s",
            "CATDESC": "Radial velocity",
            "FIELDNAM": "Radial velocity",
            "AXISLABEL": "Vr"
        }
    )


def ingest_interpolated_radial_magnetic_field(h5file, sc_data):
    """Ingest the interpolated radial magnetic field and transform as needed.

    Ingest the interpolated radial magnetic field and transform as needed.

    Parameters
    ----------
    XXX

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Fetch the Cartesian components of the spacecraft position in the
    # GH(MJDc) frame used by the gamhelio model. Each is shape (n,).
    # Then create a shape (n, 3) array for further use.
    # Then compute the radius.
    X = h5file["X"]
    Y = h5file["Y"]
    Z = h5file["Z"]
    R = np.vstack([X, Y, Z]).T
    radius = np.sqrt(X[:]**2 + Y[:]**2 + Z[:]**2)

    # Count the number of interpolated values.
    n = len(X)

    # Fetch the Cartesian components of the interpolated model magnetic field
    # in the GH(t0) frame. These variables were interpolated from the
    # gamhelio results, in the GH(MJDc) frame.
    Bx = h5file["Bx"]
    By = h5file["By"]
    Bz = h5file["Bz"]

    # Create an array of the interpolated model magnetic field components as
    # a shape (n, 3) array.
    B = np.vstack([Bx, By, Bz]).T

    # Compute the radial component of the interpolated model magnetic field
    # in the GH(MJDc) frame.
    Br = np.empty(n)
    for i in range(n):
        Br[i] = B[i].dot(R[i])/radius[i]

    # Add the radial component of the interpolated model magnetic field as a
    # new variable. The negative of Br is needed to reverse the x-axis from
    # GH(t0) to GSE(t).
    sc_data["GAMHELIO_Br"] = dm.dmarray(
        Br,
        attrs = {
            "UNITS": "nT",
            "CATDESC": "Radial magnetic field",
            "FIELDNAM": "Radial magnetic field",
            "AXISLABEL": "Br"
        }
    )


def ingest_interpolated_number_density(h5file, sc_data):
    """Ingest the interpolated number density and transform as needed.

    Ingest the interpolated number density and transform as needed.

    Parameters
    ----------
    XXX

    Returns
    -------
    None

    Raises
    ------
    None
    """
    density = h5file["D"]
    sc_data["GAMHELIO_Density"] = dm.dmarray(
        density[:],
        attrs={
            "UNITS": "#/cc",
            "CATDESC": "Density",
            "FIELDNAM": "Density",
            "AXISLABEL": "n"
        }
    )


def ingest_interpolated_temperature(h5file, sc_data):
    """Ingest the interpolated temperature and transform as needed.

    Ingest the interpolated temperature and transform as needed.

    Parameters
    ----------
    XXX

    Returns
    -------
    None

    Raises
    ------
    None
    """
    # Read the interpolated model pressure and number density.
    pressure = h5file["P"]
    density = h5file["D"]

    # Compute the gamhelio temperature from interpolated model pressure
    # and interpolated model density using the ideal gas law, and add as a
    # new variable.
    # Pressure is in cgs units (erg/cm**3).
    # Density is in 1/cm**3.
    # The CGS Boltzmann constant kbltz is erg/K.
    # The factor of 2 is needed since we have a neutral 2-component plasma.
    # The factor of 1e-8 converts nPa to erg/cm**3.
    temperature = pressure[:]*1e-8/(2*kaipy.kdefs.kbltz*density[:])
    sc_data["GAMHELIO_Temperature"] = dm.dmarray(
        temperature[:],
        attrs={
            "UNITS": "K",
            "CATDESC": "Temperature",
            "FIELDNAM": "Temperature",
            "AXISLABEL": "T"
        }
    )


def ingest_interpolated_inDomain(h5file, sc_data):
    """Ingest the interpolated inDomain flag.

    Ingest the interpolated inDomain flag.

    Parameters
    ----------
    XXX

    Returns
    -------
    None

    Raises
    ------
    None
    """
    inDom = h5file["inDom"]
    sc_data["GAMHELIO_inDom"] = dm.dmarray(
        inDom[:],
        attrs={
            "UNITS": inDom.attrs["Units"],
            "CATDESC": "In GAMERA Domain",
            "FIELDNAM": "InDom",
            "AXISLABEL": "In Domain"
        }
    )


def ingest_interpolated_variables(
    sc_data, sc_metadata, interpolated_results_path
):
    """Copy the interpolated model results and transform as needed.

    Copy the variables in the input HDF5 file into new variables in the HDF5
    file with GAMERA_-prefixed descriptive names, and more metadata.

    Convert vector values from the gamhelio frame (GH(t0)) to the frame of
    the spacecraft.

    Parameters
    ----------
    sc_data: spacepy.datamodel.SpaceData
        All of the spacecraft and interpolated model results so far.
    sc_metadata: dict
        Spacecraft descriptive information.
    interpolated_results_path: str
        Path to HDF5 file containing gamhelio model results interpolated to the
        spacecraft positions, all in the GH(t0) frame.

    Returns
    -------
    None

    Raises
    ------
    TypeError
        If target frame is not GSE.
    """
    # Open the HDF5 file containing the gamhelio model results interpolated
    # to the spacecraft positions.
    h5file = h5py.File(interpolated_results_path, "r")

    # Ingest the interpolated radial velocity.
    ingest_interpolated_radial_velocity(h5file, sc_data)

    # Ingest the interpolated radial magnetic field.
    ingest_interpolated_radial_magnetic_field(h5file, sc_data)

    # Ingest the interpolated number density.
    ingest_interpolated_number_density(h5file, sc_data)

    # Ingest the interpolated temperature.
    ingest_interpolated_temperature(h5file, sc_data)

    # Ingest the in-domain flag.
    ingest_interpolated_inDomain(h5file, sc_data)

    # At this point, we have added the interpolated model values to
    # the data object, along with extra metadata needed for comparison
    # plotting. The names of the interpolated variables are prefixed
    # with "GAMHELIO_".


def interpolate_gamhelio_results_to_trajectory(
    sc_data, sc_metadata, sc_id, first_MJD, first_elapsed_seconds,
    gamhelio_results_directory, run_id, sctrack_cmd, num_segments, keep, MJDc
):
    """Interpolate gamhelio data to the spacecraft trajectory.

    Interpolate gamhelio data to the spacecraft trajectory. This code assumes
    T, X, Y, Z are already in the GH(MJDc) frame.

    Parameters
    ----------
    XXX

    Returns
    -------
    None

    Raises
    ------
    XXX
    """
    # Create the file containing the spacecraft trajectory, and the associated
    # XML file read by sctrack.x.
    (trajectory_path, sctrack_xml_path) = create_sctrack_helio_input_files(
        sc_data, sc_metadata, sc_id,
        first_MJD, first_elapsed_seconds,
        gamhelio_results_directory, run_id,
        num_segments, MJDc
    )

    # Perform the interpolation.
    if num_segments == 1:
        # Perform a serial interpolation.
        sctrack = subprocess.run(
            [sctrack_cmd, sctrack_xml_path],
            cwd=gamhelio_results_directory,
            capture_output=True,
            text=True
        )

        # Save the interpolator output in the results directory.
        with open(os.path.join(gamhelio_results_directory,
                               "%s_sctrack.out" % sc_id), "w") as f:
            f.write(sctrack.stdout)

    else:
        # Perform a parallel interpolation.
        raise TypeError("Parallel interpolation currently unsupported!")
        # process = []
        # for seg in range(1,num_segments+1):
        #     process.append(subprocess.Popen([sctrack_cmd, sctrack_xml_path,str(seg)],
        #                     cwd=gamhelio_results_directory,
        #                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        #                     text=True))
        # for proc in process:
        #     proc.communicate()
        # h5name = mergeFiles(sc_id,gamhelio_results_directory,num_segments)

    # Ingest the interpolated gamhelio data.
    interpolated_results_path = trajectory_path
    ingest_interpolated_variables(sc_data, sc_metadata, interpolated_results_path)

    # Remove temporary files if needed.
    # if not keep:
    #     os.remove(sctrack_xml_path)
        # if num_segments > 1:
        #     h5parts = os.path.join(gamhelio_results_directory,sc_id+'.*.sc.h5')
        #     subprocess.run(['rm',h5parts])


def write_helio_error_report(error_file_path, sc_id, sc_data):
    """Save an error report for the current data.

    Compute and save an error report for the current data.

    Parameters
    ----------
    errorName: str
        Path to file to hold error report.
    scId: str
        ID string for spacecraft.
    data: spacepy.datamodel.SpaceData
        The current spacecraft and model data

    Returns
    -------
    None
    """
    # Determine which data are available for error computations.
    keysToCompute = ["Speed", "Br", "Density", "Temperature"]

    # Compute and save the error report for each variable.
    with open(error_file_path, "w") as f:
        for key in keysToCompute:
            if key in sc_data:
                maskedData = np.ma.masked_where(
                    sc_data["GAMHELIO_inDom"][:] == 0.0, sc_data[key][:]
                )
                maskedGamera = np.ma.masked_where(
                    sc_data["GAMHELIO_inDom"][:] == 0.0, sc_data["GAMHELIO_" + key][:]
                )
                MAE, MSE, RMSE, MAPE, RSE, PE = computeErrors(
                    maskedData, maskedGamera
                )
                f.write(f"Errors for: {key}\n")
                f.write(f"MAE: {MAE}\n")
                f.write(f"MSE: {MSE}\n")
                f.write(f"RMSE: {RMSE}\n")
                f.write(f"MAPE: {MAPE}\n")
                f.write(f"RSE: {RSE}\n")
                f.write(f"PE: {PE}\n")
