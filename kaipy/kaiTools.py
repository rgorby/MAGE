import numpy as np
import datetime
import os
import glob
import sys
import subprocess
from xml.dom import minidom

#import spacepy and cdasws
import spacepy
from spacepy.coordinates import Coords
from spacepy.time import Ticktock
import spacepy.datamodel as dm
from cdasws import CdasWs
from astropy.time import Time

#import hdf5
import h5py

def MJD2UT(mjd):
	astroT = Time(mjd,format='mjd').iso
	utall = []
	for ut in astroT:
		utall.append(datetime.datetime.strptime(ut,'%Y-%m-%d %H:%M:%S.%f'))

	return utall

def genSCXML(fdir,ftag,
	scid="sctrack_A",h5traj="sctrack_A.h5",numSegments=1):

	(fname,isMPI,Ri,Rj,Rk) = getRunInfo(fdir,ftag)
	root = minidom.Document()
	xml = root.createElement('Chimp')
	root.appendChild(xml)
	scChild = root.createElement("sim")
	scChild.setAttribute("runid",scid)
	xml.appendChild(scChild)
	fieldsChild = root.createElement("fields")
	fieldsChild.setAttribute("doMHD","T")
	fieldsChild.setAttribute("grType","LFM")
	fieldsChild.setAttribute("ebfile",ftag)
	if isMPI:
		fieldsChild.setAttribute("isMPI","T")
	xml.appendChild(fieldsChild)
	if isMPI:
		parallelChild = root.createElement("parallel")
		parallelChild.setAttribute("Ri","%d"%Ri)
		parallelChild.setAttribute("Rj","%d"%Rj)
		parallelChild.setAttribute("Rk","%d"%Rk)
		xml.appendChild(parallelChild)
	unitsChild = root.createElement("units")
	unitsChild.setAttribute("uid","EARTH")
	xml.appendChild(unitsChild)
	trajChild = root.createElement("trajectory")
	trajChild.setAttribute("H5Traj",h5traj)
	trajChild.setAttribute("doSmooth","F")
	xml.appendChild(trajChild)
	if numSegments > 1:
		parInTimeChild = root.createElement("parintime")
		parInTimeChild.setAttribute("NumB","%d"%numSegments)
		xml.appendChild(parInTimeChild)

	return root

def getRunInfo(fdir,ftag):
	idStr = "_0000_0000_0000.gam.h5"
	isMPI = False
	Ri = 0
	Rj = 0
	Rk = 0
	fOld = os.path.join(fdir,ftag+'.h5')
	fNew = os.path.join(fdir,ftag+'.gam.h5')
	try:
		if (os.path.exists(fOld)):
			return fOld,isMPI,Ri,Rj,Rk
		if (os.path.exists(fNew)):
			return fNew,isMPI,Ri,Rj,Rk
		fIns = glob.glob(os.path.join(fdir,ftag)+'_*'+idStr)
		if (len(fIns) > 1):
			raise ValueError('Should not find more that one parallel file')
		if (len(fIns) == 0):
			raise ValueError('No MPI database found')
		else:
			isMPI = True
			fName = fIns[0]
			Ns = [int(s) for s in fName.split('_') if s.isdigit()]
			Ri = Ns[-5]
			Rj = Ns[-4]
			Rk = Ns[-3]
			return fName,isMPI,Ri,Rj,Rk
	except ValueError as ve:
			print(ve)
			sys.exit()

def getScIds():
	scDic = {}
	scDic['GOES11'] = {'ephemId':'GOES11_K0_MAG',
						'ephemData':'SC_pos_sm',
						'ephemCoordSys':'SM',
						'denId':None,
						'denData':None,
						'presId': None,
						'presData': None,
						'velId': None,
						'velData':None,
						'velCoordSys':None,
						'magId': 'GOES11_K0_MAG',
						'magData': 'B_GSM_c',
						'magCoordSys': 'GSM'
						}
	scDic['GOES12'] = {'ephemId':'GOES12_K0_MAG',
						'ephemData':'SC_pos_sm',
						'ephemCoordSys':'SM',
						'denId':None,
						'denData':None,
						'presId': None,
						'presData': None,
						'velId': None,
						'velData':None,
						'velCoordSys':None,
						'magId': 'GOES12_K0_MAG',
						'magData': 'B_GSM_c',
						'magCoordSys': 'GSM'
						}
	scDic['GEOTAIL'] = {'ephemId':'GE_K0_MGF',
						'ephemData':'POS',
						'ephemCoordSys':'GSE',
						'denId':'GE_H0_CPI',
						'denData':'SW_P_Den',
						'presId': None,
						'presData': None,
						'velId': 'GE_H0_CPI',
						'velData':'SW_Vc',
						'velCoordSys':'GSE',
						'magId': 'GE_K0_MGF',
						'magData': 'IB_vector',
						'magCoordSys': 'GSE'
						}
	scDic['RBSPA'] = {'ephemId':'RBSP-A_MAGNETOMETER_1SEC-GSM_EMFISIS-L3',
						'ephemData':'coordinates',
						'ephemCoordSys':'GSM',
						'denId':'RBSPA_REL04_ECT-HOPE-MOM-L3',
						'denData':'Dens_p_30',
						'presId': None,
						'presData': None,
						'velId': None,
						'velData':None,
						'velCoordSys':None,
						'magId': 'RBSP-A_MAGNETOMETER_1SEC-GSM_EMFISIS-L3',
						'magData': 'Mag',
						'magCoordSys': 'GSM'
						}
	scDic['RBSPB'] = {'ephemId':'RBSP-B_MAGNETOMETER_1SEC-GSM_EMFISIS-L3',
						'ephemData':'coordinates',
						'ephemCoordSys':'GSM',
						'denId':'RBSPB_REL04_ECT-HOPE-MOM-L3',
						'denData':'Dens_p_30',
						'presId': None,
						'presData': None,
						'velId': None,
						'velData':None,
						'velCoordSys':None,
						'magId': 'RBSP-B_MAGNETOMETER_1SEC-GSM_EMFISIS-L3',
						'magData': 'Mag',
						'magCoordSys': 'GSM'
						}
	scDic['CLUSTER1'] = {'ephemId':'C1_CP_FGM_SPIN',
						'ephemData':'sc_pos_xyz_gse__C1_CP_FGM_SPIN',
						'ephemCoordSys':'GSE',
						'denId':'C1_PP_CIS',
						'denData':'N_p__C1_PP_CIS',
						'presId': None,
						'presData': None,
						'velId': 'C1_PP_CIS',
						'velData':'V_p_xyz_gse__C1_PP_CIS',
						'velCoordSys':'GSE',
						'magId': 'C1_CP_FGM_SPIN',
						'magData': 'B_vec_xyz_gse__C1_CP_FGM_SPIN',
						'magCoordSys': 'GSE'
						}
	scDic['CLUSTER2'] = {'ephemId':'C2_CP_FGM_SPIN',
						'ephemData':'sc_pos_xyz_gse__C2_CP_FGM_SPIN',
						'ephemCoordSys':'GSE',
						'denId':'C2_PP_CIS',
						'denData':'N_p__C2_PP_CIS',
						'presId': None,
						'presData': None,
						'velId': 'C2_PP_CIS',
						'velData':'V_p_xyz_gse__C2_PP_CIS',
						'velCoordSys':'GSE',
						'magId': 'C2_CP_FGM_SPIN',
						'magData': 'B_vec_xyz_gse__C2_CP_FGM_SPIN',
						'magCoordSys': 'GSE'
						}
	scDic['CLUSTER3'] = {'ephemId':'C3_CP_FGM_SPIN',
						'ephemData':'sc_pos_xyz_gse__C3_CP_FGM_SPIN',
						'ephemCoordSys':'GSE',
						'denId':'C3_PP_CIS',
						'denData':'N_p__C3_PP_CIS',
						'presId': None,
						'presData': None,
						'velId': 'C3_PP_CIS',
						'velData':'V_p_xyz_gse__C3_PP_CIS',
						'velCoordSys':'GSE',
						'magId': 'C3_CP_FGM_SPIN',
						'magData': 'B_vec_xyz_gse__C3_CP_FGM_SPIN',
						'magCoordSys': 'GSE'
						}
	scDic['CLUSTER4'] = {'ephemId':'C4_CP_FGM_SPIN',
						'ephemData':'sc_pos_xyz_gse__C4_CP_FGM_SPIN',
						'ephemCoordSys':'GSE',
						'denId':'C4_PP_CIS',
						'denData':'N_p__C4_PP_CIS',
						'presId': None,
						'presData': None,
						'velId': 'C4_PP_CIS',
						'velData':'V_p_xyz_gse__C4_PP_CIS',
						'velCoordSys':'GSE',
						'magId': 'C4_CP_FGM_SPIN',
						'magData': 'B_vec_xyz_gse__C4_CP_FGM_SPIN',
						'magCoordSys': 'GSE'
						}
	scDic['THEMISA'] = {'ephemId':'THA_OR_SSC',
						'ephemData':'XYZ_GSM',
						'ephemCoordSys':'GSM',
						'denId':'THA_L2_MOM',
						'denData':'tha_peim_density',
						'presId': 'THA_L2_MOM',
						'presData': 'tha_peim_ptot',
						'velId': 'THA_L2_MOM',
						'velData':'tha_peim_velocity_gsm',
						'velCoordSys':'GSM',
						'magId': 'THA_L2_FGM',
						'magData': 'tha_fgs_gsm',
						'magCoordSys': 'GSM'
						}
	scDic['THEMISB'] = {'ephemId':'THB_OR_SSC',
						'ephemData':'XYZ_GSM',
						'ephemCoordSys':'GSM',
						'denId':'THB_L2_MOM',
						'denData':'thb_peim_density',
						'presId': 'THB_L2_MOM',
						'presData': 'thb_peim_ptot',
						'velId': 'THB_L2_MOM',
						'velData':'thb_peim_velocity_gsm',
						'velCoordSys':'GSM',
						'magId': 'THB_L2_FGM',
						'magData': 'thb_fgs_gsm',
						'magCoordSys': 'GSM'
						}
	scDic['THEMISC'] = {'ephemId':'THC_OR_SSC',
						'ephemData':'XYZ_GSM',
						'ephemCoordSys':'GSM',
						'denId':'THC_L2_MOM',
						'denData':'thc_peim_density',
						'presId': 'THC_L2_MOM',
						'presData': 'thc_peim_ptot',
						'velId': 'THC_L2_MOM',
						'velData':'thc_peim_velocity_gsm',
						'velCoordSys':'GSM',
						'magId': 'THC_L2_FGM',
						'magData': 'thc_fgs_gsm',
						'magCoordSys': 'GSM'
						}
	scDic['THEMISD'] = {'ephemId':'THD_OR_SSC',
						'ephemData':'XYZ_GSM',
						'ephemCoordSys':'GSM',
						'denId':'THD_L2_MOM',
						'denData':'thd_peim_density',
						'presId': 'THD_L2_MOM',
						'presData': 'thd_peim_ptot',
						'velId': 'THD_L2_MOM',
						'velData':'thd_peim_velocity_gsm',
						'velCoordSys':'GSM',
						'magId': 'THD_L2_FGM',
						'magData': 'thd_fgs_gsm',
						'magCoordSys': 'GSM'
						}
	scDic['THEMISE'] = {'ephemId':'THE_OR_SSC',
						'ephemData':'XYZ_GSM',
						'ephemCoordSys':'GSM',
						'denId':'THE_L2_MOM',
						'denData':'the_peim_density',
						'presId': 'THE_L2_MOM',
						'presData': 'the_peim_ptot',
						'velId': 'THE_L2_MOM',
						'velData':'the_peim_velocity_gsm',
						'velCoordSys':'GSM',
						'magId': 'THE_L2_FGM',
						'magData': 'the_fgs_gsm',
						'magCoordSys': 'GSM'
						}
	scDic['MMS1'] = {'ephemId':'MMS1_FGM_SRVY_L2',
						'ephemData':'mms1_fgm_r_gsm_srvy_l2',
						'ephemCoordSys':'GSM',
						'denId':'MMS1_FPI_FAST_L2_DIS-MOMS',
						'denData':'mms1_dis_numberdensity_fast',
						'presId': 'MMS1_FPI_FAST_L2_DIS-MOMS',
						'presData': 'mms1_dis_pres_bg_fast',
						'velId': 'MMS1_FPI_FAST_L2_DIS-MOMS',
						'velData':'mms1_dis_bulkv_gse_fast',
						'velCoordSys':'GSE',
						'magId': 'MMS1_FGM_SRVY_L2',
						'magData': 'mms1_fgm_b_gsm_srvy_l2',
						'magCoordSys': 'GSM'
						}

	scDic['MMS2'] = {'ephemId':'MMS2_FGM_SRVY_L2',
						'ephemData':'mms2_fgm_r_gsm_srvy_l2',
						'ephemCoordSys':'GSM',
						'denId':'MMS2_FPI_FAST_L2_DIS-MOMS',
						'denData':'mms2_dis_numberdensity_fast',
						'presId': 'MMS2_FPI_FAST_L2_DIS-MOMS',
						'presData': 'mms2_dis_pres_bg_fast',
						'velId': 'MMS2_FPI_FAST_L2_DIS-MOMS',
						'velData':'mms2_dis_bulkv_gse_fast',
						'velCoordSys':'GSE',
						'magId': 'MMS2_FGM_SRVY_L2',
						'magData': 'mms2_fgm_b_gsm_srvy_l2',
						'magCoordSys': 'GSM'
						}
	scDic['MMS3'] = {'ephemId':'MMS3_FGM_SRVY_L2',
						'ephemData':'mms3_fgm_r_gsm_srvy_l2',
						'ephemCoordSys':'GSM',
						'denId':'MMS3_FPI_FAST_L2_DIS-MOMS',
						'denData':'mms3_dis_numberdensity_fast',
						'presId': 'MMS3_FPI_FAST_L2_DIS-MOMS',
						'presData': 'mms3_dis_pres_bg_fast',
						'velId': 'MMS3_FPI_FAST_L2_DIS-MOMS',
						'velData':'mms3_dis_bulkv_gse_fast',
						'velCoordSys':'GSE',
						'magId': 'MMS3_FGM_SRVY_L2',
						'magData': 'mms3_fgm_b_gsm_srvy_l2',
						'magCoordSys': 'GSM'
						}
	scDic['MMS4'] = {'ephemId':'MMS4_FGM_SRVY_L2',
						'ephemData':'mms4_fgm_r_gsm_srvy_l2',
						'ephemCoordSys':'GSM',
						'denId':'MMS4_FPI_FAST_L2_DIS-MOMS',
						'denData':'mms4_dis_numberdensity_fast',
						'presId': 'MMS4_FPI_FAST_L2_DIS-MOMS',
						'presData': 'mms1_dis_pres_bg_fast',
						'velId': 'MMS4_FPI_FAST_L2_DIS-MOMS',
						'velData':'mms4_dis_bulkv_gse_fast',
						'velCoordSys':'GSE',
						'magId': 'MMS1_FGM_SRVY_L2',
						'magData': 'mms4_fgm_b_gsm_srvy_l2',
						'magCoordSys': 'GSM'
						}

	return scDic

def pullVar(cdaObsId,cdaDataId,t0,t1,deltaT):
	cdas = CdasWs()
	status,data =  cdas.get_data(cdaObsId,cdaDataId,t0,t1,
								 binData={
									 'interval': deltaT,
									 'interpolateMissingValues': True,
									 'sigmaMultipler': 4})
	return status,data

def addVar(mydata,scDic,varname,idname,dataname,t0,t1,deltaT):
	#print(scDic,varname,idname,dataname,scDic[idname])
	if scDic[idname] is not None:
		status,data = pullVar(scDic[idname],scDic[dataname],t0,t1,deltaT)
		#print(status)
		if status['http']['status_code'] == 200:
			mydata[varname] = dm.dmarray(data[scDic[dataname]],
													 attrs=data[scDic[dataname]].attrs)
			#mydata.tree(attrs=True)
	else:
		#Mimic the cdasws return code for case when id isn't provided
		status = {'http': {'status_code': 404}}
	return status

def getSatData(scDic,t0,t1,deltaT):
	#First get the empheris data if it doesn't exist return the failed status code and
	#go no further
	status,data = pullVar(scDic['ephemId'],scDic['ephemData'],t0,t1,deltaT)
	if status['http']['status_code'] != 200:
		print('Unable to get data for ', scDic['ephemId'])
		return status,data
	else:
		#data.tree(attrs=True)
		mydata = dm.SpaceData(attrs={'Satellite':data.attrs['Source_name']})
		if 'Epoch_bin' in data.keys():
			#print('Using Epoch_bin')
			mytime = data['Epoch_bin']
		elif 'Epoch' in data.keys():
			#print('Using Epoch')
			mytime = data['Epoch']
		elif ([key for key in data.keys() if key.endswith('_state_epoch')]):
			mytime = data[[key for key in data.keys()
			if key.endswith('_state_epoch')][0]]
		else:
			print('Unable to determine time type')
			status = {'http': {'status_code': 404}}
			return status,data
		mydata['Epoch_bin'] = dm.dmarray(mytime,
										 attrs=mytime.attrs)
		mydata['Ephemeris'] = dm.dmarray(data[scDic['ephemData']],
										 attrs= data[scDic['ephemData']].attrs)
		status1 = addVar(mydata,scDic,'MagneticField','magId','magData',
			t0,t1,deltaT)
		status1 = addVar(mydata,scDic,'Velocity','velId','velData',
			t0,t1,deltaT)
		status1 = addVar(mydata,scDic,'Density','denId','denData',
			t0,t1,deltaT)
		status1 = addVar(mydata,scDic,'Pressure','presId','presData',
			t0,t1,deltaT)
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

def convertGameraVec(x,y,z,ut,fromSys,fromType,toSys,toType):
	invec = Coords(np.column_stack((x,y,z)),fromSys,fromType)
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
	if 'SM' == scDic['ephemCoordSys']:
		smpos = Coords(data['Ephemeris'][:,0:3]*toRe,'SM','car')
		smpos.ticks = Ticktock(data['Epoch_bin'])
	elif 'GSM' == scDic['ephemCoordSys'] :
		scpos = Coords(data['Ephemeris'][:,0:3]*toRe,'GSM','car')
		scpos.ticks = Ticktock(data['Epoch_bin'])
		smpos = scpos.convert('GSE','car')
		scpos = Coords(data['Ephemeris'][:,0:3]*toRe,'GSM','car')
		scpos.ticks = Ticktock(data['Epoch_bin'])
		smpos = scpos.convert('SM','car')
	elif 'GSE'== scDic['ephemCoordSys']:
		scpos = Coords(data['Ephemeris'][:,0:3]*toRe,'GSE','car')
		scpos.ticks = Ticktock(data['Epoch_bin'])
		smpos = scpos.convert('SM','car')
	else:
		print('Coordinate system transformation failed')
		return
	elapsedSecs = (smpos.ticks.getMJD()-mjd0)*86400.0+sec0
	scTrackName = os.path.join(fdir,scId+".sctrack.h5")
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
	ut = MJD2UT(h5file['MJDs'][:])
	
	bx = h5file['Bx']
	by = h5file['By']
	bz = h5file['Bz']

	if scDic['magCoordSys'] is None:
		toCoordSys = 'GSM'
	else:
		toCoordSys = scDic['magCoordSys']
	lfmb_out = convertGameraVec(bx[:],by[:],bz[:],ut,
		'SM','car',toCoordSys,'car')
	data['GAMERA_MagneticField'] = dm.dmarray(lfmb_out.data,
		attrs={'UNITS':bx.attrs['Units'],
		'CATDESC':'Magnetic Field, cartesian'+toCoordSys,
		'FIELDNAM':"Magnetic field",'AXISLABEL':'B'})
	vx = h5file['Vx']
	vy = h5file['Vy']
	vz = h5file['Vz']
	if scDic['velCoordSys'] is None:
		toCoordSys = 'GSM'
	else:
		toCoordSys = scDic['velCoordSys']
	lfmv_out = convertGameraVec(vx[:],vy[:],vz[:],ut,
		'SM','car',scDic['magCoordSys'],'car')
	data['GAMERA_Velocity'] = dm.dmarray(lfmv_out.data,
		attrs={'UNITS':vx.attrs['Units'],
		'CATDESC':'Velocity, cartesian'+toCoordSys,
			   'FIELDNAM':"Velocity",'AXISLABEL':'V'})
	den = h5file['D']
	data['GAMERA_Density'] = dm.dmarray(den[:],
		attrs={'UNITS':den.attrs['Units'],
		'CATDESC':'Density','FIELDNAM':"Density",'AXISLABEL':'n'})
	pres = h5file['P']
	data['GAMERA_Pressure'] = dm.dmarray(pres[:],
		attrs={'UNITS':pres.attrs['Units'],
		'CATDESC':'Pressure','FIELDNAM':"Pressure",'AXISLABEL':'P'})
	temp = h5file['T']
	data['GAMERA_Temperature'] = dm.dmarray(temp[:],
		attrs={'UNITS':pres.attrs['Units'],
		'CATDESC':'Temperature','FIELDNAM':"Temperature",
		'AXISLABEL':'T'})
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
		sctrack = subprocess.run([cmd, xmlFileName], cwd=fdir,
							stdout=subprocess.PIPE, stderr=subprocess.PIPE,
							text=True)
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





