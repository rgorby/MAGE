#!/usr/bin/env python

#standard python
import sys
import os
import datetime
import subprocess
from xml.dom import minidom
import argparse
from argparse import RawTextHelpFormatter


#import spacepy and cdasws
import spacepy
from spacepy.coordinates import Coords
from spacepy.time import Ticktock
import spacepy.datamodel as dm
from cdasws import CdasWs


#import numpy and matplotlib
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#Kaipy and related
from   astropy.time import Time
import h5py
import kaipy.kaiH5 as kaiH5
import kaipy.kaiViz as kv
import kaipy.kaiTools as kaiTools


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
						'magData': 'B_vec_xyz_gse__C1_CP_FGM_SPIN',
						'magCoordSys': 'GSE'
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

def addVar(mydata,scDic,varname,idname,dataname):
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
		print('Unable to get data for ', scId)
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
		mydata['Epoch_bin'] = dm.dmarray(mytime,
										 attrs=mytime.attrs)
		mydata['Ephemeris'] = dm.dmarray(data[scDic['ephemData']],
										 attrs= data[scDic['ephemData']].attrs)
		status1 = addVar(mydata,scDic,'MagneticField','magId','magData')
		status1 = addVar(mydata,scDic,'Velocity','velId','velData')
		status1 = addVar(mydata,scDic,'Density','denId','denData')
		status1 = addVar(mydata,scDic,'Pressure','presId','presData')

	return status,mydata

def labelStr(data, key, vecComp):
	vecLabel=[ 'x', 'y', 'z' ]
	if (vecComp < 3) and (vecComp > -1):
		label=(data['GAMERA_'+key].attrs['AXISLABEL']+
		vecLabel[vecComp]+' ['+
		data['GAMERA_'+key].attrs['UNITS'].decode()+']')
	else:
		label=(data['GAMERA_'+key].attrs['AXISLABEL']+
		' ['+data['GAMERA_'+key].attrs['UNITS'].decode()+']')
	return label

def itemPlot(Ax,data,key,plotNum,numPlots,vecComp=-1):
	#print(key,vecComp)
	if -1 == vecComp:
		Ax.plot(data['Epoch_bin'],data[key][:])
		Ax.plot(data['Epoch_bin'],data['GAMERA_'+key][:])
	else:
		Ax.plot(data['Epoch_bin'],data[key][:,vecComp])
		Ax.plot(data['Epoch_bin'],data['GAMERA_'+key][:,vecComp])
	if (plotNum % 2) == 0:
		left = True
	else:
		left = False
	label = labelStr(data, key,vecComp)
	if plotNum == (numPlots-1):
		kv.SetAxLabs(Ax,'UT',label,doLeft=left)
		kv.SetAxDate(Ax)
	else:
		kv.SetAxLabs(Ax,None,label,doLeft=left)
	return

def compPlot(plotname,scId,data):

	numPlots = 0
	keysToPlot = []
	keys = data.keys()
	#print(keys)
	if 'Density' in keys:
		numPlots = numPlots + 1
		keysToPlot.append('Density')
	if 'Pressue' in keys:
		numPlots = numPlots + 1
		keysToPlot.append('Pressue')
	if 'Temperature' in keys:
		numPlots = numPlots + 1
		keysToPlot.append('Temperature')
	if 'MagneticField' in keys:
		numPlots = numPlots + 3
		keysToPlot.append('MagneticField')
	if 'Velocity' in keys:
		numPlots = numPlots + 3
		keysToPlot.append('Velocity')

	figsize = (10,10)
	fig = plt.figure(figsize=figsize)
	gs = fig.add_gridspec(numPlots,1)
	plotNum = 0
	for key in keysToPlot:
		#print('Plotting',key)
		if 'MagneticField' == key or 'Velocity' == key:
			doVecPlot = True
		else:
			doVecPlot = False
		if 0 == plotNum:
			Ax1 = fig.add_subplot(gs[plotNum,0])
			if doVecPlot:
				itemPlot(Ax1,data,key,plotNum,numPlots,vecComp=0)
				plotNum = plotNum + 1
				Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=1)
				plotNum = plotNum + 1
				Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=2)
				plotNum = plotNum + 1
			else:
				itemPlot(Ax1,data,key,plotNum,numPlots)
				plotNum = plotNum + 1
		else:
			Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
			if doVecPlot:
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=0)
				plotNum = plotNum + 1
				Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=1)
				plotNum = plotNum + 1
				Ax = fig.add_subplot(gs[plotNum,0],sharex=Ax1)
				itemPlot(Ax,data,key,plotNum,numPlots,vecComp=2)
				plotNum = plotNum + 1
			else:
				itemPlot(Ax,data,key,plotNum,numPlots)
				plotNum = plotNum + 1

	Ax1.legend([scId,'GAMERA'],loc='best')
	Ax1.set_title(plotname)
	plt.subplots_adjust(hspace=0)

	kv.savePic(plotname)

def convertGameraVec(x,y,z,ut,fromSys,fromType,toSys,toType):
	invec = Coords(np.column_stack((x,y,z)),fromSys,fromType)
	invec.ticks = Ticktock(ut)
	outvec = invec.convert(toSys,toType)
	return outvec

def extractGAMERA(data,scDic,mjd0,fdir,cmd):
	Re = 6380.0
	toRe = 1.0
	if "km" in data['Ephemeris'].attrs['UNITS']:
		toRe = 1.0/Re
		if 'SM' == scDic['ephemCoordSys']:
			smpos = Coords(data['Ephemeris']*toRe,'SM','car')
			smpos.ticks = Ticktock(data['Epoch_bin'])
		elif 'GSM' == scDic['ephemCoordSys'] :
			scpos = Coords(data['Ephemeris']*toRe,'GSM','car')
			scpos.ticks = Ticktock(data['Epoch_bin'])
			smpos = scpos.convert('SM','car')
		elif 'GSE'== scDic['ephemCoordSys']:
			scpos = Coords(data['Ephemeris']*toRe,'GSE','car')
			scpos.ticks = Ticktock(data['Epoch_bin'])
			smpos = scpos.convert('SM','car')
		else:
			print('Coordinate system transformation failed')
			return
		elapsedSecs = (smpos.ticks.getMJD()-mjd0)*86400.0
		scId = scDic['ephemId']
		fOut = os.path.join(fdir,scId+"sctrack.h5")
		with h5py.File(fOut,'w') as hf:
			hf.create_dataset("T" ,data=elapsedSecs)
			hf.create_dataset("X" ,data=smpos.x)
			hf.create_dataset("Y" ,data=smpos.y)
			hf.create_dataset("Z" ,data=smpos.z)
		chimpxml = kaiTools.genSCXML(scid=scId,h5traj=os.path.basename(fOut))
		xmlfile = os.path.join(fdir,scId+'.xml')
		with open(xmlfile,"w") as f:
			f.write(chimpxml.toprettyxml())

		sctrack = subprocess.run([cmd, xmlfile], cwd=fdir,
								stdout=subprocess.PIPE, stderr=subprocess.PIPE,
								text=True)

		h5name = os.path.join(fdir, scId + '.sc.h5')
		h5file = h5py.File(h5name, 'r')
		ut = kaiTools.MJD2UT(h5file['MJDs'][:])

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

		subprocess.run(['rm',h5name])
		subprocess.run(['rm',xmlfile])
		subprocess.run(['rm',fOut])
		return

if __name__ == '__main__':
	MainS = """Extracts information from satellite trajectory for various
	spacecraft.  Space craft data is pulled from CDAWeb.  Output CDF files
	contain data pulled from CDAWeb along with data extracted from GAMERA.
	Image files of satellite comparisons are also produced.
	"""


	parser = argparse.ArgumentParser(description=MainS,
		formatter_class=RawTextHelpFormatter)
	parser.add_argument('-id',type=str,metavar='runid',default='msphere',
		help='RunID of data (default: %(default)s)')
	parser.add_argument('-path',type=str,metavar='path',default='.',
		help='Path to directory containing REMIX files (default: %(default)s)')
	parser.add_argument('-cmd',type=str,metavar='command',default=None,
		help='Full path to sctrack.x command')

	args = parser.parse_args()

	fdir = args.path
	ftag = args.id
	cmd = args.cmd

	if fdir == '.':
		fdir = os.getcwd()

	if None == cmd:
		my_env = os.environ.copy()
		cmd = os.path.join(os.getenv('KAIJUDIR'),'build','bin','sctrack.x')
	if not (os.path.isfile(cmd) and os.access(cmd, os.X_OK)):
		print(cmd,'either not found or not executable')
		sys.exit()

	scIds = getScIds()

	#print(cmddir)
	#print('Extracting GAMERA data along',scId, 'trajectory')
	#cmd = "/Users/wiltbemj/src/kaiju/build/bin/sctrack.x"
	#Pull the timestep information from the magnetosphere files
	fname = os.path.join(fdir,ftag+'.mix.h5')
	nsteps,sIds=kaiH5.cntSteps(fname)
	gamMJD=kaiH5.getTs(fname,sIds,aID='MJD')
	gamT=kaiH5.getTs(fname,sIds,aID='time')
	gamUT = kaiTools.MJD2UT(gamMJD)
	## Deal with startup by using the first non zero time as the inital
	## MJD
	loc = np.argwhere(gamT > 0.0)[0][0]
	t0 = gamUT[loc]
	t1 = gamUT[-1]
	deltaT = np.round(gamT[loc+1]-gamT[loc])
	mjd0 = gamMJD[loc]

	#scToDo =['CLUSTER1']
	scToDo = scIds
	for scId in scToDo:
		print('Getting spacecraft data for', scId)
		status,data = getSatData(scIds[scId],
			t0.strftime("%Y-%m-%dT%H:%M:%SZ"),
			t1.strftime("%Y-%m-%dT%H:%M:%SZ"),deltaT)

		if status['http']['status_code'] != 200:
			print('No data available for', scId)
		else:
			print('Extracting GAMERA data')
			extractGAMERA(data,scIds[scId],mjd0,fdir,cmd)
			cdfname = os.path.join(fdir, scId + '.comp.cdf')
			if os.path.exists(cdfname):
				print('Deleting %s' % cdfname)
				os.system('rm %s' % cdfname)
			print('Creating CDF file',cdfname,'with',scId,'and GAMERA data')
			dm.toCDF(cdfname,data)
			plotname = os.path.join(fdir,scId+'.png')
			print('Plotting results to',plotname)
			compPlot(plotname,scId,data)

















