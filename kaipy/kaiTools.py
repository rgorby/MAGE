import numpy as np
import datetime
import os
import glob
from astropy.time import Time
from xml.dom import minidom

def MJD2UT(mjd):
	astroT = Time(mjd,format='mjd').iso
	utall = []
	for ut in astroT:
		utall.append(datetime.datetime.strptime(ut,'%Y-%m-%d %H:%M:%S.%f'))

	return utall

def genSCXML(fdir,ftag,
	scid="sctrack_A",h5traj="sctrack_A.h5"):

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
	trajChild.setAttribute("doSmooth","T")
	xml.appendChild(trajChild)

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
						'magData': 'B_vec_xyz_gse__C1_CP_FGM_SPIN',
						'magCoordSys': 'GSE'
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
	scDic['THEMISA'] = {'ephemId':'THA_L1_STATE',
						'ephemData':'tha_pos_gsmV',
						'ephemCoordSys':'GSM',
						'denId':'THA_L2_MOM',
						'denData':'tha_peim_density',
						'presId': 'THA_L2_MOM',
						'presData': 'tha_peim_ptot',
						'velId': 'THA_L2_MOM',
						'velData':'tha_peim_velocity_gsm',
						'velCoordSys':'GSM',
						'magId': 'THA_L2_MOM',
						'magData': 'tha_peim_mag',
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



