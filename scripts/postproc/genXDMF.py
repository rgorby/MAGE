#!/usr/bin/env python
import argparse
import os
import h5py as h5
import kaipy.kaiH5 as kh5
import kaipy.kaixdmf as kxmf
import xml.etree.ElementTree as et
import xml.dom.minidom
import numpy as np

presets = {"gam", "mhdrcm_iono", "mhdrcm_eq", "mhdrcm_bmin", "rcm3D"}

def getDimInfo(h5fname,s0IDstr,preset):

	result = {}

	if preset == "mhdrcm_eq":
		gridVars = ['xMin', 'yMin']
		with h5.File(h5fname, 'r') as h5f:
			gDims = np.asarray(h5f[s0IDstr][gridVars[0]].shape)
		result['gridVars'] = gridVars
		result['gDims'] = gDims
		result['vDims'] = gDims
		result['Nd'] = len(gDims)
		result['geoStr'] = "X_Y"
		result['topoStr'] = "2DSMesh"
		result['doAppendStep'] = True
	elif preset == "mhdrcm_bmin":
		gridVars = ['xMin','yMin','zMin']
		with h5.File(h5fname, 'r') as h5f:
			gDims = np.asarray(h5f[s0IDstr][gridVars[0]].shape)
		#gDims = np.append(gDims, 1)
		result['gridVars'] = gridVars
		result['gDims'] = gDims
		result['vDims'] = gDims
		result['Nd'] = len(gDims)
		result['geoStr'] = "X_Y_Z"
		result['topoStr'] = "3DSMesh"
		result['doAppendStep'] = True
	elif preset=="rcm3D":
		gridVars = ["rcmxmin_kji", "rcmymin_kji", "rcmalamc_kji"]
		with h5.File(h5fname, 'r') as h5f:
			gDims = np.asarray(h5f[s0IDstr][gridVars[0]].shape)
		result['gridVars'] = gridVars
		result['gDims'] = gDims
		result['vDims'] = gDims
		result['Nd'] = len(gDims)
		result['geoStr'] = "X_Y_Z"
		result['topoStr'] = "3DSMesh"
		result['doAppendStep'] = True
	else:  # gam, mhdrcm_iono, etc.
		#Get root-level XY(Z) dimensions
		#First check to see if they exist
		with h5.File(h5fname,'r') as f5:
			if 'X' not in f5.keys():
				print("No X(YZ) in root vars. Maybe try a preset")
				quit()

		gDims = kh5.getDims(h5fname,doFlip=False) #KJI ordering
		Nd = len(gDims)
		if Nd == 2:
			gridVars = ['X', 'Y']
			topoStr = "2DSMesh"
		elif Nd == 3:
			gridVars = ['X', 'Y', 'Z']
			topoStr = "3DSMesh"

		result['gridVars'] = gridVars
		result['gDims'] = gDims
		result['vDims'] = gDims - 1
		result['Nd'] = len(gDims)
		result['geoStr'] = '_'.join(gridVars)
		result['topoStr'] = topoStr
		result['doAppendStep'] = False

	return result


def addRCMVars(Grid, dimInfo, rcmInfo, sID):

	sIDstr = "Step#" + str(sID) 

	mr_vDims = dimInfo['vDims']  # mhdrcm var dims
	mr_vDimStr = ' '.join([str(v) for v in mr_vDims])
	mr_nDims = len(vDims)
	rcmh5fname = rcmInfo['rcmh5fname']
	rcmVars = rcmInfo['rcmVars'] # List of rcm.h5 variables we want in mhdrcm.xmf
	rcmKs = rcmInfo['rcmKs'] # List if rcm.h5 k values for 3d rcm.h5 vars

	rcm5 = h5.File(rcmh5fname,'r')
	if 'Nk' not in rcmInfo.keys():
		rcmInfo['Nj'], rcmInfo['Ni'] = rcm5[sIDstr]['aloct'].shape
		rcmInfo['Nk'] = rcm5[sIDstr]['alamc'].shape[0]
		#print("Adding Nk to rcmInfo")
	Ni = rcmInfo['Ni']
	Nj = rcmInfo['Nj']
	Nk = rcmInfo['Nk']

	
	for vName in rcmVars:
		doHyperslab = False
		r_var = rcm5[sIDstr][vName]
		r_vShape = r_var.shape
		r_vDimStr = " ".join([str(d) for d in r_vShape])
		r_nDims = len(r_vShape)
		dimTrim = 0
	
		if (r_nDims == 2 and mr_vDims[0] < Nj) or (r_nDims == 3 and mr_vDims[1] < Nj):
			doHyperslab = True
			dimTrim = (Nj - mr_vDims[0]) if mr_nDims == 2 else (Nj - mr_vDims[1])

		if r_nDims == 2 and doHyperslab == False:  # Easy add
				#print("Adding " + vName)
				kxmf.AddData(Grid,rcmh5fname, vName,"Cell",mr_vDimStr,sID)
				continue
	
		#Add data as a hyperslab
		if doHyperslab:
			#Do 2D stuff. If 3D needed, will be added in a sec
			dimStr = "3 2"
			startStr = "{} 0".format(dimTrim)
			strideStr = "1 1"
			numStr = "{} {}".format(Nj-dimTrim, Ni)
			text = "{}:/{}/{}".format(rcmh5fname,sIDstr,vName)

			if r_nDims == 2:
				kxmf.addHyperslab(Grid,vName,mr_vDimStr,dimStr,startStr,strideStr,numStr,r_vDimStr,text)
				continue
			elif r_nDims == 3:
				dimStr = "3 3"
				strideStr = str(Nk+1) + " 1 1"
				numStr = "1 {} {}".format(Nj-dimTrim,Ni)
				for k in rcmKs:
					startStr = "{} {} 0".format(k,dimTrim)
					vName_k = vName + "_k{}".format(k)
					kxmf.addHyperslab(Grid,vName_k,mr_vDimStr,dimStr,startStr,strideStr,numStr,r_vDimStr,text)


if __name__ == "__main__":

	outfname = ''

	MainS = """Generates XDMF file from non-MPI HDF5 output
	"""

	parser = argparse.ArgumentParser(description="Generates XDMF file from Gamera HDF5 output")
	parser.add_argument('h5F',type=str,metavar='model.h5',help="Filename of Kaiju HDF5 Output")
	parser.add_argument('-outname',type=str,default=outfname,help="Name of generated XMF file")
	parser.add_argument('-preset', type=str,choices=presets,help="Tell the script what the file is (in case not derivble from name)")
	parser.add_argument('-rcmf',type=str,default="msphere.rcm.h5",help="rcm.h5 file to use with '-rcmv' and '-rcmk' args (default: %(default)s)")
	parser.add_argument('-rcmv',type=str,help="Comma-separated rcm.h5 vars to include in an mhdrcm preset (ex: rcmvm, rcmeeta)")
	parser.add_argument('-rcmk',type=str,help="Comma-separated RCM k values to pull from 3D vars specified with '-rcmv'")
	args = parser.parse_args()

	h5fname = args.h5F
	outfname = args.outname
	preset = args.preset
	rcmh5fname = args.rcmf
	rcmVars = args.rcmv
	rcmKs = args.rcmk

	pre,ext = os.path.splitext(h5fname)
	if outfname is None or outfname == "":
		fOutXML = pre + ".xmf"
	else:
		fOutXML = outfname

	#Scrape necessary data from H5 file
	nSteps,sIDs = kh5.cntSteps(h5fname)
	sIDs = np.sort(sIDs)
	sIDstrs = ['Step#'+str(s) for s in sIDs]
	s0 = sIDs[0]
	s0str = sIDstrs[0]

	#Determine grid and dimensionality
	if preset is None: preset = ""
	dimInfo = getDimInfo(h5fname, s0str, preset)
	gridVars = dimInfo['gridVars']
	gDims = dimInfo['gDims']
	vDims = dimInfo['vDims']
	Nd = dimInfo['Nd']
	geoStr = dimInfo['geoStr']
	topoStr = dimInfo['topoStr']
	doAppendStep = dimInfo['doAppendStep']
	gDimStr = ' '.join([str(v) for v in gDims])
	vDimStr = ' '.join([str(v) for v in vDims])

	doAddRCMVars = False
	#Prep to include some rcmh5 vars in mhdrcm.xmf file
	if 'mhdrcm' in preset and rcmVars is not None:
		doAddRCMVars = True
		rcmVars = rcmVars.split(',')
		rcmKs = [int(k) for k in rcmKs.split(',')] if rcmKs is not None else []
		rcmInfo = {}
		rcmInfo['rcmh5fname'] = rcmh5fname
		rcmInfo['rcmVars'] = rcmVars
		rcmInfo['rcmKs'] = rcmKs

	#Get variable information
	#T = kh5.getTs(h5fname,sIDs)
	T = np.array([])
	with h5.File(h5fname,'r') as f5:
		#T = np.array([f5[k].attrs['time'] for k in f5.keys() if "Step" in k])
		for sIDstr in sIDstrs:
			if 'time' in f5[sIDstr].attrs.keys():
				T = np.append(T, f5[sIDstr].attrs['time'])
			else:
				T = np.append(T, int(sIDstr.split("#")[1]))
		#steps = np.array([k for k in f5.keys() if "Step" in k])
	Nt = len(T)
	vIds ,vLocs  = kxmf.getVars(h5fname,s0str,gDims)
	rvIds,rvLocs = kxmf.getRootVars(h5fname,gDims)

	Nv = len(vIds)
	Nrv = len(rvIds)


	print("Generating XDMF from %s"%(h5fname))
	print("Writing to %s"%(fOutXML))
	print("\t%d Time Slices / %d Variables"%(nSteps,len(vIds)))
	print("\tGrid: %s"%str(gDims))
	print("\tSlices: %d -> %d"%(sIDs.min(),sIDs.max()))
	print("\tTime: %3.3f -> %3.3f"%(T.min(),T.max()))
	

	#Construct XDMF XML file
	#-----------------------
	Xdmf = et.Element("Xdmf")
	Xdmf.set("Version","2.0")
	
	Dom = et.SubElement(Xdmf,"Domain")
	
	TGrid = et.SubElement(Dom,"Grid")
	TGrid.set("Name","tMesh")
	TGrid.set("GridType","Collection")
	TGrid.set("CollectionType","Temporal")

	#Loop over time slices
	for n in range(Nt):
		nStp = sIDs[n]
		Grid = et.SubElement(TGrid,"Grid")
		mStr = "gMesh"#+str(nStp)
		Grid.set("Name",mStr)
		Grid.set("GridType","Uniform")

		Topo = et.SubElement(Grid,"Topology")
		Topo.set("TopologyType",topoStr)
		Topo.set("NumberOfElements",gDimStr)
		Geom = et.SubElement(Grid,"Geometry")
		Geom.set("GeometryType",geoStr)
		
		#Add grid info to each step
		if doAppendStep:
			stepStr = sIDstrs[n]
			sgVars = [os.path.join(stepStr, v) for v in gridVars]
			kxmf.AddGrid(h5fname,Geom,gDimStr,sgVars)
		else:
			kxmf.AddGrid(h5fname,Geom,gDimStr,gridVars)

		Time = et.SubElement(Grid,"Time")
		Time.set("Value","%f"%T[n])

		if preset=="rcm3D":
			with h5.File(h5fname,'r') as f5:
				other  = et.SubElement(Grid, "dtCpl")
				other.set("Value","%f"%f5[sIDstrs[n]].attrs['dtCpl'])

		#--------------------------------
		#Step variables
		for v in range(Nv):
			kxmf.AddData(Grid,h5fname,vIds[v],vLocs[v],vDimStr,nStp)
		#--------------------------------
		#Base grid variables
		for v in range(Nrv):
			kxmf.AddData(Grid,h5fname,rvIds[v],rvLocs[v],vDimStr)

		if doAddRCMVars:
			addRCMVars(Grid, dimInfo, rcmInfo, sIDs[n])

		#--------------------------------
		#Add some extra aliases
		if preset=="gam":
			kxmf.AddVectors(Grid,h5fname,vIds,cDims,vDims,Nd,nStp)

	#Finished creating XML tree, now write
	xmlStr = xml.dom.minidom.parseString(et.tostring(Xdmf)).toprettyxml(indent="    ")
	print("Saving as {}".format(fOutXML))
	with open(fOutXML,"w") as f:
		f.write(xmlStr)
		