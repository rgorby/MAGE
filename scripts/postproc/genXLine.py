#!/usr/bin/env python
import argparse
import os
import h5py
import kaipy.kaiH5 as kh5
#import lxml.etree as et
import xml.etree.ElementTree as et
import xml.dom.minidom
import numpy as np
import kaipy.kaiH5 as kh5

def cntX(fname,gID=None,StrX="/Step#"):
	with h5py.File(fname,'r') as hf:
		if (gID is not None):
			grps = hf[gID].values()
		else:
			grps = hf.values()
		grpNames = [str(grp.name) for grp in grps]
		#Steps = [stp if "/Step#" in stp for stp in grpNames]
		Steps = [stp for stp in grpNames if StrX in stp]
		nSteps = len(Steps)

		sIds = np.array([str.split(s,"#")[-1] for s in Steps],dtype=np.int)
		return nSteps,sIds

def getVars(fname,smin):
	#Get variable names from Step#0/Line#0
	with h5py.File(fname,'r') as hf:
		gID = "/Step#%d/Line#0"%(smin)
		vIDs = []
		for k in hf[gID].keys():
			vIDs.append(str(k))
	#Remove coordinate vars
	xyzS = ["xyz","LCon"]
	for s in xyzS:
		if s in vIDs:
			vIDs.remove(s)
	Nv = len(vIDs)
	return Nv,vIDs

def getAtts(fIn,n,m):
	#Get attribute names from Step#0/Line#0
	with h5py.File(fIn,'r') as hf:
		gId = "Step#%d"%(n)
		lId = "Line#%d"%(m)
		Atts = hf[gId][lId].attrs.keys()

		aNull = ["Np","n0"]
		aIDs = [x for x in Atts if x not in aNull]
		aVs = []
		Na = len(aIDs)
		for n in range(Na):
			aVs.append(np.float(hf[gId][lId].attrs[aIDs[n]]))

		aVs.append(np.float(m))
		aIDs.append("ID")
		
	return aIDs,aVs


def getNum(fIn,n,m):
	with h5py.File(fIn,'r') as hf:
		gId = "Step#%d"%(n)
		lId = "Line#%d"%(m)
		Np = hf[gId][lId].attrs["Np"]
	return Np
if __name__ == "__main__":
    #Set defaults
	parser = argparse.ArgumentParser(description="Generates XDMF file from CHIMP tracer HDF5 output")
	parser.add_argument('h5F',nargs=1,type=str,metavar='tracer.h5',help="Filename of CHIMP tracer HDF5 Output")
	parser.add_argument('-noatts', action='store_true', default=False,help="Don't add XDMF scalars (default: %(default)s)")
	#Finished getting arguments, parse and move on
	args = parser.parse_args()

	fIn = args.h5F[0]
	doAtts =  not args.noatts
	
	#Create XML filename
	pre,ext = os.path.splitext(fIn)
	fOutXML = pre + ".xmf"

	print("Reading from %s"%(fIn))
	kh5.CheckOrDie(fIn)
	
	#Count steps and lines
	
	Nstp,sIds = cntX(fIn)
	s0 = sIds.min()
	if (s0 == 0):
		s0 = s0+1
		Nstp = Nstp-1
	gID = "Step#%d"%(s0)
	Nl,lIds = cntX(fIn,gID=gID,StrX="Line#")
	Nv,vIds = getVars(fIn,s0)

	
	print("\tFound %d steps"%(Nstp))
	print("\tFound %d lines/step"%(Nl))
	print("\tFound %d vars/line"%(Nv))
	#Get times
	T = kh5.getTs(fIn,sIds)

	#Now build XDMF file
	#-----------------------
	Xdmf = et.Element("Xdmf")
	Xdmf.set("Version","2.0")

	Dom = et.SubElement(Xdmf,"Domain")

	TGrid = et.SubElement(Dom,"Grid")
	TGrid.set("Name","tlMesh")
	TGrid.set("GridType","Collection")
	TGrid.set("CollectionType","Temporal")

	#Loop over time slices
	for n in range(Nstp):
		lGrid = et.SubElement(TGrid,"Grid")
		lGrid.set("Name","tLines")
		lGrid.set("GridType","Collection")
		#Add time
		tLab = et.SubElement(lGrid,"Time")
		tLab.set("Value","%f"%T[n])
		nStp = n + s0
		#Loop over individual lines
		for m in range(Nl):
			#Get number of points for this step/line
			Np = getNum(fIn,nStp,m)

			#Create main grid structure
			l0G = et.SubElement(lGrid,"Grid")
			l0G.set("GridType","Uniform")
			l0G.set("Name","Line#%d"%(m))

			#Add topology/connectivity
			Topo = et.SubElement(l0G,"Topology")
			Topo.set("TopologyType","Polyline")
			Topo.set("NumberOfElements","1")
			Topo.set("NodesPerElement",str(Np-1))

			Geom = et.SubElement(l0G,"Geometry")
			Geom.set("GeometryType","XYZ")
			xC = et.SubElement(Geom,"DataItem")
			xC.set("Dimensions","%d 3"%(Np))
			xC.set("NumberType","Float")
			xC.set("Precision","4")
			xC.set("Format","HDF")
			xC.text = "%s:/Step#%d/Line#%d/xyz"%(fIn,nStp,m)


			#Now loop over variables
			for v in range(Nv):
				vAtt = et.SubElement(l0G,"Attribute")
				vAtt.set("Name",vIds[v])
				vAtt.set("AttributeType","Scalar")
				vAtt.set("Center","Node")
				vDI = et.SubElement(vAtt,"DataItem")
				vDI.set("Dimensions",str(Np))
				vDI.set("NumberType","Float")
				vDI.set("Precision","4")
				vDI.set("Format","HDF")
				vDI.text = "%s:/Step#%d/Line#%d/%s"%(fIn,nStp,m,vIds[v])
			if (doAtts):
				#Add scalar attributes in lazy XDMF way
				aIDs,aVs = getAtts(fIn,nStp,m)
				Na = len(aIDs)
				for a in range(Na):
					#Main variable
					vAtt = et.SubElement(l0G,"Attribute")
					vAtt.set("Name",aIDs[a])
					vAtt.set("AttributeType","Scalar")
					vAtt.set("Center","Node")
					#Function setup
					vDI = et.SubElement(vAtt,"DataItem")
					vDI.set("Dimensions",str(Np))
					vDI.set("Function","%e*$0/$0"%(np.float(aVs[a])))
					vDI.set("ItemType","Function")
					#Argument
					vNull = et.SubElement(vDI,"DataItem")
					vNull.set("Dimensions",str(Np))
					vNull.set("NumberType","Float")
					vNull.set("Precision","4")
					vNull.set("Format","HDF")
					vNull.text = "%s:/Step#%d/Line#%d/%s"%(fIn,nStp,m,vIds[v])

	#Finished creating XML, now write
	xmlStr = xml.dom.minidom.parseString(et.tostring(Xdmf)).toprettyxml(indent="    ")
	with open(fOutXML,"w") as f:
		f.write(xmlStr)
