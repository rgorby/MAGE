#!/usr/bin/env python
import argparse
import os
import h5py
import xml.etree.ElementTree as et
import xml.dom.minidom
import numpy as np

def cntSteps(fname):
	with h5py.File(fname,'r') as hf:
			grps = hf.values()
			grpNames = [str(grp.name) for grp in grps]
			#Steps = [stp if "/Step#" in stp for stp in grpNames]
			Steps = [stp for stp in grpNames if "/Step#" in stp]
			nSteps = len(Steps)
			
			sIds = np.array([str.split(s,"#")[-1] for s in Steps],dtype=np.int)
			return nSteps,sIds

def getTs(fname,sIds):
	tDef = 0.0 #Default time
	Nt = len(sIds)
	T = np.zeros(Nt)
	s0 = sIds.min()
	sF = sIds.max()
	with h5py.File(fname,'r') as hf:
		for n in range(s0,sF+1):
			gId = "/Step#%d"%(n)
			T[n-s0] = hf[gId].attrs.get("time",tDef)
	return T

#Get shape/dimension of grid
def getDims(fname,s0):
	with h5py.File(fname,'r') as hf:
		Dims = hf["/"]["X"].shape
	return np.array(Dims,dtype=np.int)

#Get root variables
def getRootVars(fname):
	with h5py.File(fname,'r') as hf:
		vIds = []
		for k in hf.keys():
			#Don't include stuff that starts with step
			if "Step" not in str(k):
				vIds.append(str(k))
	#Remove coordinates from list of root variables
	xyzS = ["X","Y","Z"]
	for s in xyzS:
		if s in vIds:
			vIds.remove(s)
	return vIds

#Get variables in initial Step
def getVars(fname,s0):
	Dims = getDims(fname,s0)
	Ni = Dims[0]-1
	Nj = Dims[1]-1

	with h5py.File(fname,'r') as hf:
		gId = "/Step#%d"%(s0)
		stp0 = hf[gId]
		vIds = []
		vLocs = []
		for k in stp0.keys():
			Nvi = stp0[k].shape[0]
			#Add variable
			vIds.append(str(k))
			#Decide if variable is cell-centered or node centered
			#Just using first dimension
			if (Nvi == Ni):
				vLocs.append("Cell")
			elif (Nvi == (Ni+1)):
				vLocs.append("Node")
			else:
				vLocs.append("Other")
			
	return vIds,vLocs


#Add data item to passed element
def AddDI(elt,h5F,nStp,cDims,vId):
	aDI = et.SubElement(elt,"DataItem")
	aDI.set("Dimensions",cDims)
	aDI.set("NumberType","Float")
	aDI.set("Precision","4")
	aDI.set("Format","HDF")
	if (nStp>=0):
		aDI.text = "%s:/Step#%d/%s"%(h5F,nStp,vId)
	else:
		aDI.text ="%s:/%s"%(h5F,vId)
#Add data item reference to passed element
def AddDIRef(elt,nStp,cDims,vId):
	aDI = et.SubElement(elt,"DataItem")
	aDI.set("Dimensions",cDims)
	aDI.set("Precision","4")
	aDI.set("Format","HDF")
	#First grid is for tMesh (temporal), second is step #
	xPath = '/Xdmf/Domain/Grid[1]/Grid[%d]/Attribute[@Name=\"%s\"]/DataItem[1]'%(nStp+1,vId)
	aDI.set("Reference",xPath)

if __name__ == "__main__":
	#Set defaults
	sDim = ["x","y","z"]

	parser = argparse.ArgumentParser(description="Generates XDMF file from Gamera HDF5 output")

	parser.add_argument('h5F',nargs='+',metavar='Gamera.h5',help="Filename of Gamera HDF5 Output")

	#Finished getting arguments, parse and move on
	args = parser.parse_args()
	
	for idx, h5F in enumerate(args.h5F):
		pre,ext = os.path.splitext(h5F)
		fOutXML = pre + ".xmf"
	
		#Scrape necessary data from H5 file
		nSteps,sIds = cntSteps(h5F)
		s0 = sIds.min()
	
		T = getTs(h5F,sIds)
		Nt = len(T)
		vIds,vLocs = getVars(h5F,s0)
		rvIds = getRootVars(h5F)

		Nv = len(vIds)
		Nrv = len(rvIds)

		Dims = getDims(h5F,s0)
		Dims = Dims-1 #Correct for interface vs. cell-centered
		Nd = len(Dims)


		print("Generating XDMF from %s"%(h5F))
		print("Writing to %s"%(fOutXML))
		print("\t%d Time Slices / %d Variables"%(nSteps,len(vIds)))
		print("\tGrid: %s"%str(Dims))
		print("\tSlices: %d -> %d"%(sIds.min(),sIds.max()))
		print("\tTime: %3.3f -> %3.3f"%(T.min(),T.max()))
		
		#Prepare for XDMF file
		Nx = Dims[0]
		Ny = Dims[1]
		Nz = 0
	
		if (Nd == 3):
			Nz = Dims[2]
			iDims = "%d %d %d"%(Nx+1,Ny+1,Nz+1)
			cDims = "%d %d %d"%(Nx,Ny,Nz)
			vDims = "%d %d %d %d"%(Nx,Ny,Nz,Nd)
			topoStr = "3DSMesh"
			geoStr = "X_Y_Z"
		else:
			iDims = "%d %d"%(Nx+1,Ny+1)
			cDims = "%d %d"%(Nx,Ny)
			vDims = "%d %d %d"%(Nx,Ny,Nd)
			topoStr = "2DSMesh"
			geoStr = "X_Y"

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
			nStp = n + s0
			Grid = et.SubElement(TGrid,"Grid")
			mStr = "gMesh"#+str(nStp)
			Grid.set("Name",mStr)
			Grid.set("GridType","Uniform")
	
			Topo = et.SubElement(Grid,"Topology")
			Topo.set("TopologyType",topoStr)
			Topo.set("NumberOfElements",iDims)
			Geom = et.SubElement(Grid,"Geometry")
			Geom.set("GeometryType",geoStr)
			
			xC = et.SubElement(Geom,"DataItem")
			xC.set("Dimensions",iDims)
			xC.set("NumberType","Float")
			xC.set("Precision","4")
			xC.set("Format","HDF")
			xC.text = "%s:/X"%(h5F)
			
			yC = et.SubElement(Geom,"DataItem")
			yC.set("Dimensions",iDims)
			yC.set("NumberType","Float")
			yC.set("Precision","4")
			yC.set("Format","HDF")
			yC.text = "%s:/Y"%(h5F)
	
			if (Nd == 3):
				zC = et.SubElement(Geom,"DataItem")
				zC.set("Dimensions",iDims)
				zC.set("NumberType","Float")
				zC.set("Precision","4")
				zC.set("Format","HDF")
				zC.text = "%s:/Z"%(h5F)
	
			Time = et.SubElement(Grid,"Time")
			Time.set("Value","%f"%T[n])    
	
			#--------------------------------
			#Step variables
			for v in range(Nv):
				vAtt = et.SubElement(Grid,"Attribute")
				vAtt.set("Name",vIds[v])
				vAtt.set("AttributeType","Scalar")
				vAtt.set("Center",vLocs[v])
				
				AddDI(vAtt,h5F,nStp,cDims,vIds[v])

			#--------------------------------
			#Add base grid stuff

			for v in range(Nrv):
				vAtt = et.SubElement(Grid,"Attribute")
				vAtt.set("Name",rvIds[v])
				vAtt.set("AttributeType","Scalar")
				vAtt.set("Center","Cell")
				
				AddDI(vAtt,h5F,-1,cDims,rvIds[v])

			#--------------------------------
			#Add some extra aliases
			
			#Velocity (2D)
			if ( (Nd == 2) and ("Vx" in vIds) and ("Vy" in vIds) ):
				vAtt = et.SubElement(Grid,"Attribute")
				vAtt.set("Name","VecV")
				vAtt.set("AttributeType","Vector")
				vAtt.set("Center","Cell")
				fDI = et.SubElement(vAtt,"DataItem")
				fDI.set("ItemType","Function")
				fDI.set("Dimensions",vDims)
				fDI.set("Function","JOIN($0 , $1)")
				AddDI(fDI,h5F,nStp,cDims,"Vx")
				AddDI(fDI,h5F,nStp,cDims,"Vy")

			#Velocity (3D)
			if ( (Nd == 3) and ("Vx" in vIds) and ("Vy" in vIds) and ("Vz" in vIds)):
				vAtt = et.SubElement(Grid,"Attribute")
				vAtt.set("Name","VecV")
				vAtt.set("AttributeType","Vector")
				vAtt.set("Center","Cell")
				fDI = et.SubElement(vAtt,"DataItem")
				fDI.set("ItemType","Function")
				fDI.set("Dimensions",vDims)
				fDI.set("Function","JOIN($0 , $1 , $2)")
				AddDI(fDI,h5F,nStp,cDims,"Vx")
				AddDI(fDI,h5F,nStp,cDims,"Vy")
				AddDI(fDI,h5F,nStp,cDims,"Vz")

			#Magnetic field (2D)
			if ( (Nd == 2) and ("Bx" in vIds) and ("By" in vIds) ):
				vAtt = et.SubElement(Grid,"Attribute")
				vAtt.set("Name","VecB")
				vAtt.set("AttributeType","Vector")
				vAtt.set("Center","Cell")
				fDI = et.SubElement(vAtt,"DataItem")
				fDI.set("ItemType","Function")
				fDI.set("Dimensions",vDims)
				fDI.set("Function","JOIN($0 , $1)")
				AddDI(fDI,h5F,nStp,cDims,"Bx")
				AddDI(fDI,h5F,nStp,cDims,"By")

			if ( (Nd == 3) and ("Bx" in vIds) and ("By" in vIds) and ("Bz" in vIds)):
				vAtt = et.SubElement(Grid,"Attribute")
				vAtt.set("Name","VecB")
				vAtt.set("AttributeType","Vector")
				vAtt.set("Center","Cell")
				fDI = et.SubElement(vAtt,"DataItem")
				fDI.set("ItemType","Function")
				fDI.set("Dimensions",vDims)
				fDI.set("Function","JOIN($0 , $1 , $2)")
				AddDI(fDI,h5F,nStp,cDims,"Bx")
				AddDI(fDI,h5F,nStp,cDims,"By")
				AddDI(fDI,h5F,nStp,cDims,"Bz")

		#Finished creating XML tree, now write
		xmlStr = xml.dom.minidom.parseString(et.tostring(Xdmf)).toprettyxml(indent="    ")
		with open(fOutXML,"w") as f:
			f.write(xmlStr)
			
		#Below uses lxml instead of xml python module
		#xTree = et.ElementTree(Xdmf)
		#xTree.write(fOutXML,pretty_print=True)
