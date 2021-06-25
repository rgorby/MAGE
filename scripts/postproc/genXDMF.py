#!/usr/bin/env python
import argparse
import os
import h5py
import kaipy.kaiH5 as kh5
import kaipy.kaixdmf as kxmf
import xml.etree.ElementTree as et
import xml.dom.minidom
import numpy as np

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Generates XDMF file from Gamera HDF5 output")
	parser.add_argument('h5F',nargs='+',metavar='Gamera.h5',help="Filename of Gamera HDF5 Output")
	#Finished getting arguments, parse and move on
	args = parser.parse_args()
	
	for idx, h5F in enumerate(args.h5F):
		pre,ext = os.path.splitext(h5F)
		fOutXML = pre + ".xmf"
	
		#Scrape necessary data from H5 file
		nSteps,sIds = kh5.cntSteps(h5F)
		s0 = sIds.min()
	
		T = kh5.getTs(h5F,sIds)
		Nt = len(T)
		vIds ,vLocs  = kxmf.getVars(h5F,s0)
		rvIds,rvLocs = kxmf.getRootVars(h5F)

		Nv = len(vIds)
		Nrv = len(rvIds)

		Dims = kh5.getDims(h5F,doFlip=False) #KJI ordering
		Dims = Dims-1 #Correct for interface vs. cell-centered
		Nd = len(Dims)

		print("Generating XDMF from %s"%(h5F))
		print("Writing to %s"%(fOutXML))
		print("\t%d Time Slices / %d Variables"%(nSteps,len(vIds)))
		print("\tGrid: %s"%str(Dims))
		print("\tSlices: %d -> %d"%(sIds.min(),sIds.max()))
		print("\tTime: %3.3f -> %3.3f"%(T.min(),T.max()))
		
		#Prepare for XDMF file
		Nx1 = Dims[0]
		Nx2 = Dims[1]
		Nx3 = 0
	
		if (Nd == 3):
			Nx3 = Dims[2]
			iDims = "%d %d %d"%(Nx1+1,Nx2+1,Nx3+1)
			cDims = "%d %d %d"%(Nx1,Nx2,Nx3)
			vDims = "%d %d %d %d"%(Nx1,Nx2,Nx3,Nd)
			topoStr = "3DSMesh"
			geoStr = "X_Y_Z"
		else:
			iDims = "%d %d"%(Nx1+1,Nx2+1)
			cDims = "%d %d"%(Nx1,Nx2)
			vDims = "%d %d %d"%(Nx1,Nx2,Nd)
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
			
			#Add grid info to each step
			kxmf.AddGrid(h5F,Geom,iDims,Nd)
	
			Time = et.SubElement(Grid,"Time")
			Time.set("Value","%f"%T[n])    
	
			#--------------------------------
			#Step variables
			for v in range(Nv):
				kxmf.AddData(Grid,h5F,vIds[v],vLocs[v],cDims,nStp)
			#--------------------------------
			#Base grid variables
			for v in range(Nrv):
				kxmf.AddData(Grid,h5F,rvIds[v],rvLocs[v],cDims)

			#--------------------------------
			#Add some extra aliases
			kxmf.AddVectors(Grid,h5F,vIds,cDims,vDims,Nd,nStp)

		#Finished creating XML tree, now write
		xmlStr = xml.dom.minidom.parseString(et.tostring(Xdmf)).toprettyxml(indent="    ")
		with open(fOutXML,"w") as f:
			f.write(xmlStr)
			