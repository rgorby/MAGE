#!/usr/bin/env python
#Make XMF files from an MPI-decomposed gamera run

import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import kaipy.gamera.gampp as gampp
import kaipy.kaixdmf as kxmf
import xml.etree.ElementTree as et
import xml.dom.minidom
import kaipy.kaiH5 as kh5

import os

if __name__ == "__main__":
	#Defaults
	fdir = os.getcwd()
	ftag = "msphere"
	outid = "sim"
	sStride = 1
	MainS = """Creates series of XMF files from MPI-decomposed Gamera run
	"""

	parser = argparse.ArgumentParser(description=MainS, formatter_class=RawTextHelpFormatter)
	parser.add_argument('-d',type=str,metavar="directory",default=fdir,help="Directory to read from (default: %(default)s)")
	parser.add_argument('-id',type=str,metavar="runid",default=ftag,help="RunID of data (default: %(default)s)")
	parser.add_argument('-outid',type=str,metavar="outid",default=outid,help="RunID of output XMF files (default: %(default)s)")
	parser.add_argument('-sS',type=int,metavar="stride",default=sStride,help="Output cadence (default: %(default)s)")


	#Finalize parsing
	args = parser.parse_args()
	fdir = args.d
	ftag = args.id
	outid = args.outid
	sStride = args.sS

	#---------------------
	#Init data
	
	gamData = gampp.GameraPipe(fdir,ftag)
	
	#---------------------
	#Do work
	Ri = gamData.Ri
	Rj = gamData.Rj
	Rk = gamData.Rk

	Ni = gamData.dNi
	Nj = gamData.dNj
	Nk = gamData.dNk

	iDims = "%d %d %d"%(Nk+1,Nj+1,Ni+1)
	cDims = "%d %d %d"%(Nk+0,Nj+0,Ni+0)
	iDimA = "%d %d %d"%(Nk*Rk+1,Nj*Rj+1,Ni*Ri+1)
	tOut = 0

	topoStr = "3DSMesh"
	geoStr = "X_Y_Z"	

	#for n in range(gamData.s0,gamData.sFin+1,sStride):
	for n in gamData.sids:
		if (n-gamData.s0)%sStride != 0: continue
		nslc = n-gamData.s0
		#print(n,gamData.T[nslc])

		#Create XMF tree
		Xdmf = et.Element("Xdmf")
		Xdmf.set("Version","2.0")
		
		Dom = et.SubElement(Xdmf,"Domain")


		#Spatial collection
		gCol = et.SubElement(Dom,"Grid")
		gCol.set("Name","gMesh")
		gCol.set("GridType","Collection")
		gCol.set("CollectionType","Spatial")

		Time = et.SubElement(gCol,"Time")
		#Time.set("Value","%s"%(str(gamData.T[np.where(gamData.sids == n)][0])))
		Time.set("Value","%s"%(gamData.T[nslc]))
		
		Cyc = et.SubElement(gCol,"Cycle")
		Cyc.set("Value","%d"%(n))

		#Now loop over MPI decomposition
		for i in range(Ri):
			for j in range(Rj):
				for k in range(Rk):
					nMPI = j + i*Rj + k*Ri*Rj
					h5F = kh5.genName(ftag,i,j,k,Ri,Rj,Rk)
					h5F = os.path.join(fdir, h5F)
					
					#Get variable info
					gDims = np.array([gamData.dNk+1,gamData.dNj+1,gamData.dNi+1])
					vDims = np.array([gamData.dNk,gamData.dNj,gamData.dNi])
					vDimStr = ' '.join([str(v) for v in vDims])
					if nMPI == 0 and tOut == 0:  # Only do this the first time
						vIds ,vLocs  = kxmf.getVars(h5F,'Step#'+str(n),gDims)
						rvIds,rvLocs = kxmf.getRootVars(h5F,gDims)
						Nv = len(vIds)
						Nrv = len(rvIds)

					#Create new subgrid
					gName = "gMesh%d"%(nMPI)
					Grid = et.SubElement(gCol,"Grid")
					Grid.set("GridType","Uniform")
					Grid.set("Name",gName)

					# Time = et.SubElement(Grid,"Time")
					# Time.set("TimeType","Single")
					# Time.set("Value","%s"%(str(gamData.T[nslc])))

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

					zC = et.SubElement(Geom,"DataItem")
					zC.set("Dimensions",iDims)
					zC.set("NumberType","Float")
					zC.set("Precision","4")
					zC.set("Format","HDF")
					zC.text = "%s:/Z"%(h5F)

					#Create variables
					"""
					for v in range(Nv):
						vID = vIDs[v]
						vAtt = et.SubElement(Grid,"Attribute")
						vAtt.set("Name",vID)
						vAtt.set("AttributeType","Scalar")
						vAtt.set("Center","Cell")
						aDI = et.SubElement(vAtt,"DataItem")
						aDI.set("Dimensions",cDims)
						aDI.set("NumberType","Float")
						aDI.set("Precision","4")
						aDI.set("Format","HDF")
						aDI.text = "%s:/Step#%d/%s"%(h5F,n,vID)
					"""
					for v in range(Nv):
						kxmf.AddData(Grid,h5F,vIds[v],vLocs[v],vDimStr,n)
					for rv in range(Nrv):
						kxmf.AddData(Grid,h5F,rvIds[rv],rvLocs[rv],vDimStr)

		#Write output
		fOut = "%s/%s.%06d.xmf"%(fdir,outid,tOut)
		print("Writing %s"%(fOut))
		#xTree = et.ElementTree(Xdmf)
		#xTree.write(fOut,pretty_print=True,xml_declaration=True,encoding='UTF-8')
		xmlStr = xml.dom.minidom.parseString(et.tostring(Xdmf)).toprettyxml(indent="    ")
		with open(fOut,"w") as f:
			f.write(xmlStr)

		#Prep for next step
		tOut = tOut+1
