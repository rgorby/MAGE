#Various routines to help scripts that create XMF files from H5 data
import numpy as np
import h5py
import kaipy.kaiH5 as kh5
import xml.etree.ElementTree as et

#Add grid info to step
#Geom is topology subelement, iDims is grid size string
def AddGrid(fname,Geom,iDims,Nd):
	xC = et.SubElement(Geom,"DataItem")
	xC.set("Dimensions",iDims)
	xC.set("NumberType","Float")
	xC.set("Precision","4")
	xC.set("Format","HDF")
	xC.text = "%s:/X"%(fname)

	yC = et.SubElement(Geom,"DataItem")
	yC.set("Dimensions",iDims)
	yC.set("NumberType","Float")
	yC.set("Precision","4")
	yC.set("Format","HDF")
	yC.text = "%s:/Y"%(fname)

	if (Nd == 3):
		zC = et.SubElement(Geom,"DataItem")
		zC.set("Dimensions",iDims)
		zC.set("NumberType","Float")
		zC.set("Precision","4")
		zC.set("Format","HDF")
		zC.text = "%s:/Z"%(fname)

#Add data to slice
def AddData(Grid,fname,vID,vLoc,xDims,s0=None):
	if (vLoc != 'Other'):
		#Add attribute
		vAtt = et.SubElement(Grid,"Attribute")
		vAtt.set("Name",vID)			
		vAtt.set("AttributeType","Scalar")			
		vAtt.set("Center",vLoc)			
		#Add data item
		aDI = et.SubElement(vAtt,"DataItem")
		aDI.set("Dimensions",xDims)
		aDI.set("NumberType","Float")
		aDI.set("Precision","4")
		aDI.set("Format","HDF")
		if (s0 is None):
			aDI.text ="%s:/%s"%(fname,vID)
		else:
			aDI.text = "%s:/Step#%d/%s"%(fname,s0,vID)

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

#Get root variables
def getRootVars(fname):
	Dims = kh5.getDims(fname,doFlip=False) #Do kji ordering
	with h5py.File(fname,'r') as hf:
		vIds = []
		vLocs = []
		for k in hf.keys():
			#Don't include stuff that starts with step or X,Y,Z
			vID = str(k)
			doV = True
			if ("Step" in vID):
				doV = False
			if ((vID == "X") or (vID=="Y") or (vID=="Z")):
				doV = False
			if (doV):
				Nv = hf[k].shape
				vLoc = getLoc(Dims,Nv)
				if (vLoc == "Cell"):
					vIds.append(vID)
					vLocs.append(vLoc)
				else:
					print("Excluding %s"%(vID))

	return vIds,vLocs

#Get variables in initial Step
def getVars(fname,s0):
	Dims = kh5.getDims(fname,doFlip=False) #Do kji ordering

	with h5py.File(fname,'r') as hf:
		gId = "/Step#%d"%(s0)
		stp0 = hf[gId]
		vIds = []
		vLocs = []
		for k in stp0.keys():
			vID = str(k)
			Nv = stp0[k].shape
			vLoc = getLoc(Dims,Nv)
			if (vLoc == "Cell"):
				vIds.append(vID)
				vLocs.append(vLoc)
			else:
				print("Excluding %s"%(vID))
	return vIds,vLocs

def AddVectors(Grid,fname,vIds,cDims,vDims,Nd,nStp):
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
		AddDI(fDI,fname,nStp,cDims,"Vx")
		AddDI(fDI,fname,nStp,cDims,"Vy")

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
		AddDI(fDI,fname,nStp,cDims,"Vx")
		AddDI(fDI,fname,nStp,cDims,"Vy")
		AddDI(fDI,fname,nStp,cDims,"Vz")

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
		AddDI(fDI,fname,nStp,cDims,"Bx")
		AddDI(fDI,fname,nStp,cDims,"By")

	if ( (Nd == 3) and ("Bx" in vIds) and ("By" in vIds) and ("Bz" in vIds)):
		vAtt = et.SubElement(Grid,"Attribute")
		vAtt.set("Name","VecB")
		vAtt.set("AttributeType","Vector")
		vAtt.set("Center","Cell")
		fDI = et.SubElement(vAtt,"DataItem")
		fDI.set("ItemType","Function")
		fDI.set("Dimensions",vDims)
		fDI.set("Function","JOIN($0 , $1 , $2)")
		AddDI(fDI,fname,nStp,cDims,"Bx")
		AddDI(fDI,fname,nStp,cDims,"By")
		AddDI(fDI,fname,nStp,cDims,"Bz")	
#Decide on centering
def getLoc(gDims,vDims):
	vDims = np.array(vDims,dtype=np.int)
	
	Nd = len(vDims)
	if (Nd == 3):
		Ngx,Ngy,Ngz = gDims-1
		Nvx,Nvy,Nvz = vDims
	elif (Nd==2):
		Ngx,Ngy = gDims-1
		Nvx,Nvy = vDims
		Ngz = 0
		Nvz = 0
	else:
		print("Not enough dimensions!")
		quit()
	#For now just testing for cell centers
	if ( (Ngx == Nvx) and (Ngy == Nvy) and (Ngz == Nvz) ):
		vLoc = "Cell" #Cell-centered data
	elif ( (Ngx == Nvx+1) and (Ngy == Nvy+1) and (Ngz == Nvz+1) ):
		vLoc = "Node"
	else:
		vLoc = "Other"
	return vLoc