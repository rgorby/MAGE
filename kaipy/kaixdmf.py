#Various routines to help scripts that create XMF files from H5 data
import numpy as np
import h5py
import kaipy.kaiH5 as kh5
import xml.etree.ElementTree as et

#Add grid info to step
#Geom is topology subelement, iDims is grid size string
def AddGrid(fname,Geom,iDims,coordStrs):
	
	for coordStr in coordStrs:
		xC = et.SubElement(Geom,"DataItem")
		xC.set("Dimensions",iDims)
		xC.set("NumberType","Float")
		xC.set("Precision","4")
		xC.set("Format","HDF")
		
		text = fname+":/"+coordStr
		xC.text = text

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
def getRootVars(fname,gDims):
	#Dims = kh5.getDims(fname,doFlip=False) #Do kji ordering
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
				vLoc = getLoc(gDims,Nv)
				if (vLoc != "Other"):
					vIds.append(vID)
					vLocs.append(vLoc)
				else:
					print("Excluding %s"%(vID))

	return vIds,vLocs

#Get variables in initial Step
def getVars(fname,gId,gDims,recursive=False):

	def getVarName(gId, k, recursive):
		if recursive:
			# Want to keep the group name, since all vars need to be at the same level in xdmf
			spl = gId.split('/')
			return '/'.join(spl[1:]) + '/' + str(k)
		else:
			return str(k)


	with h5py.File(fname,'r') as hf:
		stp0 = hf[gId]
		vIds = []
		vLocs = []
		for k in stp0.keys():
			vID = getVarName(gId,k,recursive)
			
			# If its a group, recursively call this function to get to the vars and add them to the list
			if type(stp0[k]) == h5py._hl.group.Group:
				vIdR, vLR = getVars(fname, gId+'/'+k, gDims, recursive=True)
				for vi, vl in zip(vIdR, vLR):
					vIds.append(vi)
					vLocs.append(vl)
			else:
				Nv = stp0[k].shape
				vLoc = getLoc(gDims,Nv)
				if (vLoc != "Other"):
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
	dimLocs = []
	if len(gDims) != len(vDims):
		return "Other"
	for d in range(len(gDims)):
		Ngd = gDims[d]-1
		Nvd = vDims[d]
		if Ngd == Nvd:
			dimLocs.append("Cell")
		elif Ngd == Nvd-1:
			dimLocs.append("Node")
		else:
			dimLocs.append("Other")

	if "Other" in dimLocs:
		return "Other"
	#If all the same, we have consensus
	if all(x == dimLocs[0] for x in dimLocs):
		return dimLocs[0]
	else:
		return "Other"

def addHyperslab(Grid,vName,dSetDimStr,vdimStr,startStr,strideStr,numStr,origDSetDimStr,fileText):
	vAtt = et.SubElement(Grid, "Attribute")
	vAtt.set("Name",vName)
	vAtt.set("AttributeType","Scalar")
	vAtt.set("Center","Node")
	slabDI = et.SubElement(vAtt, "DataItem")
	slabDI.set("ItemType","HyperSlab")
	slabDI.set("Dimensions",dSetDimStr)
	slabDI.set("Type","HyperSlab")  # Not sure if redundant, but it works
	cutDI = et.SubElement(slabDI,"DataItem")
	cutDI.set("Dimensions",vdimStr)
	cutDI.set("Format","XML")
	cutDI.text = "\n{}\n{}\n{}\n".format(startStr,strideStr,numStr)
	datDI = et.SubElement(slabDI,"DataItem")
	datDI.set("Dimensions",origDSetDimStr)
	datDI.set("DataType","Float")
	datDI.set("Precision","4")
	datDI.set("Format","HDF")
	datDI.text = fileText
