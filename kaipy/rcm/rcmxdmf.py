"""
RCM-specific tools to include h5 data in xmf files
"""
import h5py as h5
import kaipy.kaixdmf as kxmf


"""
def AddGrid(fname,Geom,iDims,coordStrs):
	
	for coordStr in coordStrs:
		xC = et.SubElement(Geom,"DataItem")
		xC.set("Dimensions",iDims)
		xC.set("NumberType","Float")
		xC.set("Precision","4")
		xC.set("Format","HDF")
		
		text = fname+":/"+coordStr
		xC.text = text
"""

def addRCMGeom(Geom, dimInfo, rcmInfo, sID):
	#Unpack
	rcmh5fname = rcmInfo['rcmh5fname']

	sIDstr = "Step#" + str(sID)

	gridVars = dimInfo['gridVars']  # Expect something like rcmxmin,rcmymin,alamc
	gDims = dimInfo['gDims']  # Someone before us should set kji dims
	#gDims[1] -= 2  # Remove 2 j vals (jWrap)

	#all
	dSetDimStr = ' '.join([str(v) for v in gDims])

	hsDimStr = "3 2"  # row, col
	startStr = "2 0"
	strideStr = "1 1"
	numStr = "359 180"
	origDSetDimStr = str(gDims[0])

	fileText = "{}:/{}/{}".format(rcmh5fname,sIDstr,'rxmxmin')
	kxmf.addHyperslab_DI(Geom, dSetDimStr, hsDimStr, startStr, strideStr, numStr, origDSetDimStr, fileText)	

	fileText = "{}:/{}/{}".format(rcmh5fname,sIDstr,'rcmymin')
	kxmf.addHyperslab_DI(Geom, dSetDimStr, hsDimStr, startStr, strideStr, numStr, origDSetDimStr, fileText)	

	hsDimStr = "3 1"
	startStr = "0"
	strideStr = "1"
	numStr = "140"
	fileText = "{}:/{}/{}".format(rcmh5fname,sIDstr,'alamc')
	kxmf.addHyperslab_DI(Geom, dSetDimStr, hsDimStr, startStr, strideStr, numStr, origDSetDimStr, fileText)	


def addRCMVars(Grid, dimInfo, rcmInfo, sID):  # Used with the mhdrcm presets to add links to given RCM variables

	sIDstr = "Step#" + str(sID) 

	trg_vDims = dimInfo['vDims']  # target file var dims (probably mhdrcm or rcm)
	trg_vDimStr = ' '.join([str(v) for v in trg_vDims])
	trg_nDims = len(trg_vDims)
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
	
		if (r_nDims == 2 and trg_vDims[0] < Nj) or (r_nDims == 3 and trg_vDims[1] < Nj):
			doHyperslab = True
			dimTrim = (Nj - trg_vDims[0]) if trg_nDims == 2 else (Nj - trg_vDims[1])  # Shortening j-dir

		if r_nDims == 2 and doHyperslab == False:  # Easy add
				#print("Adding " + vName)
				kxmf.AddData(Grid,rcmh5fname, vName,"Cell",trg_vDimStr,sID)
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
				kxmf.addHyperslab_Attr(Grid,vName,trg_vDimStr,dimStr,startStr,strideStr,numStr,r_vDimStr,text)
				continue
			elif r_nDims == 3:
				dimStr = "3 3"
				strideStr = str(Nk+1) + " 1 1"
				numStr = "1 {} {}".format(Nj-dimTrim,Ni)
				for k in rcmKs:
					startStr = "{} {} 0".format(k,dimTrim)
					vName_k = vName + "_k{}".format(k)
					kxmf.addHyperslab_Attr(Grid,vName_k,trg_vDimStr,dimStr,startStr,strideStr,numStr,r_vDimStr,text)
