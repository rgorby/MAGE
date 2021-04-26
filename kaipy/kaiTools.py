import numpy as np
import datetime
from astropy.time import Time
from xml.dom import minidom

def MJD2UT(mjd):
	astroT = Time(mjd,format='mjd').iso
	utall = []
	for ut in astroT:
		utall.append(datetime.datetime.strptime(ut,'%Y-%m-%d %H:%M:%S.%f'))

	return utall

def genSCXML(scid="sctrack_A",ebfile="msphere",h5traj="sctrack_A.h5"):

	root = minidom.Document()
	xml = root.createElement('Chimp')
	root.appendChild(xml)
	scChild = root.createElement("sim")
	scChild.setAttribute("runid",scid)
	xml.appendChild(scChild)
	fieldsChild = root.createElement("fields")
	fieldsChild.setAttribute("doMHD","T")
	fieldsChild.setAttribute("grType","LFM")
	fieldsChild.setAttribute("ebfile",ebfile)
	xml.appendChild(fieldsChild)
	unitsChild = root.createElement("units")
	unitsChild.setAttribute("uid","EARTH")
	xml.appendChild(unitsChild)
	trajChild = root.createElement("trajectory")
	trajChild.setAttribute("H5Traj",h5traj)
	trajChild.setAttribute("doSmooth","T")
	xml.appendChild(trajChild)

	return root


