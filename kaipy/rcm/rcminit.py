#Various routines to generate RCM data
import numpy as np

def LoadLAS1(fIn="rcmlas1"):
	import os
	__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

	fIn = os.path.join(__location__,fIn)

	print("Reading %s"%(fIn))
		
	Q = np.loadtxt(fIn,skiprows=2)

	alamc = Q[:,0]
	etac = Q[:,1]
	ikflavc = Q[:,2]
	fudgec = Q[:,3]

	return alamc,etac,ikflavc,fudgec

def LoadDKTab(fIn="dktable"):
	import os
	__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

	fIn = os.path.join(__location__,fIn)

	print("Reading %s"%(fIn))

	Q = np.loadtxt(fIn)
	
	return Q.flatten()
	
def LoadEnchan(fIn="enchan.dat"):
	import os
	__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

	fIn = os.path.join(__location__,fIn)

	print("Reading %s"%(fIn))

	Q = np.loadtxt(fIn,skiprows=1)

	iflavin = Q[:,0]
	alamin = Q[:,1]

	return iflavin,alamin
		