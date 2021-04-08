import numpy as np
import datetime
from astropy.time import Time

def MJD2UT(mjd):
	astroT = Time(mjd,format='mjd').iso
	utall = []
	for ut in astroT:
		utall.append(datetime.datetime.strptime(ut,'%Y-%m-%d %H:%M:%S.%f'))

	return utall


