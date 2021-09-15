import os
import json

#Packages needed to handle variable types
import numpy as np
import datetime

dtformat = "%Y-%m-%dT%H:%M:%SZ"

#======
#Custom handlers for non-standard types
#======

#TODO: Handle saving/loading attributes for hdf5 data
class CustomEncoder(json.JSONEncoder):
	def default(self, obj):
		if isinstance(obj, (datetime.time, datetime.datetime)):
			return obj.strftime(dtformat)

		if isinstance(obj, np.ndarray):
			return {'shape': obj.shape, 'data': obj.tolist()}

		#Now you don't need to care if you have floats or numpy floats
		if isinstance(obj, np.float32):
			return "_f32_{:.16f}".format(obj)

		return json.JSONEncoder.default(self, obj)

def customhook(dct):
	for key in dct.keys():
		#Handle datetime
		try:
			datetime.datetime.strptime(dct[key][0], dtformat)
			#If we're still here, it worked
			#So go ahead and replace this whole list with the proper datetime objects
			newlist = [datetime.datetime.strptime(dtStr, dtformat) for dtStr in dct[key]]
			dct[key] = newlist
		except:
			pass

		#Handle numpy arrays
		try:
			if 'shape' in dct[key].keys():
				shape = tuple(dct[key]['shape'])
				newdata = np.array(dct[key]['data']).reshape(shape)
				dct[key] = newdata
		except:
			pass

		#Handle floats
		try:
			if type(dct[key]) == str and '_f32_' in dct[key]:
				dct[key] = float(dct[key].split('_f32_')[1])
		except:
			pass
		
	return dct

#======
#Main functions
#======
def dump(fname, data, action='w'):
	"""Store data [dict] in file fname
		action: (over)'w'rite or 'a'ppend
	"""
	with open(fname, action) as jfile:
		json.dump(data, jfile, indent=4, cls=CustomEncoder)

def load(fname):
	#Pull data[dict] from file

	if not os.path.exists(fname):
		print("File " + fname + " doesn't exist, can't load json")
		return

	with open(fname, 'r') as jfile:
		data = json.load(jfile, object_hook=customhook)

	return data
