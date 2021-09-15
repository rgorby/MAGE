
import datetime
import kaipy.satcomp.scutils as scutils

#TODO: Need to add "epoch" str for each dataset

def test_getscIds():
	print("Testing ability to grab spacecraft Ids from json file")
	scIdDict = scutils.getScIds()
	assert type(scIdDict) == dict, "Returned type is {}, but should be type dict".format(type(scIdDict))
	assert len(scIdDict.keys()) != 0, "Dictionary has zero entries"
	#Check if every spacecraft entry has any data at all
	#Check if every spacecraft entry has an "Ephem" data product
	#Check if every spacefraft data product has at least an "Id" and "data" k-v pair

def test_getCdasData():
	print("Testing if all data in scId dict is retrievable from cdasws")

	scIdDict = scutils.getScIds()

	for scName in scIdDict.keys():
		print(" " + scName)
		scStrs = scIdDict[scName]
		
		for dpStr in scStrs.keys(): 
			if dpStr == '_testing': continue
			print("  " + dpStr, end=" : ")
		
			#Get valid time interval for dataset
			tStart, tEnd = scutils.getCdasDsetInterval(scStrs[dpStr]['Id'])
			#assert tStart != None, "getCdasDsetInterval returned null for dataset " + dpStr
			if tStart is None: 
				print("Bad dset")
				continue

			t0 = tStart
			t0dt = datetime.datetime.strptime(t0, "%Y-%m-%dT%H:%M:%S.%fZ")
			t1 = (t0dt + datetime.timedelta(days=1)).strftime("%Y-%m-%dT%H:%M:%SZ")

			dset_id = scStrs[dpStr]['Id']
			dset_vname = scStrs[dpStr]['Data']

			if "EpochStr" in scStrs[dpStr].keys():
				cdasResult = scutils.getCdasData(dset_id, dset_vname, t0,t1, epochStr=scStrs[dpStr]["EpochStr"])
			else:
				cdasResult = scutils.getCdasData(dset_id, dset_vname, t0,t1)
			assert cdasResult != {}, "getCdasData returned with no information"

			print("Good")

if __name__ == "__main__":
	test_getscIds()

	test_getCdasData()
