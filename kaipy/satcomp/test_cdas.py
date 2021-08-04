
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
	t0, t1 = ['2013-03-18T05:00:00Z','2013-03-19T05:00:00Z']

	for scName in scIdDict.keys():
		print(" " + scName)
		scStrs = scIdDict[scName]
		for dpStr in scStrs.keys(): 
			print("  " + dpStr, end=" : ")
			dset_id = scStrs[dpStr]['Id']
			dset_vname = scStrs[dpStr]['Data']
			cdasResult = scutils.getCdasData(dset_id, dset_vname, t0,t1)

			assert cdasResult != {}, "getCdasData returned with no information"

			print("Good")

if __name__ == "__main__":
	test_getscIds()

	test_getCdasData()
