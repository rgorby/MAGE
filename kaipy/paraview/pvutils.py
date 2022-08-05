import paraview.simple as pvs

"""
Notes:
    It looks like everything in the pipeline is ran for each individual timestep. 
      Meaning anything that uses the time value (like addTime below) is getting a single value at a time,
      not the entire array of times

"""
defaultUT = '2000-01-01T00:00:00'

def addTime(base, fmt="h", ut0Str=defaultUT):
	"""Adds a label showing the simulation time
		fmt options: s,m,h,UT
		ut0Str: If fmt=UT, associates this argument with simulation t=0 and uses it to calculate subsequent UT values

		TODO: Use the run's MJD0 to get ut0Str. Need to manually read the h5 file, meaning more arguments to provide.
	"""

	# "Therefore unlike other programming languages, Python gives meaning to the indentation rule, which is very simple and also helps to make the code readable."
	# https://www.educba.com/indentation-in-python/

	scriptStr = r"""
timeVal = vtk.vtkStringArray()
timeVal.SetName('{:s}')
t = inputs[0].GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())

{:s}

timeVal.InsertNextValue(timeAsString)
self.GetOutput().GetFieldData().AddArray(timeVal)

"""
	varName = "time_"+fmt

	#Here we choose how to take the variable <t> and turn it into our desired <timeAsString> to be displayed
	if fmt=="UT":
		if ut0Str==defaultUT:
			print("No ut0Str specified, defaulting to " + defaultUT)
		prefStr = ''
		tScr = r"""
import datetime
isotfmt = '%Y-%m-%dT%H:%M:%S'
ut0Str = '{}'
ut0 = datetime.datetime.strptime(ut0Str,isotfmt)
t = int(t)  # Remove decimals
timeAsString = str(ut0+datetime.timedelta(seconds=t))
		""".format(ut0Str)

	elif fmt=="s":
		prefStr = 'Second '
		tScr = r"""timeAsString = "{:4.2f}".format(int(t)) """
	elif fmt=="m":
		prefStr = 'Minute '
		tScr = r"""timeAsString = "{:4.2f}".format(t/60) """
	elif fmt=="h":
		prefStr = 'Hour '
		tScr = r"""timeAsString = "{:4.2f}".format(t/3600) """



	#Create the programmable filter
	pf = pvs.ProgrammableFilter(Input=base)
	#Now insert our specific code into the script block
	scriptStr = scriptStr.format(varName, tScr)
	pf.Script = scriptStr
	#Set other parameters
	pf.RequestInformationScript = ''
	pf.RequestUpdateExtentScript = ''
	pf.PythonPath = ''
	pf.CopyArrays=1
	pvs.RenameSource('pf_Time', pf)

	renderView1 = pvs.GetActiveViewOrCreate('RenderView')
	pvs.Show(pf, renderView1, 'GeometryRepresentation')

	aGD = pvs.AnnotateGlobalData(Input=pf)
	aGD.SelectArrays = varName
	aGD.Prefix = prefStr
	aGD.Format = '%s'
	pvs.RenameSource('annotateGD_Time', aGD)
	pvs.Show(aGD, renderView1, 'TextSourceRepresentation')

	pvs.Hide(pf, renderView1)  # Turn off programmable filter view, no need to view after we have the annotation
