import paraview.simple as pvs


def pfs_UT(ut0Str):
	print(ut0Str)
	retStr = r"""
import datetime
isotfmt = '%Y-%m-%dT%H:%M:%S'
ut0Str = '{}'
ut0 = datetime.datetime.strptime(ut0Str,isotfmt)
timeUT = vtk.vtkStringArray()
timeUT.SetName('UT')
t = inputs[0].GetInformation().Get(vtk.vtkDataObject.DATA_TIME_STEP())
t = int(t)  # Remove decimals
timeAsString = str(ut0+datetime.timedelta(seconds=t))
timeUT.InsertNextValue(timeAsString)

self.GetOutput().GetFieldData().AddArray(timeUT)

""".format(ut0Str)
	return retStr

def pf_UT(base, ut0Str):
	pf = pvs.ProgrammableFilter(Input=base)
	pf.Script = pfs_UT(ut0Str)
	pf.RequestInformationScript = ''
	pf.RequestUpdateExtentScript = ''
	pf.PythonPath = ''
	pf.CopyArrays=1

	pvs.RenameSource('pf_UT', pf)

def f_calcSpeeds(base):
	""" Clearly not implemented yet
		Will eventually calc sound speed, alfven speed, and magnetosonic mach
		  and add to dataset
	"""



	scriptStr = """input0 = inputs[0]
bby = input0.CellData['By']*2
output.CellData.append(bby,'Bby')"""

	pf = pvs.ProgrammableFilter(Input=base)
	pf.Script = scriptStr
	pf.RequestInformationScript = ''
	pf.RequestUpdateExtentScript = ''
	pf.PythonPath = ''
	pf.CopyArrays=1

	# Show
	renderView1 = pvs.GetActiveViewOrCreate('RenderView')
	programmableFilter2Display = pvs.Show(pf, renderView1, 'GeometryRepresentation')