import paraview.simple as pvs

"""
Notes:
    It looks like everything in the pipeline is ran for each individual timestep. 
      Meaning anthing that uses the time value (like addUT below) is getting a single value at a time,
      not the entire array of times

"""


def addUT(base, ut0Str):
	scriptStr = r"""
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

	pf = pvs.ProgrammableFilter(Input=base)
	pf.Script = scriptStr
	pf.RequestInformationScript = ''
	pf.RequestUpdateExtentScript = ''
	pf.PythonPath = ''
	pf.CopyArrays=1
	pvs.RenameSource('pf_UT', pf)

	renderView1 = pvs.GetActiveViewOrCreate('RenderView')
	programmableFilter2Display = pvs.Show(pf, renderView1, 'GeometryRepresentation')

	aGD = pvs.AnnotateGlobalData(Input=pf)
	aGD.SelectArrays = 'UT'
	aGD.Prefix = ''
	aGD.Format = '%s'
	pvs.RenameSource('annotateGD_UT', aGD)
	aGDDisplay = pvs.Show(aGD, renderView1, 'TextSourceRepresentation')
