import paraview.simple as pvs

def pf_calcSpeeds(base):

	scriptStr = """
in0 = inputs[0]

bmag = sqrt(in0.CellData['Bx']**2+in0.CellData['By']**2+in0.CellData['Bz']**2)
vmag = sqrt(in0.CellData['Vx']**2+in0.CellData['Vy']**2+in0.CellData['Vz']**2)
temp_J = in0.CellData['P']*1E-9/(in0.CellData['D']*1E6)
temp_kT = temp_J/1.380649E-23

cs = sqrt(temp_J/1.6726219E-27)
va = bmag*1E-9 / sqrt(12.566*1E-7 * in0.CellData['D']*1E6 * 1.6726219E-27)
msmach = vmag*1E3 / sqrt(cs**2 + va**2)

output.CellData.append(bmag,'Bmag [nT]')
output.CellData.append(vmag,'Vmag [km/s]')
output.CellData.append(temp_kT,'Temp [kT]')
output.CellData.append(cs,'Cs')
output.CellData.append(va,'Va')
output.CellData.append(msmach,'Magnetosonic Mach')

"""

	pf = pvs.ProgrammableFilter(Input=base)
	pf.Script = scriptStr
	pf.RequestInformationScript = ''
	pf.RequestUpdateExtentScript = ''
	pf.PythonPath = ''
	pf.CopyArrays=1
	pvs.RenameSource('pf_calcSpeeds', pf)
	# Show
	renderView1 = pvs.GetActiveViewOrCreate('RenderView')
	programmableFilter2Display = pvs.Show(pf, renderView1, 'GeometryRepresentation')