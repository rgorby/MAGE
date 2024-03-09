import paraview.simple as pvs


# Helper
def makePF(base, scriptStr, filterName):
	pf = pvs.ProgrammableFilter(Input=base)
	pf.Script = scriptStr
	pf.RequestInformationScript = ''
	pf.RequestUpdateExtentScript = ''
	pf.PythonPath = ''
	pf.CopyArrays=1
	pvs.RenameSource(scriptStr, pf)
	# Show
	renderView1 = pvs.GetActiveViewOrCreate('RenderView')
	programmableFilter2Display = pvs.Show(pf, renderView1, 'GeometryRepresentation')

def pf_calcSpeeds(base):

	scriptStr = """
import kaipy.kdefs as kdefs
in0 = inputs[0]

bmag = sqrt(in0.CellData['Bx']**2+in0.CellData['By']**2+in0.CellData['Bz']**2)
vmag = sqrt(in0.CellData['Vx']**2+in0.CellData['Vy']**2+in0.CellData['Vz']**2)
temp_J = in0.CellData['P']*1E-9/(in0.CellData['D']*1E6)
temp_kT = temp_J/(kdefs.kbltz*kdefs.erg2J)
temp_eV = temp_J/kdefs.ev2J

cs = sqrt(temp_J/(kdefs.Mp_cgs*1e-3))
va = bmag*1E-9 / sqrt(kdefs.Mu0 * in0.CellData['D']*1E6 * kdefs.Mp_cgs*1e-3)
msmach = vmag*1E3 / sqrt(cs**2 + va**2)

output.CellData.append(bmag,'Bmag [nT]')
output.CellData.append(vmag,'Vmag [km/s]')
output.CellData.append(temp_kT,'Temp [kT]')
output.CellData.append(temp_eV,'Temp [eV]')
output.CellData.append(cs,'Cs')
output.CellData.append(va,'Va')
output.CellData.append(msmach,'Magnetosonic Mach')

"""
	makePF(base, scriptStr, 'pf_calcSpeeds')


def pf_calcEfield(base):

	scriptStr="""
import kaipy.kdefs as kdefs
import kaipy.kaiTools as kt
in0 = inputs[0]

Bx = in0.CellData['Bx']
By = in0.CellData['By']
Bz = in0.CellData['Bz']
Vx = in0.CellData['Vx']
Vy = in0.CellData['Vy']
Vz = in0.CellData['Vz']

#E=-VxB
Ex = -(Vy*Bz-Vz*By)*0.001 # [mV/m]
Ey =  (Vx*Bz-Vz*Bx)*0.001
Ez = -(Vx*By-Vy*Bx)*0.001

xxx

output.CellData.append(Ex,'Ex')
output.CellData.append(Ey,'Ey')
output.CellData.append(Ez,'Ez')
"""
	makePF(base, scriptStr, 'pf_Efield')
	
#------

"""
Bx = gsph.EggSlice("Bx",nStp,doEq=True)
	By = gsph.EggSlice("By",nStp,doEq=True)
	Bz = gsph.EggSlice("Bz",nStp,doEq=True)
	Vx = gsph.EggSlice("Vx",nStp,doEq=True)
	Vy = gsph.EggSlice("Vy",nStp,doEq=True)
	Vz = gsph.EggSlice("Vz",nStp,doEq=True)

	# calculating some variables to to plot
	#E=-VxB
	Ex = -(Vy*Bz-Vz*By)*0.001 # [mV/m]
	Ey =  (Vx*Bz-Vz*Bx)*0.001
	Ez = -(Vx*By-Vy*Bx)*0.001

	# coordinate transform
	ppc = np.arctan2(gsph.yyc,gsph.xxc)
	theta = np.pi #eq plane
	Er,Et,Ep = kt.xyz2rtp(ppc,theta,Ex,Ey,Ez)
"""