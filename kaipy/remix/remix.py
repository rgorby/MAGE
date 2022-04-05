import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys

Ri = 6.5  # radius of ionosphere in 1000km
Re = 6.38 # radius of Earth in 1000km
mu0o4pi = 1.e-7 # mu0=4pi*10^-7 => mu0/4pi=10^-7

#Setting globals here to grab them from other modules
facMax = 1.5
facCM = cm.RdBu_r
flxCM = cm.inferno

class remix:
	def __init__(self,h5file,step):
		# create the ion object to store data and coordinates
		self.ion = self.get_data(h5file,step)
		self.Initialized=False

		# DEFINE DATA LIMITS
		self.variables = { 'potential' : {'min':-100,
										'max': 100},
						 'current'   : {'min':-facMax,
										'max':facMax},
						 'sigmap'    : {'min':1,
										'max':20},
						 'sigmah'    : {'min':2,
										'max':40},
						 'energy'    : {'min':0,
										'max':20},
						 'flux'      : {'min':0,
										'max':1.e9},
						 'eflux'     : {'min':0,
										'max':10.},
						 'efield'    : {'min':-1,
										'max':1},
						 'joule'     : {'min':0,
										'max':10},
						 'jhall'     : {'min':-2,
						 				'max':2},
						 'Menergy'    : {'min':0,
										'max':20},
						 'Mflux'      : {'min':0,
										'max':1.e10},
						 'Meflux'     : {'min':0,
										'max':10.},
						 'Denergy'    : {'min':0,
										'max':20},
						 'Dflux'      : {'min':0,
										'max':1.e10},
						 'Deflux'     : {'min':0,
										'max':10.},
						 'Penergy'    : {'min':0,
										'max':240},
						 'Pflux'      : {'min':0,
										'max':1.e7},
						 'Peflux'     : {'min':0,
										'max':1.}
						 }

	def get_data(self,h5file,step):
		ion = {}

		with h5py.File(h5file,'r') as f:
			ion['X'] = f['X'][:]
			ion['Y'] = f['Y'][:]
			for h in f['Step#%d'%step].keys():
				ion[h] = f['Step#%d'%step][h][:]

		# Get spherical coords
		ion['R'],ion['THETA'] = self.get_spherical(ion['X'],ion['Y'])
						
		return ion

	# TODO: check for variable names passed to plot
	def init_vars(self,hemisphere):
		h = hemisphere # for shortness

		if (h.lower()=='north'):
			self.variables['potential']['data'] = self.ion['Potential '+h]
			self.variables['current']['data']   =-self.ion['Field-aligned current '+h]  # note, converting to common convention (upward=positive)
			self.variables['sigmap']['data']    = self.ion['Pedersen conductance '+h]
			self.variables['sigmah']['data']    = self.ion['Hall conductance '+h]
			self.variables['energy']['data']    = self.ion['Average energy '+h]
			self.variables['flux']['data']      = self.ion['Number flux '+h]
			# variables['efield']['data']    = efield_n*1.e6
			# variables['joule']['data']     = sigmap_n*efield_n**2*1.e-3
			if 'Zhang average energy '+h in self.ion.keys():
				self.variables['Menergy']['data'] = self.ion['Zhang average energy '+h]
				self.variables['Mflux']['data']   = self.ion['Zhang number flux '+h]
				self.variables['Meflux']['data'] = self.variables['Menergy']['data']*self.variables['Mflux']['data']*1.6e-9
			if 'IM average energy '+h in self.ion.keys():
				self.variables['Denergy']['data'] = self.ion['IM average energy '+h]
				self.variables['Deflux']['data']  = self.ion['IM Energy flux '+h]
				self.variables['Denergy']['data'][self.variables['Denergy']['data']==0] = 1.e-20
				self.variables['Denergy']['data'][np.isnan(self.variables['Denergy']['data'])] = 1.e-20
				self.variables['Dflux']['data']   = self.variables['Deflux']['data']/self.variables['Denergy']['data']/(1.6e-9)
				self.variables['Dflux']['data'][self.variables['Denergy']['data']==1.e-20] = 0.
			if 'IM average energy proton '+h in self.ion.keys():
				self.variables['Penergy']['data'] = self.ion['IM average energy proton '+h]
				self.variables['Peflux']['data']  = self.ion['IM Energy flux proton '+h]
				self.variables['Penergy']['data'][self.variables['Penergy']['data']==0] = 1.e-20
				self.variables['Penergy']['data'][np.isnan(self.variables['Penergy']['data'])] = 1.e-20
				self.variables['Pflux']['data']   = self.variables['Peflux']['data']/self.variables['Penergy']['data']/(1.6e-9)
				self.variables['Pflux']['data'][self.variables['Penergy']['data']==1.e-20] = 0.
		else:  # note flipping the y(phi)-axis
			self.variables['potential']['data'] = self.ion['Potential '+h][:,::-1]
			self.variables['current']['data']   = self.ion['Field-aligned current '+h][:,::-1]
			self.variables['sigmap']['data']    = self.ion['Pedersen conductance '+h][:,::-1]
			self.variables['sigmah']['data']    = self.ion['Hall conductance '+h][:,::-1]
			self.variables['energy']['data']    = self.ion['Average energy '+h][:,::-1]
			self.variables['flux']['data']      = self.ion['Number flux '+h][:,::-1]
			if 'Zhang average energy '+h in self.ion.keys():
				self.variables['Menergy']['data'] = self.ion['Zhang average energy '+h][:,::-1]
				self.variables['Mflux']['data']   = self.ion['Zhang number flux '+h][:,::-1]
				self.variables['Meflux']['data'] = self.variables['Menergy']['data']*self.variables['Mflux']['data']*1.6e-9
			if 'IM average energy '+h in self.ion.keys():
				self.variables['Denergy']['data'] = self.ion['IM average energy '+h][:,::-1]
				self.variables['Deflux']['data']  = self.ion['IM Energy flux '+h][:,::-1]
				self.variables['Denergy']['data'][self.variables['Denergy']['data']==0] = 1.e-20
				self.variables['Denergy']['data'][np.isnan(self.variables['Denergy']['data'])] = 1.e-20
				self.variables['Dflux']['data']   = self.variables['Deflux']['data']/self.variables['Denergy']['data']/(1.6e-9)
				self.variables['Dflux']['data'][self.variables['Denergy']['data']==1.e-20] = 0.
			if 'IM average energy proton '+h in self.ion.keys():
				self.variables['Penergy']['data'] = self.ion['IM average energy proton '+h][:,::-1]
				self.variables['Peflux']['data']  = self.ion['IM Energy flux proton '+h][:,::-1]
				self.variables['Penergy']['data'][self.variables['Penergy']['data']==0] = 1.e-20
				self.variables['Penergy']['data'][np.isnan(self.variables['Penergy']['data'])] = 1.e-20
				self.variables['Pflux']['data']   = self.variables['Peflux']['data']/self.variables['Penergy']['data']/(1.6e-9)
				self.variables['Pflux']['data'][self.variables['Penergy']['data']==1.e-20] = 0.

		# convert energy flux to erg/cm2/s to conform to Newell++, doi:10.1029/2009JA014326, 2009
		self.variables['eflux']['data'] = self.variables['energy']['data']*self.variables['flux']['data']*1.6e-9 
		self.Initialized=True


	def get_spherical(self,x,y):
		# note, because of how the grid is set up in the h5 file,
		# the code below produces the first theta that's just shy of 2pi.
		# this is because the original grid is staggered at half-cells from data.
		# empirically, this is OK for pcolormesh plots under remix.plot.
		# but I still fix it manually.
		theta=np.arctan2(y,x)
		theta[theta<0]=theta[theta<0]+2*np.pi
		theta[:,0] -= 2*np.pi  # fixing the first theta point to just below 0
		r=np.sqrt(x**2+y**2)

		return(r,theta)
	def distance(self,p0, p1):
	  """ Calculate the distance between p0 & p1--two points in R**3 """
	  return( np.sqrt( (p0[0]-p1[0])**2 + 
                      (p0[1]-p1[1])**2 + 
                      (p0[2]-p1[2])**2 ) )

	def calcFaceAreas(self,x,y):
	  """ 
	  Calculate the area of each face in a quad mesh.
	 
	  >>> calcFaceAreas(numpy.array([[ 0., 1.],[ 1., 0.]]), numpy.array([[ 0.,1.],[ 1.,0.]]), numpy.array([[ 0., 0.],[ 0.,0.]]))
	  array([[ 2.]])
	  """
	  (nLonP1, nLatP1) = x.shape
	  (nLon, nLat) = (nLonP1-1, nLatP1-1)
	  z=np.sqrt(1.0-x**2-y**2)

	  area = np.zeros((nLon, nLat))

	  for i in range(nLon):
	    for j in range(nLat):
	      left  = self.distance( (x[i,j],   y[i,j],   z[i,j]),   (x[i,j+1],   y[i,j+1],   z[i,j+1]) )
	      right = self.distance( (x[i+1,j], y[i+1,j], z[i+1,j]), (x[i+1,j+1], y[i+1,j+1], z[i+1,j+1]) )
	      top =   self.distance( (x[i,j+1], y[i,j+1], z[i,j+1]), (x[i+1,j+1], y[i+1,j+1], z[i+1,j+1]) )
	      bot =   self.distance( (x[i,j],   y[i,j],   z[i,j]),   (x[i+1,j],   y[i+1,j],   z[i+1,j]) )
	      
	      area[i,j] = 0.5*(left+right) * 0.5*(top+bot)

	  return area

	# TODO: define and consolidate allowed variable names
	def plot(self,varname,
			 ncontours=16,   # default number of potential contours
			 addlabels={},
			 gs=None,doInset=False,doCB=True,doCBVert=True):

		# define function for potential contour overplotting
		# to keep code below clean and compact
		def potential_overplot(doInset=False):
			tc = 0.25*(theta[:-1,:-1]+theta[1:,:-1]+theta[:-1,1:]+theta[1:,1:])
			rc = 0.25*(r[:-1,:-1]+r[1:,:-1]+r[:-1,1:]+r[1:,1:])

			# trick to plot contours smoothly across the periodic boundary:
			# wrap around: note, careful with theta -- need to add 2*pi to keep it ascending
			# otherwise, contours mess up
			tc = np.hstack([tc,2.*np.pi+tc[:,[0]]])
			rc = np.hstack([rc,rc[:,[0]]])
			tmp=self.variables['potential']['data']
			tmp = np.hstack([tmp,tmp[:,[0]]])

			# similar trick to make contours go through the pole
			# add pole
			tc = np.vstack([tc[[0],:],tc])
			rc = np.vstack([0.*rc[[0],:],rc])			
			tmp = np.vstack([tmp[0,:].mean()*np.ones_like(tmp[[0],:]),tmp])						

			# finally, plot
			if (doInset):
				LW = 0.25
				alpha = 1
				tOff = 0.0
			else:
				LW = 0.5
				alpha = 1
				tOff = np.pi/2.
			ax.contour(tc+tOff,rc,tmp,15,colors='black',linewidths=LW,alpha=alpha)

			if (not doInset):
				# also, print min/max values of the potential
				ax.text(73.*np.pi/180.,1.03*r.max(),('min: '+format_str+'\nmax: ' +format_str) % 
					  (tmp.min() ,tmp.max()))


		if not self.Initialized:
			sys.exit("Variables should be initialized for the specific hemisphere (call init_var) prior to plotting.")

		# Aliases to keep things short
		x = self.ion['X']
		y = self.ion['Y']
		r = self.ion['R']
		theta = self.ion['THETA']

		# List all possible variable names here, add more from the input parameter if needed
		cblabels = {'potential' : r'Potential [kV]',
					'current'   : r'Current density [$\mu$A/m$^2$]',
					'sigmap'    : r'Pedersen conductance [S]',
					'sigmah'    : r'Hall conductance [S]',
					'energy'    : r'Energy [keV]',
					'flux'      : r'Flux [1/cm$^2$s]',
					'eflux'     : r'Energy flux [erg/cm$^2$s]',
					'ephi'      : r'$E_\phi$ [mV/m]',
					'etheta'    : r'$E_\theta$ [mV/m]',
					'efield'    : r'|E| [mV/m]',
					'joule'     : r'Joule heating [mW/m$^2$]',
					'jped'      : r'Pedersen current [$\mu$A/m]',
					'magpert'   : r'Magnetic perturbation [nT]',
					'Menergy'    : r'Mono Energy [keV]',
					'Mflux'      : r'Mono Flux [1/cm$^2$s]',
					'Meflux'     : r'Mono Energy flux [erg/cm$^2$s]',
					'Denergy'    : r'Diffuse Energy [keV]',
					'Dflux'      : r'Diffuse Flux [1/cm$^2$s]',
					'Deflux'     : r'Diffuse Energy flux [erg/cm$^2$s]',
					'Penergy'    : r'Diffuse Proton Energy [keV]',
					'Pflux'      : r'Diffuse Proton Flux [1/cm$^2$s]',
					'Peflux'     : r'Diffuse Proton Energy flux [erg/cm$^2$s]'
					}
		cblabels.update(addlabels)  # a way to add cb labels directly through function arguments (e.g., for new variables)

		# if limits are given use them, if not use the variables min/max values
		if ('min' in self.variables[varname]):
			lower = self.variables[varname]['min']
		else:
			lower = self.variables[varname]['data'].min()

		if ('max' in self.variables[varname]):
			upper = self.variables[varname]['max']
		else:
			upper = self.variables[varname]['data'].max()

		# define number format string for min/max labels
		if varname !='flux' and varname !='Pflux': 
			if varname == 'jped':
				format_str = '%.2f'
			else:
				format_str = '%.1f'
		else:
			format_str = '%.1e'

		if (varname == 'potential') or (varname == 'current'):
			cmap=facCM
		elif varname in ['flux','energy','eflux','joule','Mflux','Menergy','Meflux','Dflux','Denergy','Deflux','Pflux','Penergy','Peflux']:
			cmap=flxCM
			latlblclr = 'white'
		elif (varname == 'velocity'): 
			cmap=cm.YlOrRd
		elif (varname == 'efield'): 
			cmap=None #cm.RdBu_r#None # cm.hsv
		elif (varname == 'magpert'): 
			cmap = cm.YlOrRd
		else:
			cmap=None # default is used
		
		if varname in ['eflux','Meflux','Deflux','Peflux']:
			ri=6500.e3
			areaMixGrid = self.calcFaceAreas(x,y)*ri*ri
			hp = areaMixGrid*self.variables[varname]['data'][:,:]/(1.6e-9)
			power = hp.sum()*1.6e-21
		if (varname == 'current'):
			ri=6500.e3
			areaMixGrid = self.calcFaceAreas(x,y)*ri*ri
			fac = self.variables[varname]['data'][:,:]
			dfac = areaMixGrid[fac>0.]*fac[fac>0.]/1.0e12
			pfac = dfac.sum()

		# DEFINE GRID LINES AND LABELS
		if (doInset):
			circle_list = [15,30,45]
			lbls = ["","",str(45)+u'\xb0']
			hour_labels = ["","","",""]
			LW = 0.25
			latlblclr = 'silver'
			tOff = 0.0
		else:
			circle_list = [10,20,30,40]
		# convert to string and add degree symbol
			lbls = [str(elem)+u'\xb0' for elem in circle_list] 
			hour_labels = ['06','12','18','00']
			LW = 1.0
			latlblclr = 'black'
			tOff = np.pi/2.
		circles = np.sin(np.array(circle_list)*np.pi/180.)

		

		if varname == 'joule':
			self.variables[varname]['data'] = self.joule()*1.e3  # convert to mW/m^2

		variable = self.variables[varname]['data']

		fig = plt.gcf()
		
		# Now plotting
		if gs != None:
			ax=fig.add_subplot(gs,polar=True)
		else:
			ax=fig.add_subplot(polar=True) 

		p=ax.pcolormesh(theta+tOff,r,variable,cmap=cmap,vmin=lower,vmax=upper)

		if (not doInset):
			if (doCB):
				if (doCBVert):
					cb=plt.colorbar(p,ax=ax,pad=0.1,shrink=0.85,orientation="vertical")
				else:
					cb=plt.colorbar(p,ax=ax,pad=0.1,shrink=0.85,orientation="horizontal")
				cb.set_label(cblabels[varname])

				ax.text(-75.*np.pi/180.,1.2*r.max(),('min: '+format_str+'\nmax: ' +format_str) % 
					  (variable.min() ,variable.max()))
			
		lines, labels = plt.thetagrids((0.,90.,180.,270.),hour_labels)
		lines, labels = plt.rgrids(circles,lbls,fontsize=8,color=latlblclr)
		# if (doInset):
		# 	lines, labels = plt.rgrids(circles,lbls,fontsize=8,color=latlblclr)

		ax.grid(True,linewidth=LW)
		# ax.axis([0,2*np.pi,0,r.max()],'tight')
		ax.axis([0,2*np.pi,0,r.max()])
		if (not doInset):
			ax.text(-75.*np.pi/180.,1.2*r.max(),('min: '+format_str+'\nmax: ' +format_str) % 
				  (variable.min() ,variable.max()))
			if varname in ['eflux','Meflux','Deflux','Peflux']:
				ax.text(75.*np.pi/180.,1.3*r.max(),('GW: '+format_str) % power)
			if (varname == 'current'):
				ax.text(-(73.+45.)*np.pi/180.,1.3*r.max(),('MA: '+format_str) % pfac)
		ax.grid(True)

		if varname=='current': 
			potential_overplot(doInset)

		return ax
	# mpl.rcParams['contour.negative_linestyle'] = 'solid'
	# if (varname == 'efield' or varname == 'velocity' or varname =='joule'): 
	#     contour(theta+pi/2.,r,variables['potential']['data'][:,2:-1],21,colors='black')
#                                      arange(variables['potential']['min'],variables['potential']['max'],21.),colors='purple')

	# FIXME: MAKE WORK FOR SOUTH (I THINK IT DOES BUT MAKE SURE)
	def efield(self,returnDeltas=False,ri=Ri*1e3): 
		if not self.Initialized:
			sys.exit("Variables should be initialized for the specific hemisphere (call init_var) prior to efield calculation.")

		Psi = self.variables['potential']['data']  # note, these are numbers of cells. self.ion['X'].shape = Nr+1,Nt+1
		Nt,Np = Psi.shape

		# Aliases to keep things short
		x = self.ion['X']
		y = self.ion['Y']

		# note the change in naming convention from above
		# i.e., theta is now the polar angle
		# and phi is the azimuthal (what was theta)
		# TODO: make consistent throughout
		theta = np.arcsin(self.ion['R'])
		phi   = self.ion['THETA']

		# interpolate Psi to corners
		Psi_c = np.zeros(x.shape)
		Psi_c[1:-1,1:-1] = 0.25*(Psi[1:,1:]+Psi[:-1,1:]+Psi[1:,:-1]+Psi[:-1,:-1])

		# fix up periodic
		Psi_c[1:-1,0]  = 0.25*(Psi[1:,0]+Psi[:-1,0]+Psi[1:,-1]+Psi[:-1,-1])
		Psi_c[1:-1,-1] = Psi_c[1:-1,0]

		# fix up pole
		Psi_pole = Psi[0,:].mean()
		Psi_c[0,1:-1] = 0.25*(2.*Psi_pole + Psi[0,:-1]+Psi[0,1:])
		Psi_c[0,0]    = 0.25*(2.*Psi_pole + Psi[0,-1]+Psi[0,0])		
		Psi_c[0,-1]   = 0.25*(2.*Psi_pole + Psi[0,-1]+Psi[0,0])				

		# fix up low lat boundary
		# extrapolate linearly just like we did for the coordinates
		# (see genOutGrid in src/remix/mixio.F90)
		# note, neglecting the possibly non-uniform spacing (don't care)
		Psi_c[-1,:] = 2*Psi_c[-2,:]-Psi_c[-3,:]

		# now, do the differencing
		# for each cell corner on the original grid, I have the coordinates and Psi_c
		# need to find the gradient at cell center
		# the result is the same size as Psi

		# first etheta
		tmp    = 0.5*(Psi_c[:,1:]+Psi_c[:,:-1])  # move to edge center
		dPsi   = tmp[1:,:]-tmp[:-1,:]
		tmp    = 0.5*(theta[:,1:]+theta[:,:-1])
		dtheta = tmp[1:,:]-tmp[:-1,:]
		etheta = dPsi/dtheta/ri  # this is in V/m

		# now ephi
		tmp    = 0.5*(Psi_c[1:,:]+Psi_c[:-1,:])  # move to edge center
		dPsi   = tmp[:,1:]-tmp[:,:-1]
		tmp    = 0.5*(phi[1:,:]+phi[:-1,:])
		dphi   = tmp[:,1:]-tmp[:,:-1]
		tc = 0.25*(theta[:-1,:-1]+theta[1:,:-1]+theta[:-1,1:]+theta[1:,1:]) # need this additionally 
		ephi = dPsi/dphi/np.sin(tc)/ri  # this is in V/m

		if returnDeltas:
			return (-etheta,-ephi,dtheta,dphi)  # E = -grad Psi
		else:	
			return (-etheta,-ephi)  # E = -grad Psi			

	def joule(self):
		etheta,ephi = self.efield()
		SigmaP = self.variables['sigmap']['data']
		J = SigmaP*(etheta**2+ephi**2)  # this is in W/m^2
		return(J)

	# note, here we're taking the Cartesian average for cell centers
	# this is less accurate than the angular averages that are used in efield above
	# furthermore, I'm lazily mixing this with dtheta,dphi that came from efield in hCurrents
	# TODO: this should be fixed by pulling out the cell-centered tc,pc, and dtheta,dphi
	# from the efield calculation and using throughout.
	# note, however, that the staggered grid for storage is currently create in the remix Fortran code
	# by Cartesian averaging, so if we go to angular in this python postprocessing code, 
	# we should probably start with going angular in the Fortran code.
	def cartesianCellCenters(self):
		# Aliases to keep things short
		x = self.ion['X']
		y = self.ion['Y']

		xc = 0.25*(x[:-1,:-1]+x[1:,:-1]+x[:-1,1:]+x[1:,1:])
		yc = 0.25*(y[:-1,:-1]+y[1:,:-1]+y[:-1,1:]+y[1:,1:])		

		r,phi = self.get_spherical(xc,yc)

		theta = np.arcsin(r)

		return(xc,yc,theta,phi)
		
	# FIXME: MAKE WORK FOR SOUTH
	# TODO: write a separate geom function to define cell-centered coords, deltas, and even cosDipAngle (see comment above cartesianCellCenters)
	# right now, it's just ugly below where I take dtheta,phi from efield
	# and the coordinates separately from cartesianCellCenters.
	# it was a lazy solution
	def hCurrents(self): # horizontal currents
		etheta,ephi,dtheta,dphi = self.efield(returnDeltas=True)
		SigmaH = self.variables['sigmah']['data']
		SigmaP = self.variables['sigmap']['data']		

		xc,yc,theta,phi = self.cartesianCellCenters()

		cosDipAngle = -2.*np.cos(theta)/np.sqrt(1.+3.*np.cos(theta)**2)
		Jh_theta = -SigmaH*ephi/cosDipAngle
		Jh_phi   =  SigmaH*etheta/cosDipAngle
		Jp_theta =  SigmaP*etheta/cosDipAngle**2
		Jp_phi   =  SigmaP*ephi


		# current above is in SI units [A/m]
		# i.e., height-integrated current density
		return(xc,yc,theta,phi,dtheta,dphi,Jh_theta,Jh_phi,Jp_theta,Jp_phi,cosDipAngle)

	# Main function to compute magnetic perturbations using the Biot-Savart integration
	# This includes Hall, Pedersen and FAC with the option to do Hall only
	# Rin = Inner boundary of MHD grid [Re]
	def dB(self,xyz,hallOnly=True,Rin=2.0,rsegments=10):
		# xyz = array of points where to compute dB
		# xyz.shape should be (N,3), where N is the number of points
		# xyz = (x,y,z) in units of Ri

		if len(xyz.shape)!=2:
			sys.exit("dB input assumes the array of points of (N,3) size.")			
		if xyz.shape[1]!=3: 
			sys.exit("dB input assumes the array of points of (N,3) size.")

		self.init_vars('NORTH')
		x,y,theta,phi,dtheta,dphi,jht,jhp,jpt,jpp,cosDipAngle = self.hCurrents()
		z =  np.sqrt(1.-x**2-y**2)  # ASSUME NORTH
#		z = -np.sqrt(1.-x**2-y**2)	# ASSUME SOUTH

		# fake dimensions for numpy broadcasting
		xSource = x[:,:,np.newaxis]		
		ySource = y[:,:,np.newaxis]		
		zSource = z[:,:,np.newaxis]						

		tSource = theta[:,:,np.newaxis]
		pSource = phi[:,:,np.newaxis]
		dtSource = dtheta[:,:,np.newaxis]
		dpSource = dphi[:,:,np.newaxis]		

		# Hall
		jhTheta = jht[:,:,np.newaxis]
		jhPhi   = jhp[:,:,np.newaxis]
		# Pedersen
		jpTheta = jpt[:,:,np.newaxis]
		jpPhi   = jpp[:,:,np.newaxis]

		if not hallOnly:
			jTheta = jpTheta + jhTheta
			jPhi   = jpPhi   + jhPhi
		else:
			jTheta = jhTheta
			jPhi   = jhPhi

		# convert to Cartesian for the Biot-Savart summation
		# otherwise, spherical coordinates get mixed up betwen the source and destination grids
		# theta_unit and phi_unit vectors rotate from point to point and are different on the two grids
		jx = jTheta*np.cos(tSource)*np.cos(pSource) - jPhi*np.sin(pSource)
		jy = jTheta*np.cos(tSource)*np.sin(pSource) + jPhi*np.cos(pSource)
		jz =-jTheta*np.sin(tSource)

		# x,y,z are size (ntheta,nphi)
		# make array of destination points on the mix grid
		# add fake dimension for numpy broadcasting
		xDest = xyz[np.newaxis,np.newaxis,:,0]
		yDest = xyz[np.newaxis,np.newaxis,:,1]
		zDest = xyz[np.newaxis,np.newaxis,:,2]				

		# up to here things are fast (checked for both source and destination 90x720 grids)
		# the operations below are slow and kill the memory becase (90x720)^2 (the size of each array below) is ~30GB
		# solution: break the destination grid up into pieces before passing here

		# vector between destination and source
		Rx = xDest - xSource
		Ry = yDest - ySource
		Rz = zDest - zSource
		R  = np.sqrt(Rx**2+Ry**2+Rz**2)

		# vector product with the current
		# note the multiplication by sin(tSource)*dtSource*dpSource -- area of the surface element
		#
		# note on normalization
		# Efield is computed in V/m, sigma is also in SI units
		# Current coming out of hCurrents() should be in SI units [A/m] (height integrated).
		# The Biot-Savart law is mu0/4pi*Int( j x r' dV/|r'|^3 ) = 
		# mu0/4pi Int( j x r'hat r^2 dr sin(theta) dtheta dphi /|r'|^2) = 
		# mu0/4pi Int( (j dr) x r'hat (r/Ri)^2 sin(theta) dtheta dphi /(|r'|/Ri)^2)
		# where we have combined j and dr (= j dr) which is what is coming out of hCurrents in SI units
		# and normalized everything else to Ri, which it already is in the code below
		# in other words the fields below should be in [T]		
		dA = np.sin(tSource)*dtSource*dpSource
		dBx = mu0o4pi*np.sum( (jy*Rz - jz*Ry)*dA/R**3,axis=(0,1))
		dBy = mu0o4pi*np.sum( (jz*Rx - jx*Rz)*dA/R**3,axis=(0,1))
		dBz = mu0o4pi*np.sum( (jx*Ry - jy*Rx)*dA/R**3,axis=(0,1))

		if not hallOnly:
			intx,inty,intz = self.BSFluxTubeInt(xyz,Rinner=Rin*Re/Ri,rsegments=rsegments)
			# FIXME: don't fix sign for south
			jpara = -self.variables['current']['data'] # note, the sign was inverted by the reader for north, put it back to recover true FAC 
			jpara = jpara[:,:,np.newaxis]
			cosd  = abs(cosDipAngle[:,:,np.newaxis])  # note, only need abs value of cosd regardless of hemisphere

			# note on normalization
			# jpara is in microA/m^2 -- convert to A (1.e-6)
			# further, after all is said and done and all distance-like variables are accounted for
			# the answer below should be multiplied by Ri in m (6.5e6)
			# factors of 1.e6 cancel out and we only have Ri
			dBx += Ri*mu0o4pi*np.sum(jpara*dA*cosd*intx,axis=(0,1))
			dBy += Ri*mu0o4pi*np.sum(jpara*dA*cosd*inty,axis=(0,1))
			dBz += Ri*mu0o4pi*np.sum(jpara*dA*cosd*intz,axis=(0,1))


		# finally, convert to spherical *at the destination*
		# note, this is ugly because we specified the spherical grid before passing to this function (in calcdB.py)
		# FIXME: think about how to make it less ugly
		rDest = np.sqrt(xDest**2+yDest**2+zDest**2)
		tDest = np.arccos(zDest/rDest)
		pDest = np.arctan2(yDest,xDest)

		dBr     = dBx*np.sin(tDest)*np.cos(pDest) + dBy*np.sin(tDest)*np.sin(pDest) + dBz*np.cos(tDest)
		dBtheta = dBx*np.cos(tDest)*np.cos(pDest) + dBy*np.cos(tDest)*np.sin(pDest) - dBz*np.sin(tDest)	
		dBphi   =-dBx*np.sin(pDest) + dBy*np.cos(pDest)
		return(dBr,dBtheta,dBphi)

	# FIXME: Make work for SOUTH
	# Flux-tube Biot-Savart integral \int dl bhat x r'/|r'|^3
	def BSFluxTubeInt(self,xyz,Rinner,rsegments = 10):
		# xyz = array of points where to compute dB  (same as above in dB)
		# xyz.shape should be (N,3), where N is the number of points
		# xyz = (x,y,z) in units of Ri
		#
		# Rinner is the radius of the inner boundary of the MHD domain
		# expressed in Ri
		#
		# rsegments: how many segments to break the fluxtube into between ionosphere and Rinner

		if len(xyz.shape)!=2:
			sys.exit("dB input assumes the array of points of (N,3) size.")			
		if xyz.shape[1]!=3: 
			sys.exit("dB input assumes the array of points of (N,3) size.")

		xc,yc,theta,phi = self.cartesianCellCenters()

		# radii of centers of segments of the flux tube
		Rs = np.linspace(1.,Rinner,rsegments+1)
		Rcenters = 0.5*(Rs[:-1]+Rs[1:])
		dR = Rs[1:]-Rs[:-1]

		# add fake dims to conform to theta & phi size
		Rcenters = Rcenters[:,np.newaxis,np.newaxis]
		dR = dR[:,np.newaxis,np.newaxis]		

		# theta on the field line (m for magnetosphere)
		# note, since theta is in [0,45] deg range or so,
		# Rcenters is in the range [1,2] or so, and
		# arsin is in the range [-pi/2,pi/2]
		# the result below is >0, i.e., we only take the part of the field line
		# that lies in the same hemisphere as the ionospheric footpoint
		thetam = np.arcsin(np.sqrt(Rcenters)*np.sin(theta))   

		# note order: R, theta, phi
		x = Rcenters*np.sin(thetam)*np.cos(phi)
		y = Rcenters*np.sin(thetam)*np.sin(phi)
		z = np.sqrt(Rcenters**2-x**2-y**2)

		# distance along flux tube for given dR
#		dl = dR*np.sqrt(1.+Rcenters*np.sin(theta)**2/4./(1.-Rcenters*np.sin(theta)**2))
		dl = dR*np.sqrt(1.+3.*np.cos(thetam)**2)/2./np.cos(thetam)

		# unit vector in B-direction (e.g., Baumjohann Eq. 3.1 -- note their use of latitude vs colatitude)
#		br = -2.*np.sqrt(1-Rcenters*np.sin(theta)**2)/np.sqrt(4.-3.*Rcenters*np.sin(theta)**2)
#		bt = -np.sqrt(Rcenters)*np.sin(theta)/np.sqrt(4.-3.*Rcenters*np.sin(theta)**2)
		br = -2.*np.cos(thetam)/np.sqrt(np.sin(thetam)**2+4.*np.cos(thetam)**2)
		bt = -np.sin(thetam)/np.sqrt(np.sin(thetam)**2+4.*np.cos(thetam)**2)

		# now the same in cartesian
		bx = br*np.sin(thetam)*np.cos(phi) + bt*np.cos(thetam)*np.cos(phi)
		by = br*np.sin(thetam)*np.sin(phi) + bt*np.cos(thetam)*np.sin(phi)
		bz = br*np.cos(thetam) - bt*np.sin(thetam)

		# fake dimensions for numpy broadcasting
		# remember dimenstions: R,t,p along field line + adding the destination point number
		xSource = x[:,:,:,np.newaxis]
		ySource = y[:,:,:,np.newaxis]
		zSource = z[:,:,:,np.newaxis]

		bx = bx[:,:,:,np.newaxis]
		by = by[:,:,:,np.newaxis]
		bz = bz[:,:,:,np.newaxis]
		dl = dl[:,:,:,np.newaxis]

		xDest = xyz[np.newaxis,np.newaxis,np.newaxis,:,0]
		yDest = xyz[np.newaxis,np.newaxis,np.newaxis,:,1]
		zDest = xyz[np.newaxis,np.newaxis,np.newaxis,:,2]				

		# vector between destination and source
		Rx = xDest - xSource
		Ry = yDest - ySource
		Rz = zDest - zSource
		R  = np.sqrt(Rx**2+Ry**2+Rz**2)

		# vector product with the current
		intx = np.sum( dl*(by*Rz - bz*Ry)/R**3,axis=0)
		inty = np.sum( dl*(bz*Rx - bx*Rz)/R**3,axis=0)
		intz = np.sum( dl*(bx*Ry - by*Rx)/R**3,axis=0)

		return(intx,inty,intz)



# Code below is from old mix scripts.
# Keeping for further development but commenting out for now.

# The function is deprecated as of 10-25-10 and is left here only for
# backward compatibility. Use get_time_series for extracting cpcp and
# other quantities as funcitons of time.
# def get_cpcp_time_series(directory):
# 	from pylab import sort
# 	import glob,datetime

# 	files = glob.glob(directory+'*_mix_*.hdf')

# 	time   = []
# 	cpcp_n = []
# 	cpcp_s = []
# 	for file in sort(files):
# 		print(file)
# 		(simtime,
# 		 x,y,
# 		 psi_n,psi_s,
# 		 fac_n,fac_s,
# 		 sigmap_n,sigmap_s,
# 		 sigmah_n,sigmah_s,
# 		 energy_n,energy_s,
# 		 flux_n,flux_s) = get_data(file)
	
# 		cpcp_n.append(psi_n.max()-psi_n.min())
# 		cpcp_s.append(psi_s.max()-psi_s.min())
		
# 		t = datetime.datetime(simtime[0],simtime[1],simtime[2],
# 							  simtime[3],simtime[4],simtime[5])
# 		time.append(t)


# 	return (time,cpcp_n,cpcp_s)

# def get_time_series(directory):
# 	from pylab import sort,sqrt,nonzero
# 	import glob,datetime
# 	import integrate

# 	files = glob.glob(directory+'*_mix_*.hdf')


# 	ri=6500.e3
# 	# Pre-compute some constants:
# 	(simtime,x,y,psi_n,psi_s,fac_n,fac_s,sigmap_n,sigmap_s,sigmah_n,sigmah_s,
# 	 energy_n,energy_s,flux_n,flux_s) = get_data(files[0])
# 	x[:,0]=0.0
# 	y[:,0]=0.0
# 	z=sqrt(1.0-x**2-y**2)
# 	areaSphere = integrate.calcFaceAreas(x,y,z)*ri*ri

# 	time   = []
# 	cpcp_n = []
# 	cpcp_s = []
# 	fac_pos_n = []
# 	fac_neg_n = []
# 	fac_pos_s = []
# 	fac_neg_s = []

# 	for file in sort(files):
# 		(simtime,
# 		 x,y,
# 		 psi_n,psi_s,
# 		 fac_n,fac_s,
# 		 sigmap_n,sigmap_s,
# 		 sigmah_n,sigmah_s,
# 		 energy_n,energy_s,
# 		 flux_n,flux_s) = get_data(file)
	
# 		cpcp_n.append(psi_n.max()-psi_n.min())
# 		cpcp_s.append(psi_s.max()-psi_s.min())

# 		pfac = fac_n.copy()
# 		locs = nonzero(fac_n < 0)
# 		pfac[locs] = 0.0
# 		val = areaSphere*pfac[0:180,0:26]
# 		fac_pos_n.append(val.sum()*1.e-12)

# 		pfac = fac_n.copy()
# 		locs = nonzero(fac_n > 0)
# 		pfac[locs] = 0.0
# 		val = areaSphere*pfac[0:180,0:26]
# 		fac_neg_n.append(-val.sum()*1.e-12)


# 		pfac = fac_s.copy()
# 		locs = nonzero(fac_s < 0)
# 		pfac[locs] = 0.0
# 		val = areaSphere*pfac[0:180,0:26]
# 		fac_pos_s.append(val.sum()*1.e-12)

# 		pfac = fac_s.copy()
# 		locs = nonzero(fac_s > 0)
# 		pfac[locs] = 0.0
# 		val = areaSphere*pfac[0:180,0:26]
# 		fac_neg_s.append(-val.sum()*1.e-12)

		
# 		t = datetime.datetime(simtime[0],simtime[1],simtime[2],
# 							  simtime[3],simtime[4],simtime[5])
# 		time.append(t)

# 	return (time,cpcp_n,cpcp_s,fac_pos_n,fac_neg_n,fac_pos_s,fac_neg_s)






