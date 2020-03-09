import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys

class remix:

	def __init__(self,h5file,step):
		# create the ion object to store data and coordinates
		self.ion = self.get_data(h5file,step)
		self.Initialized=False

		# DEFINE DATA LIMITS
		self.variables = { 'potential' : {'min':-100,
										'max': 100},
						 'current'   : {'min':-1,
										'max':1},
						 'sigmap'    : {'min':1,
										'max':10},
						 'sigmah'    : {'min':2,
										'max':20},
						 'energy'    : {'min':0,
										'max':10},
						 'flux'      : {'min':0,
										'max':1.e9},
						 'eflux'     : {'min':0,
										'max':1.},
						 'efield'    : {'min':-1,
										'max':1},
						 'joule'     : {'min':0,
										'max':10},
						 'jhall'	 : {'min':-2,
						 				'max':2}
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
		else:  # note flipping the y(phi)-axis
			self.variables['potential']['data'] = self.ion['Potential '+h][:,::-1]
			self.variables['current']['data']   = self.ion['Field-aligned current '+h][:,::-1]
			self.variables['sigmap']['data']    = self.ion['Pedersen conductance '+h][:,::-1]
			self.variables['sigmah']['data']    = self.ion['Hall conductance '+h][:,::-1]
			self.variables['energy']['data']    = self.ion['Average energy '+h][:,::-1]
			self.variables['flux']['data']      = self.ion['Number flux '+h][:,::-1]

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

	# TODO: define and consolidate allowed variable names
	def plot(self,varname,
			 ncontours=16,   # default number of potential contours
			 addlabels={},
			 gs=None):

		# define function for potential contour overplotting
		# to keep code below clean and compact
		def potential_overplot():
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
			ax.contour(tc+np.pi/2.,rc,tmp,15,colors='black',linewidths=0.5)

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
		if varname !='flux': 
			if varname == 'jped':
				format_str = '%.2f'
			else:
				format_str = '%.1f'
		else:
			format_str = '%.1e'

		# define red/blue colortable for potential and current
		latlblclr = 'black'
		if (varname == 'potential') or (varname == 'current'):
			cmap=cm.RdBu_r
		elif varname in ['flux','energy','eflux','joule']:
			cmap=cm.inferno			
			latlblclr = 'white'
		elif (varname == 'velocity'): 
			cmap=cm.YlOrRd
		elif (varname == 'efield'): 
			cmap=None #cm.RdBu_r#None # cm.hsv
		elif (varname == 'magpert'): 
			cmap = cm.YlOrRd
		else:
			cmap=None # default is used
		
		# DEFINE GRID LINES AND LABELS
		circle_list = [10,20,30,40]
		circles = np.sin(np.array(circle_list)*np.pi/180.)

		# convert to string and add degree symbol
		lbls = [str(elem)+u'\xb0' for elem in circle_list] 
		
		hour_labels = ['06','12','18','00']

		if varname == 'joule':
			self.variables[varname]['data'] = self.joule()*1.e3  # convert to mW/m^2

		variable = self.variables[varname]['data']

		fig = plt.gcf()
		# Now plotting
		if gs != None:
			ax=fig.add_subplot(gs,polar=True)
		else:
			ax=fig.add_subplot(polar=True) 

		p=ax.pcolormesh(theta+np.pi/2.,r,variable,cmap=cmap,vmin=lower,vmax=upper)
		cb=plt.colorbar(p,ax=ax,pad=0.1,shrink=0.85)  
		cb.set_label(cblabels[varname])

		lines, labels = plt.rgrids(circles,lbls,fontsize=8,color=latlblclr)
		lines, labels = plt.thetagrids((0.,90.,180.,270.),hour_labels)
		ax.axis([0,2*np.pi,0,r.max()],'tight')
		ax.text(-75.*np.pi/180.,1.2*r.max(),('min: '+format_str+'\nmax: ' +format_str) % 
			  (variable.min() ,variable.max()))
		ax.grid(True)

		if varname=='current': 
			potential_overplot()

	# mpl.rcParams['contour.negative_linestyle'] = 'solid'
	# if (varname == 'efield' or varname == 'velocity' or varname =='joule'): 
	#     contour(theta+pi/2.,r,variables['potential']['data'][:,2:-1],21,colors='black')
#                                      arange(variables['potential']['min'],variables['potential']['max'],21.),colors='purple')

	# FIXME: MAKE WORK FOR SOUTH (I THINK IT DOES BUT MAKE SURE)
	def efield(self,ri=6.5e3): 
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

		return (-etheta,-ephi)  # E = -grad Psi

	def joule(self):
		etheta,ephi = self.efield()
		SigmaP = self.variables['sigmap']['data']
		J = SigmaP*(etheta**2+ephi**2)  # this is in W/m^2
		return(J)

	# FIXME: MAKE WORK FOR SOUTH
	def jHall(self):
		etheta,ephi = self.efield()
		SigmaH = self.variables['sigmah']['data']

		# Aliases to keep things short
		x = self.ion['X']
		y = self.ion['Y']

		xc = 0.25*(x[:-1,:-1]+x[1:,:-1]+x[:-1,1:]+x[1:,1:])
		yc = 0.25*(y[:-1,:-1]+y[1:,:-1]+y[:-1,1:]+y[1:,1:])		

		r,phi = self.get_spherical(xc,yc)

		theta = np.arcsin(r)

		cosDipAngle = -2.*np.cos(theta)/np.sqrt(1.+3.*np.cos(theta)**2)
		Jh_theta = -SigmaH*ephi/cosDipAngle
		Jh_phi   =  SigmaH*etheta/cosDipAngle

		return(xc,yc,theta,phi,Jh_theta,Jh_phi)

	def dB(self,xyz):
		# xyz = array of points where to compute dB
		# xyz.shape should be (N,3), where N is the number of points
		# xyz = (x,y,z) in units of Ri

		mu4pi = 1. # FIXME: change to real values

		if len(xyz.shape)!=2:
			sys.exit("dB input assumes the array of points of (N,3) size.")			
		if xyz.shape[1]!=3: 
			sys.exit("dB input assumes the array of points of (N,3) size.")

		nPoints = xyz.shape[0]

		self.init_vars('NORTH')
		x,y,theta,phi,jht,jhp = self.jHall()
		z =  np.sqrt(1.-x**2-y**2)  # ASSUME NORTH
#		z = -np.sqrt(1.-x**2-y**2)	# ASSUME SOUTH

		# fake dimensions for numpy broadcasting
		xSource = x[:,:,np.newaxis]		
		ySource = y[:,:,np.newaxis]		
		zSource = z[:,:,np.newaxis]						

		tSource = theta[:,:,np.newaxis]
		pSource = phi[:,:,np.newaxis]
		jhTheta = jht[:,:,np.newaxis]
		jhPhi   = jhp[:,:,np.newaxis]		

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

		# convert R to spherical to compute the vector product with the current
		# since the current is in spherical and we want output in spherical
		Rr     = Rx*np.sin(tSource)*np.cos(pSource) + Ry*np.sin(tSource)*np.sin(pSource) + Rz*np.cos(tSource)
		Rtheta = Rx*np.cos(tSource)*np.cos(pSource) + Ry*np.cos(tSource)*np.sin(pSource) - Rz*np.sin(tSource)	
		Rphi   =-Rx*np.sin(pSource) + Ry*np.cos(pSource)
		
		# vector product with the current
		# note the multiplication by sin(tSource) -- it's the area of the surface element
		dBphi   = np.sum(-jhTheta*np.sin(tSource)*Rr/R**3,axis=(0,1))
		dBtheta = np.sum(jhPhi*np.sin(tSource)*Rr/R**3   ,axis=(0,1))
		dBr     = np.sum( (jhTheta*Rphi - jhPhi*Rtheta)*np.sin(tSource)/R**3,axis=(0,1))

		# FIXME: change to real values
		dtheta = 1.
		dphi   = 1. 

		return(dBr,dBtheta,dBphi)



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






