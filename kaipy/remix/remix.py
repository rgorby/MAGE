import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import sys

class remix:

	def __init__(self,h5file,step):
		# create the ion object to store data and coordinates
		self.ion = self.get_data(h5file,step)
		self.Initialized=True

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
										'max':100},
						 'flux'      : {'min':0,
										'max':3.e8},
						 'efield'    : {'min':-1,
										'max':1},
						  'joule'    : {'min':0,
										'max':50},
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
			 ncontours=51,
			 nticks=11,
			 addlabels={},
			 ax=None):

		# Aliases to keep things short
		x = self.ion['X']
		y = self.ion['Y']
		r = self.ion['R']
		theta = self.ion['THETA']

		# List all possible variable names here, add more from the input parameter if needed
		cblabels = {'potential' : r'Potential [kV]',
					'current'   : r'Current density [$\mu\mathrm{A/m^2}$]',
					'sigmap'    : r'Pedersen conductance [S]',
					'sigmah'    : r'Hall conductance [S]',
					'energy'    : r'Energy [keV]',
					'flux'      : r'Flux [$\mathrm{1/cm^2s}$]',
					'ephi'      : r'$E_\phi$ [mV/m]',
					'etheta'    : r'$E_\theta$ [mV/m]',
					'efield'    : r'|E| [mV/m]',
					'joule'     : r'Joule heating [mW/m$^2$]',
					'jped'      : r'Pedersen current [$\mu\mathrm{A/m^2}$]',
					'magpert'   : r'Magnetic perturbation [nT]',
					}
		cblabels.update(addlabels)

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
		if varname != 'flux': 
			if varname == 'jped':
				format_str = '%.2f'
			else:
				format_str = '%.1f'
		else:
			format_str = '%.1e'

		# define red/blue colortable for potential and current
		if (varname == 'potential') or (varname == 'current'):
			cmap=cm.RdBu_r
		elif (varname == 'flux') or (varname == 'energy'): 
			cmap=cm.inferno			
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

		contours = np.linspace(lower,upper,ncontours)
		ticks = np.linspace(lower,upper,nticks)

		variable = self.variables[varname]['data']

		p=ax.pcolormesh(theta+np.pi/2.,r,variable,cmap=cmap,vmin=lower,vmax=upper)
		cb=plt.colorbar(p,ax=ax,pad=0.1,ticks=ticks,shrink=0.85)
		cb.set_label(cblabels[varname])

		lines, labels = plt.rgrids(circles,lbls,fontsize=8)
		lines, labels = plt.thetagrids((0.,90.,180.,270.),hour_labels)
		ax.axis([0,2*np.pi,0,r.max()],'tight')
		ax.text(-75.*np.pi/180.,1.2*r.max(),('min: '+format_str+'\nmax: ' +format_str) % 
			  (variable.min() ,variable.max()))
		ax.grid(True)

		if varname=='current': 
			# use this trick to define contour levels.
			# if you just use number of levels as an argument for contour,
			# matplotlib plots a bad contour near pole (apparently a bug)
			potmax = np.max(abs(self.variables['potential']['data']))
			levels = np.linspace(-potmax,potmax,16)

			tc = 0.25*(theta[:-1,:-1]+theta[1:,:-1]+theta[:-1,1:]+theta[1:,1:])
			rc = 0.25*(r[:-1,:-1]+r[1:,:-1]+r[:-1,1:]+r[1:,1:])		
			ax.contour(tc+np.pi/2.,rc,self.variables['potential']['data'],levels=levels,colors='black',linewidths=0.5)

			# this is another trick  to avoid problems with contours across periodic boundary (noon)
			# above seems to work
#			xc = 0.25*(x[:-1,:-1]+x[1:,:-1]+x[:-1,1:]+x[1:,1:])
#			yc = 0.25*(y[:-1,:-1]+y[1:,:-1]+y[:-1,1:]+y[1:,1:])		
#			fig = plt.gcf()
#			new_axis = fig.add_axes(ax.get_position(), frameon = False)
	#		new_axis.plot()
#			new_axis.set_axis_off()
#			new_axis.contour(-yc,xc,self.variables['potential']['data'],levels=levels,colors='black',linewidths=0.5)
	
	# mpl.rcParams['contour.negative_linestyle'] = 'solid'
	# if (varname == 'efield' or varname == 'velocity' or varname =='joule'): 
	#     contour(theta+pi/2.,r,variables['potential']['data'][:,2:-1],21,colors='black')
#                                      arange(variables['potential']['min'],variables['potential']['max'],21.),colors='purple')

	def efield(self): 
		if not self.Initialized:
			sys.exit("Variables should be initialized for the specific hemisphere (call init_var) prior to efield calculation.")

		Nr,Nt = self.variables['potential']['data'].shape  # note, these are numbers of cells. self.ion['X'].shape = Nr+1,Nt+1

		# Aliases to keep things short
		x = self.ion['X']
		y = self.ion['Y']
		r = self.ion['R']
		theta = self.ion['THETA']

		# dtheta=theta[0,1:]-theta[0,:-1]
		# dphi=phi[1:,1]-phi[:-1,1]


		# get deltas

		# phi=p.arctan2(y,x)
		# phi[phi<0]=phi[phi<0]+2*p.pi
	   
		# # Fix phi at pole just in case
		# phi[:,0]=phi[:,1]

		# theta=p.arcsin(p.sqrt(x**2+y**2))


		# ephi   = p.zeros(x.shape)
		# etheta = p.zeros(x.shape)

		# phi_mid   = p.zeros(x.shape)
		# theta_mid = p.zeros(x.shape)

		# nphi   = x.shape[0]
		# ntheta = x.shape[1]

		# for j in range(nphi-1): # loop over phi
		# 	theta_mid[j,:-1] = (theta[j,1:]+theta[j,:-1])/2.
		# 	etheta[j,:-1] =( (psi[j,1:]-psi[j,:-1])+(psi[j+1,1:]-psi[j+1,:-1]))/2./dtheta/ri*1.e6

		# for i in range(ntheta-1): # loop over theta
		# 	phi_mid[:-1,i] = ( phi[1:,i]+phi[:-1,i] )/2.
		# 	ephi[:-1,i] = 1./p.sin(theta_mid[0,i])*( (psi[1:,i]-psi[:-1,i])+(psi[1:,i+1]-psi[:-1,i+1]) )/2./dphi/ri*1.e6

		# # fixup the periodic axis
		# phi_mid[nphi-1,:]   = phi_mid[0,:]+2.*p.pi
		# theta_mid[nphi-1,:]   = theta_mid[0,:]
		# ephi[nphi-1,:] = ephi[0,:]
		# etheta[nphi-1,:] = etheta[0,:]

		# # return the angular coordinates at cell centers (first tuple) and
		# # electric field in mV/m (second tuple)
		# return (phi_mid,theta_mid),(ephi,etheta)



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




# def joule(efield,sigmap):
# 	from pylab import zeros
# 	# Average conductance
# 	sigp_mid = zeros(efield[0].shape)
# 	sigp_mid[:-1,:-1] = 0.25*(sigmap[:-1,:-1]+sigmap[:-1,1:]+sigmap[1:,:-1]+sigmap[1:,1:])

# 	# fixup the periodic axis
# 	sigp_mid[-1,:]   = sigp_mid[0,:]

# 	return (sigp_mid*(efield[0]**2+efield[1]**2)*1.e-3)  # Joule heating in mW/m^2

# def jped(efield,sigmap):
# 	from pylab import zeros,sqrt
# 	# Average conductance
# 	sigp_mid = zeros(efield[0].shape)
# 	sigp_mid[:-1,:-1] = 0.25*(sigmap[:-1,:-1]+sigmap[:-1,1:]+sigmap[1:,:-1]+sigmap[1:,1:])

# 	# fixup the periodic axis
# 	sigp_mid[-1,:]   = sigp_mid[0,:]

# 	return (sigp_mid*sqrt(efield[0]**2+efield[1]**2)*1.e-3)  # should be in microA/m


