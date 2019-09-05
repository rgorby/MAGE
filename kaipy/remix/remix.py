import h5py

def get_data(h5file):
    import numpy
    from numpy import array

    f = h5py.File(h5file,'r')

    ion = {}
    # Add simtime later
    # simtime = getattr(hdffile,'UT (Year, Month, Day, Hour, Min, Sec)')
    ion['Sim time'] = array([1900,1,1,0,0,0])
    for h in f.keys():
        ion[h] = {}
        ion[h]['x'] = add_periodic(f[h]['X'][:])
        ion[h]['y'] = add_periodic(f[h]['Y'][:])
        ion[h]['psi'] = add_periodic(f[h]['Potential'][:])
        ion[h]['fac'] = add_periodic(f[h]['Field-aligned current'][:])
        ion[h]['sigmap'] = add_periodic(f[h]['Pedersen conductance'][:])
        ion[h]['sigmah'] = add_periodic(f[h]['Hall conductance'][:])
        ion[h]['energy'] = add_periodic(f[h]['Average energy'][:])
        ion[h]['flux']   = add_periodic(f[h]['Number flux'][:])

    f.close()

    return ion


def add_periodic(A):
    from numpy import hstack
    return(hstack((A,A[:,[0]])))

def plot(theta,r,variables,varname,
         ncontours=51,
         nticks=11,
         addlabels={}):
    from pylab import pi as pi
    from pylab import cm
    from pylab import sin,contour,contourf,colorbar,rgrids,thetagrids,axis,text,linspace,arange,pcolormesh

    # List all possible variable names here, add more from the input parameter if needed
    cblabels = {'potential' : r'Potential, kV',
                'current'   : r'Current density, $\mu\mathrm{A/m^2}$',
                'sigmap'    : r'Pedersen conductance, S',
                'sigmah'    : r'Hall conductance, S',
                'energy'    : r'Energy, keV',
                'flux'      : r'Flux, $\mathrm{1/cm^2s}$',
                'ephi'      : r'$E_\phi$,  mV/m',
                'etheta'    : r'$E_\theta$, mV/m',
                'efield'    : r'|E|, mV/m',
                'joule'     : r'Joule heating, mW/m$^2$',
                'jped'      : r'Pedersen current, $\mu\mathrm{A/m^2}$',
                'magpert'   : r'Magnetic perturbation, nT',
                }
    cblabels.update(addlabels)

    # if limits are given use them, if not use the variables min/max values
    if variables[varname].has_key('min'):
        lower = variables[varname]['min']
    else:
        lower = variables[varname]['data'].min()

    if variables[varname].has_key('max'):
        upper = variables[varname]['max']
    else:
        upper = variables[varname]['data'].max()


# Leaving this here for possible later use of pcolor instead of contours
    # pcolor(-y,x,psi_n*1.e-3) 
    # cb=colorbar()
    # cb.set_label('Potential, kV')
    # pcolor(-y,x,psi_n*1.e-3) 
    # axes(polar=True)

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
    elif (varname == 'velocity'): 
        cmap=cm.YlOrRd
    elif (varname == 'efield'): 
        cmap=None #cm.RdBu_r#None # cm.hsv
    elif (varname == 'magpert'): 
        cmap = cm.YlOrRd
    else:
        cmap=None # default is used
    
# DEFINE GRID LINES AND LABELS
# it's a list of circles in degrees
    circle_list = [10,20,30,40] 
# circles in radians
    circles = [sin(elem*pi/180.) for elem in circle_list] 
# convert to string and add degree symbol
    lbls = [str(elem)+u'\xb0' for elem in circle_list] 
    # latex string for the same
    # ['10$^\circ$','20$^\circ$','30$^\circ$','40$^\circ$']

    hour_labels = ['06','12','18','00']#{'north':['06','12','18','00'],'south':['18','12','06','00']}

    contours = linspace(lower,upper,ncontours)
    ticks = linspace(lower,upper,nticks)

    variable = variables[varname]['data']
#    contourf(theta+pi/2.,r,variable,contours,extend='both',cmap=cmap)
    pcolormesh(theta+pi/2.,r,variable)
    cb=colorbar(pad=0.075,ticks=ticks,shrink=0.85)
    cb.set_label(cblabels[varname])
#    contour(theta+pi/2.,r,variable,contours,cmap=cmap)
    
#    from pylab import asarray 
    if (varname == 'potential'): contour(theta+pi/2.,r,variables['potential']['data'],21,colors='black')
    import matplotlib
#    matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
    # if (varname == 'efield' or varname == 'velocity' or varname =='joule'): 
    #     contour(theta+pi/2.,r,variables['potential']['data'][:,2:-1],21,colors='black')
#                                      arange(variables['potential']['min'],variables['potential']['max'],21.),colors='purple')

#    if (varname == 'current'): 
#        contour(theta+pi/2.,r,variables['potential']['data'],21,colors='black')

    lines, labels = rgrids(circles,lbls)
    lines, labels = thetagrids((0.,90.,180.,270.),hour_labels)
    axis([0,2*pi,0,r.max()],'tight')
    text(-75.*pi/180.,1.2*r.max(),('min: '+format_str+'\nmax: ' +format_str) % 
          (variable.min() ,variable.max()))

    

# The function is deprecated as of 10-25-10 and is left here only for
# backward compatibility. Use get_time_series for extracting cpcp and
# other quantities as funcitons of time.
def get_cpcp_time_series(directory):
    from pylab import sort
    import glob,datetime

    files = glob.glob(directory+'*_mix_*.hdf')

    time   = []
    cpcp_n = []
    cpcp_s = []
    for file in sort(files):
        print(file)
        (simtime,
         x,y,
         psi_n,psi_s,
         fac_n,fac_s,
         sigmap_n,sigmap_s,
         sigmah_n,sigmah_s,
         energy_n,energy_s,
         flux_n,flux_s) = get_data(file)
    
        cpcp_n.append(psi_n.max()-psi_n.min())
        cpcp_s.append(psi_s.max()-psi_s.min())
        
        t = datetime.datetime(simtime[0],simtime[1],simtime[2],
                              simtime[3],simtime[4],simtime[5])
        time.append(t)


    return (time,cpcp_n,cpcp_s)

def get_time_series(directory):
    from pylab import sort,sqrt,nonzero
    import glob,datetime
    import integrate

    files = glob.glob(directory+'*_mix_*.hdf')


    ri=6500.e3
    # Pre-compute some constants:
    (simtime,x,y,psi_n,psi_s,fac_n,fac_s,sigmap_n,sigmap_s,sigmah_n,sigmah_s,
     energy_n,energy_s,flux_n,flux_s) = get_data(files[0])
    x[:,0]=0.0
    y[:,0]=0.0
    z=sqrt(1.0-x**2-y**2)
    areaSphere = integrate.calcFaceAreas(x,y,z)*ri*ri

    time   = []
    cpcp_n = []
    cpcp_s = []
    fac_pos_n = []
    fac_neg_n = []
    fac_pos_s = []
    fac_neg_s = []

    for file in sort(files):
        (simtime,
         x,y,
         psi_n,psi_s,
         fac_n,fac_s,
         sigmap_n,sigmap_s,
         sigmah_n,sigmah_s,
         energy_n,energy_s,
         flux_n,flux_s) = get_data(file)
    
        cpcp_n.append(psi_n.max()-psi_n.min())
        cpcp_s.append(psi_s.max()-psi_s.min())

        pfac = fac_n.copy()
        locs = nonzero(fac_n < 0)
        pfac[locs] = 0.0
        val = areaSphere*pfac[0:180,0:26]
        fac_pos_n.append(val.sum()*1.e-12)

        pfac = fac_n.copy()
        locs = nonzero(fac_n > 0)
        pfac[locs] = 0.0
        val = areaSphere*pfac[0:180,0:26]
        fac_neg_n.append(-val.sum()*1.e-12)


        pfac = fac_s.copy()
        locs = nonzero(fac_s < 0)
        pfac[locs] = 0.0
        val = areaSphere*pfac[0:180,0:26]
        fac_pos_s.append(val.sum()*1.e-12)

        pfac = fac_s.copy()
        locs = nonzero(fac_s > 0)
        pfac[locs] = 0.0
        val = areaSphere*pfac[0:180,0:26]
        fac_neg_s.append(-val.sum()*1.e-12)

        
        t = datetime.datetime(simtime[0],simtime[1],simtime[2],
                              simtime[3],simtime[4],simtime[5])
        time.append(t)

    return (time,cpcp_n,cpcp_s,fac_pos_n,fac_neg_n,fac_pos_s,fac_neg_s)



def efield(x,y,psi,ri): #Radius of ionosphere in m
    import pylab as p

    phi=p.arctan2(y,x)
    phi[phi<0]=phi[phi<0]+2*p.pi
   
    # Fix phi at pole just in case
    phi[:,0]=phi[:,1]

    theta=p.arcsin(p.sqrt(x**2+y**2))

    dtheta=theta[0,1:]-theta[0,:-1]
    dphi=phi[1:,1]-phi[:-1,1]

    ephi   = p.zeros(x.shape)
    etheta = p.zeros(x.shape)

    phi_mid   = p.zeros(x.shape)
    theta_mid = p.zeros(x.shape)

    nphi   = x.shape[0]
    ntheta = x.shape[1]

    for j in range(nphi-1): # loop over phi
        theta_mid[j,:-1] = (theta[j,1:]+theta[j,:-1])/2.
        etheta[j,:-1] =( (psi[j,1:]-psi[j,:-1])+(psi[j+1,1:]-psi[j+1,:-1]))/2./dtheta/ri*1.e6

    for i in range(ntheta-1): # loop over theta
        phi_mid[:-1,i] = ( phi[1:,i]+phi[:-1,i] )/2.
        ephi[:-1,i] = 1./p.sin(theta_mid[0,i])*( (psi[1:,i]-psi[:-1,i])+(psi[1:,i+1]-psi[:-1,i+1]) )/2./dphi/ri*1.e6

    # fixup the periodic axis
    phi_mid[nphi-1,:]   = phi_mid[0,:]+2.*p.pi
    theta_mid[nphi-1,:]   = theta_mid[0,:]
    ephi[nphi-1,:] = ephi[0,:]
    etheta[nphi-1,:] = etheta[0,:]

    # return the angular coordinates at cell centers (first tuple) and
    # electric field in mV/m (second tuple)
    return (phi_mid,theta_mid),(ephi,etheta)

def joule(efield,sigmap):
    from pylab import zeros
    # Average conductance
    sigp_mid = zeros(efield[0].shape)
    sigp_mid[:-1,:-1] = 0.25*(sigmap[:-1,:-1]+sigmap[:-1,1:]+sigmap[1:,:-1]+sigmap[1:,1:])

    # fixup the periodic axis
    sigp_mid[-1,:]   = sigp_mid[0,:]

    return (sigp_mid*(efield[0]**2+efield[1]**2)*1.e-3)  # Joule heating in mW/m^2

def jped(efield,sigmap):
    from pylab import zeros,sqrt
    # Average conductance
    sigp_mid = zeros(efield[0].shape)
    sigp_mid[:-1,:-1] = 0.25*(sigmap[:-1,:-1]+sigmap[:-1,1:]+sigmap[1:,:-1]+sigmap[1:,1:])

    # fixup the periodic axis
    sigp_mid[-1,:]   = sigp_mid[0,:]

    return (sigp_mid*sqrt(efield[0]**2+efield[1]**2)*1.e-3)  # should be in microA/m


