#Main constants
import numpy as np

#------
#Helpful conversions
#------
G2nT  = 1E5              # Gauss->nanoTesla
kev2J = 1.602176634E-16  # keV -> J
ev2J  = 1.602176634E-19  # eV  -> J
erg2J = 1e-7             # erg -> J



#------
#Physical Constants
#------
Mu0     = 4E-7*np.pi             # [Tm/A]
Me_cgs  = 9.1093837015E-28       # [g]  Electron mass
Mp_cgs  = 1.67262192369E-24      # [g]  Proton mass
eCharge = 1.602E-19              # [J]  Charge of electron
dalton  = 1.66053906660*1.0E-27  # [kg] Mass unit



#------
#Planetary Constants
#------
Re_cgs = 6.3781E8       # [cm]  Earth's radius
EarthM0g = 0.2961737    # Gauss, Olsen++ 2000
REarth = Re_cgs*1.0e-2  # m
RionE  = 6.5            # Earth Ionosphere radius in 1000 km
EarthPsi0 = 92.4        # Corotation potential [kV]
#Saturn
SaturnM0g = 0.21        # Gauss
RSaturnXE = 9.5         # Rx = X*Re
#Jupiter
JupiterM0g = 4.8        # Gauss
RJupiterXE = 11.0       # !Rx = X*Re
#Mercury
MercuryM0g = 0.00345    # Gauss
RMercuryXE = 0.31397    # Rx = X*Re
#Neptune
NeptuneM0g = 0.142      # Gauss
RNeptuneXE = 3.860      # Rx = X*Re



#------
#Helio
#------
Rsolar = 6.956E5  #[km] Solar radius
kbltz  = 1.38e-16 #Boltzmann constant [erg/K]
mp     = 1.67e-24 #Proton mass in grams
Tsolar = 25.38    #Siderial solar rotation period
JD2MJD = 2400000.5 #Conversion from JD to MJD: MJD = JD - 2400000.5

#------
#Output defaults
#------
barLen = 30
barLab = 30 