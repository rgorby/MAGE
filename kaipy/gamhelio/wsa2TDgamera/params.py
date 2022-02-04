import configparser

class params():
    def __init__(self,ConfigFileName):
        config = configparser.ConfigParser()
        config.read(ConfigFileName)
        

        #self.initLFMfile = config.get('Output','Prefix')+'_%dx%dx%d_mhd_0000000.hdf'%(self.ni,self.nj,self.nk)
    
        self.scale    = config.getfloat('Normalization','Xscale')
        self.gameraGridFile = config.get('Grid','gameraGridFile')
        self.GridDir  = config.get('Grid','GridDir')
        self.gameraIbcFile = config.get('Grid','gameraIbcFile')
        self.IbcDir   = config.get('Grid','IbcDir')
        self.tMin     = config.getfloat('Grid','tMin')
        self.tMax     = config.getfloat('Grid','tMax')
        self.Rin      = config.getfloat('Grid','Rin')
        self.Rout     = config.getfloat('Grid','Rout')
        self.Ni       = config.getint('Grid','Ni')
        self.Nj       = config.getint('Grid','Nj')
        self.Nk       = config.getint('Grid','Nk')

        
        self.gamma = config.getfloat('Constants','gamma')
        self.NO2   = config.getint('Constants','NO2')
        self.Tsolar = config.getfloat('Constants','Tsolar')
        
        if config.has_section('MAS'):
            self.masdir = config.get('MAS','masdir')
            self.masTimeLabel = config.get('MAS','masTimeLabel')
            self.masFrame = config.get('MAS','masFrame')
            self.masFakeRotation = config.getboolean('MAS','masFakeRotation')

        if config.has_section('WSA'):
            self.wsaFile = config.get('WSA','wsafile')
            self.gaussSmoothWidth = config.getint('WSA','gauss_smooth_width')
            self.plots = config.getboolean('WSA','plots')
            self.densTempInfile = config.getboolean('WSA','density_temperature_infile')
            self.normalized = config.getboolean('WSA','normalized')

        if config.has_section('ADAPT'):
            self.adaptdir = config.get('ADAPT','adaptdir')
            self.adaptWildCard = config.get('ADAPT','wild_card')
            self.adaptCadence = config.getfloat('ADAPT','cadence')
            self.gaussSmoothWidth = config.getint('ADAPT','gauss_smooth_width')
            self.plots = config.getboolean('ADAPT','plots')
            self.densTempInfile = config.getboolean('ADAPT','density_temperature_infile')
            self.normalized = config.getboolean('ADAPT','normalized')


        self.dumpInit = config.getboolean('DUMPS','init')
        self.dumpBC   = config.getboolean('DUMPS','BC')


