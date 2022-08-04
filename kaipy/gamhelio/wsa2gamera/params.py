import configparser

class params():
    def __init__(self,ConfigFileName):
        config = configparser.ConfigParser(inline_comment_prefixes=(';','#'))
        config.read(ConfigFileName)
        
        self.gameraGridFile = config['Gamera']['gameraGridFile']
        self.GridDir = config['Gamera']['GridDir']
        self.gameraIbcFile = config['Gamera']['gameraIbcFile'] 
        self.IbcDir = config['Gamera']['IbcDir']   

        self.wsaFile = config['WSA']['wsafile']
        self.gaussSmoothWidth = config.getint('WSA','gauss_smooth_width')
        #self.plots = config.getboolean('WSA','plots')
        self.densTempInfile = config.getboolean('WSA','density_temperature_infile')
        self.normalized = config.getboolean('WSA','normalized')

        self.gamma = config.getfloat('Constants','gamma')
        self.Nghost   = config.getint('Constants','Nghost')
        self.Tsolar = config.getfloat('Constants','Tsolar')
        self.TCS = config.getfloat('Constants','TCS')
        self.nCS = config.getfloat('Constants','nCS')

        self.B0 = config.getfloat('Normalization','B0')
        self.n0 = config.getfloat('Normalization','n0')
        # self.T0 = config.getfloat('Normalization','T0')
        
        
        self.tMin     = config.getfloat('Grid','tMin')
        self.tMax     = config.getfloat('Grid','tMax')
        self.Rin      = config.getfloat('Grid','Rin')
        self.Rout     = config.getfloat('Grid','Rout')
        self.Ni       = config.getint('Grid','Ni')
        self.Nj       = config.getint('Grid','Nj')
        self.Nk       = config.getint('Grid','Nk')


