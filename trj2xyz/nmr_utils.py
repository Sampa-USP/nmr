class post_processing():
    """
    classes :
    read_txt -> used internally to get txt file;
    make_fft -> uses fftw to make the fourier transform in the interval [0,2pi];
    get_t2 -> get the t2 time using equations from paper doi:10.1016/j.jmr.2017.02.001;
    """
    def __init__(self,path,dt,thermo_step,tsim):
        """
        path = 'path to the output generate from corr.x';
        dt = 'dt choose in the corr-water.in file for corr.x input';
        thermo_step = 'thermo step in the lammps file';
        tsim = 'number of frames in the .lammpstrj file';
        """
        self.path = path
        self.dt = dt
        self.thermo_step = thermo_step
        self.tsim = tsim
        self.read_txt()
        self.make_fft()

    def read_txt(self):
        file = open(self.path,"r")
        matriz = file.read().split('\n')[2:-1]
        file.close()

        #self.matriz = matriz
        linhas = len(matriz)
        self.data = np.array([matriz[i].split() for i in range(linhas)]).astype(float)
        self.lines = len(self.data[:,0])

        self.time = self.data[:,0]*(self.thermo_step*self.dt*self.tsim)
        self.norm_gr = self.data[:,1]
        self.norm_gt = self.data[:,2]
        self.gr = self.data[:,3]
        self.gt = self.data[:,4]

    def make_fft(self):
        self.jr = np.fft.fft(self.gr,norm="ortho").real
        self.jt = np.fft.fft(self.gt,norm="ortho").real

    def get_t2(self,dt,thermo):
        pass

    def help(self):
        text = """I'm need to explain some points about the implementation:
\n using only double precision float values;
\n the fourier transform is in the interval [0,2pi)"""
        #normalizer = dt*thermo
