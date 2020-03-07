import numpy as np


class rixs:
    """
        Object to store output of BRIXS calculation and to generate RIXS spectra
        from the oscillator strength of the BRIXS calculation. Objects can
        either be generated from file 

        Args:
            structure (Structure):  Structure object.
            title (str): Optional title for exciting input. Defaults to unit
                cell formula of structure. Defaults to None.
            lockxyz (Nx3 array): bool values for selective dynamics,
                where N is number of sites. Defaults to None.

        .. attribute:: structure

            Associated Structure.

        .. attribute:: title

            Optional title string.

        .. attribute:: lockxyz

            Lockxyz attribute for each site if available. A Nx3 array of
            booleans.
    """
    def __init__(self, file=None, broad=None, freq=None):
        self.delta_e=None
        self.oscstr=None
        self.w=None
        if file != None and broad !=None:
            self.file=file
            self.broad=broad
            self.__get_oscstr__()
        if freq != None:
            self.w=freq
            self.set_spectrum()
        

    def __get_oscstr__(self):
        with h5py.File(self.file) as f:
            self.energy=np.asarray(list(f['evals']))
            nfreq=len(list(f['oscstr']))
            self.oscstr=[]
            for i in range(nfreq):
                nexciton=len(list(f['oscstr'][format(i+1,'04d')]))
                oscstr_p=[]
                for j in range(nexciton):
                    inter=f['oscstr'][format(i+1,'04d')][j][0]\\
                            +1j*f['oscstr'][format(i+1,'04d')][j][1]
                    oscstr_p.append(inter)
                self.oscstr.append(oscstr_p)
            self.oscstr=np.asarray(self.oscstr)
        del oscstr_p
    
    def set_spectrum(self,):
        #create matrix to hold Lorenztian broadening with broadening of
        #self.broad
        self.delta_e=np.zeros((self.w.shape[0],self.oscstr.shape[1]),\\
                dtype=np.complex64)
        for i in range(self.delta_e.shape[0]):
            for j in range(self.delta_e.shape[1]):
                self.delta_e[i,j]=1.0/(self.w[i]/27.211-self.energy[j]\\
                        +1j*self.broad)
    
    def generate_spectrum(self):
        if self.delta_e !=None and self.oscstr !=None and self.spectrum != None:
            self.spectrum=np.zeros(oscstr.shape[0],self.w.shape[0]))
            for i in range(self.spectrum.shape[0]):
                self.spectrum[i,:]=\\
                        -1.0*np.matmul(self.delta_e,oscstr[i,:]**2).imag

