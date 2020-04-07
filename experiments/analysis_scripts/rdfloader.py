import numpy as np
import scipy.signal as signal

class RdfLoader():

    def __init__(self,file_loc,Nsamples, nbins=2000):

        self.file_loc = file_loc
        self.Nsamples = Nsamples
        self.nbins = nbins

        return
    
    def read_gs(self,fp,rho):

        gs = np.empty([self.nbins,self.Nsamples],float)
        
        for i in range(self.Nsamples):
            fname = self.file_loc
            fname += f'g_{fp}_{rho}_{i}.rdf'
            data = np.loadtxt(fname)

            rbins = data[:,1]
            gs[:,i] = data[:,2]
            
        return rbins,gs
        

    def read_gmeans(self,fp,rho,savgol=False,
                    window = 9,order = 2):

        rbins,gs = self.read_gs(fp,rho)

        gcorrs = np.mean(gs,axis=1)
        gstd = np.std(gs,axis=1)/np.sqrt(self.Nsamples)

        if savgol:

            return rbins,signal.savgol_filter(gcorrs,window,order)

        return rbins,gcorrs,gstd



if __name__ == "__main__":

    pass
