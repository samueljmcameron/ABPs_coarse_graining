import numpy as np
import scipy.signal as signal
import scipy.interpolate as interpolate

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


    
    def g2_spliner_setup(self,g2data,rs,kind='quadratic'):

        g2spline = interpolate.interp1d(rs,g2data,
                                        kind='quadratic')
        
        return g2spline

    
    def g2_spline(self,r,g2data,rs,g2tol=1e-12,r0_ret=False):

        gs = self.g2_spliner_setup(g2data,rs)(r)


        try:
            r0 = r[gs<g2tol][-1]
        except:
            print('error! spline region isn\'t going to small'
                  'enough r to find g(r)=0 !!')

        if r0_ret:
            return gs[r>r0],r[r>r0]
        else:
            return np.where(r<r0,0,gs)
    

if __name__ == "__main__":

    pass
