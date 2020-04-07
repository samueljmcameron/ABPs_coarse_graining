import numpy as np

import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.special as spec
import matplotlib.pyplot as plt

from softsphere import SoftSphere

from effectivepotentials import W2potential, Chi2potential

class Kirkwood(SoftSphere):

    def __init__(self,potential='W_2',Omega=None,l0=1,epsilon=1,
                 Pi=1,convert_units=True):

        super().__init__(Omega=Omega,l0=l0,epsilon=epsilon,Pi=Pi)

        if convert_units:
            self.prefactor = 2**(1./6.)
        else:
            self.prefactor = 1.

        if potential == 'W_2':
            self.potential = W2potential(Omega=Omega,l0=l0,
                                         epsilon=epsilon,Pi=Pi)
        elif potential == 'chi_2':
            self.potential = Chi2potential(Omega=Omega,l0=l0,
                                           epsilon=epsilon,Pi=Pi)
        
        return

    
    def g2_spliner_lam(self,g2data,rs):

        g2spline = interpolate.interp1d(rs/self.prefactor,g2data,
                                        kind='quadratic')
        
        return g2spline

    def g2_spline_scalar(self,r,g2data,rs):

        g = self.g2_spliner_lam(g2data,rs)(r)

        if g < 1e-3:

            g = 0

        return g
    
    def g2_spline(self,r,g2data,rs,g2tol=1e-12,r0_ret=False):

        gs = self.g2_spliner_lam(g2data,rs)(r)


        try:
            r0 = r[gs<g2tol][-1]
        except:
            print('error! spline region isn\'t going to small'
                  'enough r to find g(r)=0 !!')

        if r0_ret:
            return gs[r>r0],r[r>r0]
        else:
            return np.where(r<r0,0,gs)
    
    
    def g2_integrand(self,r,g2spline):
        
        return g2spline*self.potential.prime(r)*r*r
    
    def g2_integral(self,g2data,rs,min_r = 0.5,
                    max_r=10,g2tol=1e-12):

        # define 'exploratory' r region to determine what a suitable
        # lower integration bound will be.
        rexplore = np.linspace(min_r,max_r,num=10000,endpoint=True)

        # find the smallest r0 where g2_spline is >= g2tol
        # this step is necessary to avoid have a really small number
        # (i.e. the correlation function g(r) at r << 1) which
        # should be zero, being multiplied by the potential at r<<1
        # which is quite large, and so getting a nonsense result.
        spline,rints = self.g2_spline(rexplore,g2data,rs,g2tol=g2tol,
                                      r0_ret=True)

        f = self.g2_integrand(rints,spline)

        return integrate.simps(f,rints)
        

    def pressure_correct(self,rho,g2data,rs):


        # here rho must be in units of rmin (i.e. distance of WCA
        # potential's minimum and cutoff), NOT lj sigma units

        a = -rho**2*self.g2_integral(g2data,rs)*np.pi/2

        return a


if __name__ == "__main__":

    import sys

    sys.path.append('../analysis_scripts')
    from rdfloader import RdfLoader

    def P_of_rho(rho,g2data,rs):

        Kw = Kirkwood(l0 = 0)

        coeff = Kw.pressure_correct(g2data,rs)*rho**2

        return coeff
    
    Nsamples = 20
    fp = 0


    prefix1 = '../2020_03_19/raw_data_processing/raw_data/'
    tcut = 1000000


    rhostrings = ['0.05','0.1','0.2','0.4']

    rhos = np.array(rhostrings,float)

    fig1,axarr1 = plt.subplots(len(rhos),ncols=2,sharex=True)

    fig1.set_size_inches(8,10)
    Davs = np.empty([len(rhos)],float)
    Pavs = np.empty([len(rhos)],float)
    Perrs = np.empty([len(rhos)],float)

    Pcalcs = np.empty([len(rhos)],float)

    kw = Kirkwood(l0 = 0)

    for i,rho in enumerate(rhostrings):

        data = np.loadtxt(prefix1 + f'pressure_{fp}_{rho}.txt')

        ts = data[:,0]
        Ts = data[:,1]
        Ds = data[ts>tcut,2]
        Ps = data[ts>tcut,3]

        Davs[i] = np.mean(Ds)
        Pavs[i] = np.mean(Ps)
        Perrs[i] = np.std(Ps)/np.sqrt(len(Ps))

        dens = float(rho)

        file_loc = '../2020_03_22/correlation/corr_files/'
        
        rdf = RdfLoader(file_loc,Nsamples)

        rs,g2data = rdf.read_gmeans(fp,rho,savgol=True)

        axarr1[i][1].plot(rs*2**(-1/6),g2data)
        rxs = np.linspace(0.5,2,num=1000,endpoint=True)
        func = kw.g2_spline(rxs,g2data,rs)
        
        axarr1[i][1].plot(rxs,func,'k-')
        



        #axarr1[i][1].plot(rs[rs<rs[-1]]*2**(-1/6),gs[rs<rs[-1]],'.')
        #axarr1[i][1].plot(rxs,func,'k-')
        
        Pcalcs[i] = P_of_rho(dens,g2data,rs)
        print(Pcalcs)
        
    axarr1[i][0].set_xlim(0,4)
    axarr1[i][0].set_xlabel(r'$r/r_{min}$')

        


    if True:
        fig2,axarr2 = plt.subplots(len(rhos),sharex=True)

        axarr2[0].plot(rhos,Davs,'o-')
        axarr2[1].errorbar(rhos,Pavs,yerr=Perrs,fmt='o')

        axarr2[1].plot(rhos,Pcalcs,'k-')

        axarr2[0].set_ylabel(r'$D$')

        axarr2[1].set_ylabel(r'$P$')
        #axarr2[1].set_ylim(0,0.1)
        #axarr2[1].set_xlim(0.02,0.22)

        axarr2[1].set_xlabel(r'$\rho$')

        axarr2[0].legend(frameon=False)

    plt.show()




        
    """
    prefix1 = '../2020_03_19/raw_data_processing/raw_data/'

    prefix2 = '../2020_03_19/raw_data_processing/pickled_data/'
        
    tcut = 1000000

    fp = 0

    rhostrings = ['0.05','0.1','0.2','0.4']

    rhos = np.array(rhostrings,float)

    fig1,axarr1 = plt.subplots(len(rhos),ncols=2,sharex=True)

    fig1.set_size_inches(8,10)
    Davs = np.empty([len(rhos)],float)
    Pavs = np.empty([len(rhos)],float)
    Perrs = np.empty([len(rhos)],float)

    Pcalcs = np.empty([len(rhos)],float)

    kw = Kirkwood(l0 = 0)

    for i,rho in enumerate(rhostrings):

        data = np.loadtxt(prefix1 + f'pressure_{fp}_{rho}.txt')


        ts = data[:,0]
        Ts = data[:,1]
        Ds = data[ts>tcut,2]
        Ps = data[ts>tcut,3]

        Davs[i] = np.mean(Ds)
        Pavs[i] = np.mean(Ps)
        Perrs[i] = np.std(Ps)/np.sqrt(len(Ps))

        dens = float(rho)

        dc_name = prefix2 + f'ret_o_{fp}_{rho}'

        ret_o = load_obj(dc_name)

        gmatrix = ret_o['sum_g']
        Nsamples = ret_o['g_cnt']
        
        rs = np.array(gmatrix[:,0]/Nsamples).flatten()
        gs = np.array(gmatrix[:,1]/Nsamples).flatten()

        print(rs[2]-rs[1])
        rxs = np.linspace(0.5,2,num=1000,endpoint=True)
        func = kw.g2_spline(rxs,gs,rs)

        axarr1[i][1].plot(rs[rs<rs[-1]]*2**(-1/6),gs[rs<rs[-1]],'.')
        axarr1[i][1].plot(rxs,func,'k-')
        
        Pcalcs[i] = P_of_rho(dens,gs,rs)
        
    axarr1[i][0].set_xlim(0,4)
    axarr1[i][0].set_xlabel(r'$r/r_{min}$')

        


    if True:
        fig2,axarr2 = plt.subplots(len(rhos),sharex=True)

        axarr2[0].plot(rhos,Davs,'o-')
        axarr2[1].errorbar(rhos,Pavs,yerr=Perrs,fmt='o-')

        axarr2[1].plot(rhos,Pcalcs,'k-')

        axarr2[0].set_ylabel(r'$D$')

        axarr2[1].set_ylabel(r'$P$')
        #axarr2[1].set_ylim(0,0.1)
        #axarr2[1].set_xlim(0.02,0.22)

        axarr2[1].set_xlabel(r'$\rho$')

        axarr2[0].legend(frameon=False)

    plt.show()
    """
