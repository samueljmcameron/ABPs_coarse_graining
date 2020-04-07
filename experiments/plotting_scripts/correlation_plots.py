import numpy as np
import matplotlib.pyplot as plt

import sys


from jupyterplots import JupyterPlots


def plot_correlations_and_pressure(fp,rhostrings,Nsamples,prefix1,prefix2,
                                   tcut,Kirkwood,RdfLoader,potential='W_2'):
    """
    Utility plotting function of local and global simulation quantities.
    
    Figure 1 shows local radial distribution function g_2(r) as a
    function of pair separation r, as well as g_2(r)*W_2'(r), where
    W_2'(r) is the effective two-body interaction of ABPs (equal to
    the WCA potential when activity is turned off, fp = 0). This
    figure has two columns and four rows, for eight subplots total.
    Column 1 is the g_2(r) (raw data and Savitzky-Golay filtered),
    column 2 is the g_2(r)W_2'(r). Rows correspond to the four
    different number densities that these quantities were measured
    at.
    
    Figure 2 shows global thermodynamic properties, namely the
    diffusion coefficient, D, defined via D = <r^2>'(t)/4 (row 1); 
    the non-ideal pressure correction, P, measured directly from
    simulation and the pressure correction calculated via
    g_2(r)W'(r)r^2, P_{calc} (row2); the fractional error between P
    and P_{calc} (row 3).
    
    """
    # make nice plots by changing matplotlib defaults
    figsize_default = JupyterPlots()    
    
    # convert to WCA (r_min) units
    rhos = np.array(rhostrings,float)*2**(1./3.)

    # first figure will plot the product of W_2'(r) and
    # g(r), where g(r) is the rdf
    # it will also plot the bare g(r) and its smoothed
    # version (smoothed using savgol_filter from scipy)


    
    fig1size = [figsize_default[0]*1.5,figsize_default[1]*1.5]
    fig1,axarr1 = plt.subplots(len(rhos)//2,ncols=4,sharex=True,
                              figsize = fig1size)

    
    # define arrays to store averaged thermodynamic data,
    # which will be the diffusion coefficient and the pressure
    
    # note that diffusion coefficient is determined via
    # the relationship derivative(<R^2>)/(4*dt), and so is
    # not the same as the translational or rotational
    # diffusivity.
    Davs = np.empty([len(rhos)],float)
    Pavs = np.empty([len(rhos)],float)
    # also include error bars
    Derrs = np.copy(Davs)
    Perrs = np.copy(Pavs)

    # load in class to compute kirkwood stress via rdf.
    # l0 = fp_d*rmin*beta, where fp_d is the activity in
    # real units.
    # However, our analytic calculations (e.g. potentials)
    # are defined in a length scale rmin, where as
    # the simulations are defined in terms of sigma.
    # so fp = fp_d*sigma*beta in the simulations. Therefore,
    # l0/fp = rmin/sigma = 2**(1./6.), so

    l0 = fp*2**(1./6.)

    # Pi = D^r*rmin**2/(2D^t) in analytic calculations, but
    # Pi_s = D^r*sigma**2/(2D^t) in the simulations. So
    # Pi/Pi_s = (rmin/sigma)**2 = 2**(1./3.). In the
    # simulations, we hold D^r = 3D^t/sigma**2. Therefore,
    # Pi_s = 3./2. and

    Pi = 3./2.*2**(1./3.)

    kw = Kirkwood(l0 = l0,Pi=Pi,potential = potential)

    # finally, write an array for storing the pressure values
    # calculated via the rdf.
    
    Pcalcs = np.empty([len(rhos)],float)

    for i, rho in enumerate(rhostrings):
        
        data = np.loadtxt(prefix1 + f'pressure_{fp}_{rho}.txt')
        
        # load in observables as functions of simulation time
        ts = data[:,0] # timesteps
        Ts = data[:,1] # ABP mean kinetic energy
        Ds = data[ts>tcut,2] # diffusion coefficient 
        Ps = data[ts>tcut,3] # pressure
    
        # convert to WCA (r_min) units
        Davs[i] = np.mean(Ds)
        Derrs[i] = np.std(Ds)/np.sqrt(len(Ds))
        Pavs[i] = 2**(1./3.)*np.mean(Ps)
        Perrs[i] = 2**(1./3.)*np.std(Ps)/np.sqrt(len(Ps))
    
        dens = float(rho)*2**(1./3.)
            
        rdf = RdfLoader(prefix2,Nsamples)

        rs,g2data = rdf.read_gmeans(fp,rho,savgol=True)
    
        axarr1[0][i].plot(rs*2**(-1./6.),g2data,'o-',ms=1)
    
        rxs = np.linspace(0.5,5,num=1000)
    
        g2spline = kw.g2_spline(rxs,g2data,rs)
    
        g2timesW2 = kw.potential.prime(rxs)*g2spline
        
        axarr1[0][i].plot(rxs,g2spline,'k-')
    
        axarr1[1][i].plot(rxs,g2timesW2,'k-')
    
        if i == 0:
            axarr1[1][i].set_xlabel(r'$\tilde{r}\equiv r/r_{min}$')
        else:
            axarr1[1][i].set_xlabel(r'$\tilde{r}$')
    
        Pcalcs[i] = kw.pressure_correct(dens,g2data,rs)
    
    axarr1[1][i].set_xlim(0,4)
    axarr1[0][0].set_ylabel(r'$g_2(\tilde{r})$')
    axarr1[1][0].set_ylabel(r'$g_2(\tilde{r})\tilde{W}_2^{\prime}(\tilde{r})$')

    fig1.subplots_adjust(hspace=0.25,wspace=0.5)

    fig2size = [figsize_default[0],figsize_default[1]*2]

    # now plot the global thermodynamic quantities
    fig2,axarr2 = plt.subplots(3,sharex=True,figsize=fig2size)


    axarr2[0].errorbar(rhos,Davs,yerr=Derrs,fmt='o-')
    axarr2[1].errorbar(rhos,Pavs,yerr=Perrs,fmt='o',ms=2)
    axarr2[1].plot(rhos,Pcalcs,'k-')
    
    axarr2[2].plot(rhos,np.abs(Pcalcs - Pavs)/Pavs,'o-')
    axarr2[2].set_yscale('log')

    axarr2[0].set_ylabel(r'$\tilde{D}\equiv D/D^t$')
    axarr2[1].set_ylabel(r'$\tilde{P}\equiv P\beta (r_{min})^2$')
    
    axarr2[2].set_ylabel(r'$\frac{|\tilde{P}-\tilde{P}_{calc}|}{\tilde{P}}$')

    axarr2[2].set_xlabel(r'$\tilde{\rho}\equiv\rho (r_{min})^2$')

    axarr2[0].legend(frameon=False)
    plt.show()

    return

if __name__ == "__main__":


    sys.path.append('../analysis_scripts')
    from rdfloader import RdfLoader

    sys.path.append('../scripts')
    from kirkwood import Kirkwood


    Nsamples = 20
    fp = 0

    # when computing average pressure, only consider 
    # pressure data from times later than tcut
    tcut = 1000000

    # location of time dependent pressure data measured via lammps
    prefix1 = '../2020_03_19/raw_data_processing/raw_data/'

    # location of rdf data computed via lammps rerun command
    prefix2 = '../2020_03_22/correlation/correlations/'

    # list of densities (in LJ units)
    rhostrings = ['0.05','0.1','0.2','0.4']

    plot_correlations_and_pressure(fp,rhostrings,Nsamples,prefix1,prefix2,
                                   tcut,Kirkwood,RdfLoader,potential='W_2')
    
