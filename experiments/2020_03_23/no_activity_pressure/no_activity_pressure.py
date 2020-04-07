import numpy as np

import matplotlib.pyplot as plt

import sys

sys.path.append('../../scripts')

from virialcoeffs import VirialCoeffs
from kirkwood import Kirkwood

sys.path.append('../../analysis_scripts')
from rdfloader import RdfLoader


Nsamples = 20
fp = 0


prefix1 = '../../2020_03_19/raw_data_processing/raw_data/'
tcut = 1000000


rhostrings = ['0.05','0.1','0.2','0.4']

# density in WCA r_min units (not lj units)
rhos = np.array(rhostrings,float)*2**(1./3.)

fig1,axarr1 = plt.subplots(len(rhos),ncols=2,sharex=True)

fig1.set_size_inches(8,10)
Davs = np.empty([len(rhos)],float)
Pavs = np.empty([len(rhos)],float)
Perrs = np.empty([len(rhos)],float)

Pcalcs = np.empty([len(rhos)],float)

kw = Kirkwood(l0 = 0)

for i,rho in enumerate(rhostrings):

    data = np.loadtxt(prefix1 + f'pressure_{fp}_{rho}.txt.fixprint')

    ts = data[:,0]
    Ts = data[:,1]
    Ds = data[ts>tcut,2]
    Ps = data[ts>tcut,3]

    Davs[i] = np.mean(Ds)
    Pavs[i] = 2**(1./3.)*np.mean(Ps)
    Perrs[i] = 2**(1./3.)*np.std(Ps)/np.sqrt(len(Ps))

    dens = float(rho)*2**(1./3.)

    file_loc = '../../2020_03_22/correlation/corr_files/'

    rdf = RdfLoader(file_loc,Nsamples)

    rs,g2data = rdf.read_gmeans(fp,rho,savgol=True)

    axarr1[i][1].plot(rs*2**(-1/6),g2data)
    rxs = np.linspace(0.5,2,num=1000,endpoint=True)
    g2spline = kw.g2_spline(rxs,g2data,rs)

    g2integrand = kw.W_2_prime(rxs)*g2spline
    
    axarr1[i][1].plot(rxs,g2spline,'k-')


    
    axarr1[i][0].plot(rxs,g2integrand,'k-')


    Pcalcs[i] = kw.pressure_correct(dens,g2data,rs)
    print(Pcalcs)

axarr1[i][0].set_xlim(0,4)
axarr1[i][0].set_xlabel(r'$r/r_{min}$')




if True:
    fig2,axarr2 = plt.subplots(2,sharex=True)

    axarr2[0].plot(rhos,Davs,'o-')
    axarr2[1].errorbar(rhos,Pavs,yerr=Perrs,fmt='o')

    axarr2[1].plot(rhos,Pcalcs,'k-')

    axarr2[0].set_ylabel(r'$D/D^t$')

    axarr2[1].set_ylabel(r'$P \beta r_{min}^2$')
    #axarr2[1].set_ylim(0,0.1)
    #axarr2[1].set_xlim(0.02,0.22)

    axarr2[1].set_xlabel(r'$\rho r_{min}^2$')

    axarr2[0].legend(frameon=False)

plt.show()



