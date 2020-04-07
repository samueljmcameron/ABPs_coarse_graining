"""
MEASUREMENT OF DIFFUSION COEFFICIENT IN THIS CODE IS INCORRECT!!!
"""


import numpy as np

import matplotlib.pyplot as plt

import sys

sys.path.append('../../scripts')

from virialcoeffs import VirialCoeffs
from kirkwood import Kirkwood

sys.path.append('../../analysis_scripts')
from rdfloader import RdfLoader
from logloader import LogLoader

Nsamples = 20

def diffusion(tsim,fp,dt = 1e-5):

    # this equation is wrong!!!!!!!!!!!!!!!!!
    
    Dt = 1.
    taur = 1./(3*Dt)

    t = tsim*dt
    
    return Dt + fp**2*taur**2/12.*(2/taur)*(1-np.exp(-2*t/taur/10.))

def D0(fp):

    Dt = 1
    taur = 1./(3*Dt)

    return fp**2*taur/2.

prefix1 = 'data/'
tcut = 1000000

fps = np.array([1,5,10,20,40,60,80,100],int)


# density in WCA r_min units (not lj units)
rho = '0.001'

#Davs = np.empty([len(rhos)],float)
#Pavs = np.empty([len(rhos)],float)
#Perrs = np.empty([len(rhos)],float)

fig,axarr = plt.subplots(2,sharex=True)

for i,fp in enumerate(fps):

    data = np.loadtxt(prefix1 + f'pressure_{fp}_{rho}.txt.fixprint')

    ts = data[:,0]
    Ts = data[:,1]
    Ds = data[:,2]
    Ps = data[:,3]


    #Davs[i] = np.mean(Ds)
    #Pavs[i] = 2**(1./3.)*np.mean(Ps)
    #Perrs[i] = 2**(1./3.)*np.std(Ps)/np.sqrt(len(Ps))
    axarr[0].plot(ts,Ds,'o',label=rf'$f_p={fp}$')
    axarr[0].plot(ts,D0(fp)+0*ts,'k-')
    axarr[1].plot(ts,Ps,'o',label=rf'$f_p={fp}$')

plt.show()


