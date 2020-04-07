import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../scripts')

from virialcoeffs import VirialCoeffs

prefix = '../raw_data_processing/raw_data/'


tcut = 1000000

fps = np.array([0,1,5])#,10,20,40,60,80,100])

fig,axarr = plt.subplots(2,sharex=True)

for rho in ['0.05','0.1','0.2','0.4']:

    Davs = np.empty([len(fps)],float)
    Pavs = np.empty([len(fps)],float)

    Pcalcs = np.zeros([len(fps)],float)

    for i,fp in enumerate(fps):

        vc = VirialCoeffs(l0=fp,epsilon=1,Pi=3./2.)

        data = np.loadtxt(prefix + f'pressure_{fp}_{rho}.txt.fixprint')


        ts = data[:,0]
        Ts = data[:,1]
        Ds = data[ts>tcut,2]
        Ps = data[ts>tcut,3]

        Davs[i] = np.mean(Ds)
        Pavs[i] = np.mean(Ps)

        dens = float(rho)
        Pcalcs[i] = dens**2*vc.B_2()
        
    axarr[0].plot(fps,Davs,'o-',label=rf'$\rho={rho}$')
    axarr[1].plot(fps,Pavs,'o-',label=rf'$\rho={rho}$')
    axarr[1].plot(fps,Pcalcs,'k-',label=rf'$\rho={rho}$')

#axarr[0].set_yscale('log') 
axarr[0].set_ylabel(r'$D$')
axarr[1].set_yscale('log') 
axarr[1].set_ylabel(r'$P$')

#axarr[1].set_xscale('log')
axarr[1].set_xlabel(r'$f_p$')

axarr[0].legend(frameon=False)


fig.savefig(f'results/pressure.pdf')

plt.show()
