import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../scripts')

from virialcoeffs import VirialCoeffs

prefix = '../raw_data_processing/raw_data/'


tcut = 1000000

fps = np.array([0,1,5,10])

fig,axarr = plt.subplots(4,sharex=True)

rs = np.linspace(0.5,1.2,num=301,endpoint=True)

for i,rho in enumerate(['0.05','0.1','0.2','0.4']):

    Davs = np.empty([len(fps)],float)
    Pavs = np.empty([len(fps)],float)

    Pcalcs = np.empty([len(fps)],float)

    for j,fp in enumerate(fps):

        vc = VirialCoeffs(l0=fp,epsilon=1,Pi=3./2.)

        axarr[i].plot(rs,vc.f_12(rs),'o-',label=rf'$f_p={fp}$')
        

axarr[0].set_ylabel(r'$f_{12}$')
axarr[1].set_ylabel(r'$f_{12}$')
axarr[2].set_ylabel(r'$f_{12}$')
axarr[3].set_ylabel(r'$f_{12}$')

axarr[3].set_xlabel(r'$r$')

axarr[0].legend(frameon=False)


#fig.savefig(f'results/pressure.pdf')

plt.show()
