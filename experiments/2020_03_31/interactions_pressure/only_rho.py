"""
Plot steady state pressure vs activity and vs density

"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

sys.path.append('../../plotting_scripts')
from jupyterplots import JupyterPlots

sys.path.append('../no_interactions_pressure')
from single_pressure import D0, swim


def pid(rho,fp):

    return rho*(1+fp**2/(2*3))

figsize = JupyterPlots()

colors = sns.color_palette()


prefix1 = 'data/'

fps = np.array([0,1,5,10,20,40,80],int)


# density (in lj units)
rho = '0.4'
tcut = 200000

fig,ax = plt.subplots()



fps = np.array([0,1,5,10,20,40,80],int)

for j,fp in enumerate(fps):

    rhos = ['0.05','0.1','0.2','0.3',
            '0.4','0.5','0.6','0.7']

    Pavs = np.empty([len(rhos)],float)
    Perrs = np.empty([len(rhos)],float)
    Dlasts = np.empty([len(rhos)],float)
    prefix1 = 'data/'

    for i,rho in enumerate(rhos):

        data = np.loadtxt(prefix1 + f'pressure_{fp}_{rho}.txt.fixprint')


        ts = data[:,0]
        Ts = data[:,1]
        Ds = data[ts>tcut,2]
        Ps = data[ts>tcut,3]

        Pavs[i] = np.mean(Ps)
        Perrs[i] = np.std(Ps)/np.sqrt(len(Ps))
        Dlasts[i] = Ds[-1]


    rhonum = np.array(rhos,float)

    ax.plot(rhonum,(Pavs+rhonum)*rhonum/pid(rhonum,fp),'.-',
                label=rf'$Pe={fp/3.:.2f}$')

ax.plot(rhonum,rhonum,'k:')


ax.set_xlabel(r'$\rho$')
ax.set_ylabel(r'$P\rho/fp^2$')
ax.set_ylim(0,0.6)
ax.legend(frameon=False)

fig.subplots_adjust(hspace=0.05,wspace=0.05)
fig.savefig('results/pressure_rho_only.pdf')

plt.show()


