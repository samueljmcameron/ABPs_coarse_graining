
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../plotting_scripts')
from jupyterplots import JupyterPlots
sys.path.append('../no_interactions_pressure')
from single_pressure import D0,swim

figsize = JupyterPlots()

prefix1 = 'data/'
tcut = -1
fps = np.array([1,5,10,20,40,60,80,100],int)
rhos = ['0.05','0.1','0.2','0.3',
        '0.4','0.5','0.6','0.7']

# density (in lj units)
rho = sys.argv[1]

fig,axarr = plt.subplots(2,sharex=True,figsize=[figsize[0],figsize[0]*2])

fp = int(sys.argv[2])


fname = f'pressure_{fp}_{rho}'


data = np.loadtxt(prefix1 + fname + '.txt.fixprint')


ts = data[:,0]
Ts = data[:,1]
Ds = data[:,2]
Ps = data[:,3]


axarr[0].plot(ts,Ds,'o',label=rf'$f_p={fp}$')
axarr[0].plot(ts,D0(fp)+0*ts,'k-')
axarr[1].plot(ts,Ps,'o',label=rf'$f_p={fp}$')
axarr[1].plot(ts,swim(fp,float(rho))+ts*0,'k-',label=rf'$f_p={fp}$')

axarr[0].set_ylabel(r'$D_{eff}$')
axarr[1].set_ylabel(r'$P_{eff}$')
axarr[1].set_xlabel(r'$t$')
fig.savefig('results_single_pressure/' + fname + '.pdf')
plt.show()


