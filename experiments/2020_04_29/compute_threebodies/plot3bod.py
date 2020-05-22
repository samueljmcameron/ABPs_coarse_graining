import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../plotting_scripts')

from jupyterplots import JupyterPlots

figsize = JupyterPlots()

nbins = 200
fname = f'g3bod_av_0_0.05_{nbins}bins'

if nbins == 100:
    rij_choices = [1.0625,4.0625,7.8875]
elif nbins == 200:
    rij_choices = [1.08125,4.00625,7.86875]
else:
    print("Need to specify rij_choices in the script for nbins = "
          f"{nbins}")
    exit(1)

data = np.loadtxt('data/' + fname + '.g3bod',skiprows=4)

thetas_full = data[:,1]
rijs_full = data[:,2]
riks_full = data[:,3]

g3s_full = data[:,4]

theta_choose = 0.958186
mask1 = np.isclose(thetas_full,theta_choose)
rijs_full = rijs_full[mask1]

riks_full = riks_full[mask1]

g3s_full = g3s_full[mask1]



fig,axarr = plt.subplots(3,sharex=True,figsize=[figsize[0],figsize[1]*3])

for i,rij_choose in enumerate(rij_choices):

    mask1 = np.isclose(rij_choose,rijs_full)

    riks = riks_full[mask1]

    g3s = g3s_full[mask1]
    
    label=rf'$\theta={theta_choose},\:r_{{12}}={rij_choose}$'
    
    axarr[i].plot(riks,g3s,label=label)

    axarr[i].legend(frameon=False)
    axarr[i].set_ylabel(r'$g^{(3)}(r_{12},r_{13},\theta)$')
axarr[2].set_xlabel(r'$r_{13}$')

fig.subplots_adjust(left = 0.2)
fig.savefig('results/' + fname + '.pdf')
plt.show()
