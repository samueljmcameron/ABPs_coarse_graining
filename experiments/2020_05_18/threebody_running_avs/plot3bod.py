import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../plotting_scripts')

from jupyterplots import JupyterPlots



figsize = JupyterPlots()

fp,rho = sys.argv[1],sys.argv[2]

nbins = 200
theta_max = np.pi

r_min,r_max = 0.5,8.0

dtheta = theta_max/nbins

dr = (r_max-r_min)/nbins


fname = f'g3bod_av_{fp}_{rho}'

if nbins == 100:
    rij_choices = [1.0625,4.0625,7.8875]
elif nbins == 200:
    rij_choices = [1.08125,4.00625,7.86875]
else:
    print("Need to specify rij_choices in the script for nbins = "
          f"{nbins}")
    exit(1)

theta_choose = dtheta*66+0.5*dtheta
print(theta_choose)
    
data = np.loadtxt('data/' + fname + '.g3bod',skiprows=4)

thetas_full = data[:,1]
rijs_full = data[:,2]
riks_full = data[:,3]

g3s_full = data[:,4]


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
    axarr[i].plot(riks,g3s*0+1,'k--')

    axarr[i].legend(frameon=False)
    axarr[i].set_ylabel(r'$g^{(3)}(r_{12},r_{13},\theta)$')
axarr[2].set_xlabel(r'$r_{13}$')

fig.subplots_adjust(left = 0.2)
fig.savefig('results/' + fname + '.pdf')
plt.show()
