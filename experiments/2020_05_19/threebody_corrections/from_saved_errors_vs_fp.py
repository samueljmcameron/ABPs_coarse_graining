import numpy as np
import matplotlib.pyplot as plt
import sys

import time

sys.path.append('../../analysis_scripts')
import thermodata_analysis as thermo
sys.path.append('../../plotting_scripts')
from jupyterplots import JupyterPlots
import seaborn as sns


figsize = JupyterPlots()

colors = sns.color_palette()

fig,axarr = plt.subplots(2,sharex=True,figsize=[figsize[0],2*figsize[1]])

rho = '0.05'

prefactor = 2**(1./6.)
rho_num = float(rho)*prefactor**2

data = np.loadtxt(f"saved_data/rho_is_{rho}.txt")
fps = data[:,0]*prefactor
Ploads = data[:,1]
Pcalcs_rdf = data[:,2]
Pcalcs_total = data[:,3]+data[:,2]

Pswim = thermo.ideal_swim_press(fps,rho_num)

axarr[0].plot(fps,Ploads*1000,'>-',color=colors[0],
              label=r'$P_L-P_{\mathrm{swim}}-P_{\mathrm{id}}$')
axarr[0].plot(fps,Pcalcs_rdf*1000,'s-',color=colors[2],
              label=r'$P_{ss}^{(2)}$')
axarr[0].plot(fps,Pcalcs_total*1000,'^-',color=colors[3],
              label=r'$P_{ss}^{(2)}+P_{ss}^{(3)}$')
ybottom,ytop=axarr[0].get_ylim()
axarr[0].plot(fps,(Ploads+Pswim)*1000,'o-',color=colors[1],
              label=r'$P_L-P_{\mathrm{id}}$')

axarr[0].set_ylabel(r'$P\times\num{e3}$')
axarr[0].set_ylim(ybottom,ytop)
axarr[0].text(3,6,rf'$\rho={rho_num}$')
#axarr[0].legend(frameon=False)

rho = '0.5'
rho_num = float(rho)*prefactor**2

data = np.loadtxt(f"saved_data/rho_is_{rho}.txt")
fps = data[:,0]*prefactor
Ploads = data[:,1]
Pcalcs_rdf = data[:,2]
Pcalcs_total = data[:,3]+data[:,2]

Pswim = thermo.ideal_swim_press(fps,rho_num)


axarr[1].plot(fps,Ploads+Pswim,'o-',color=colors[1],
              label=r'$P_L-P_{\mathrm{id}}$')
axarr[1].plot(fps,Ploads,'>-',color=colors[0],
              label=r'$P_L-P_{\mathrm{swim}}-P_{\mathrm{id}}$')
axarr[1].plot(fps,Pcalcs_rdf,'s-',color=colors[2],
              label=r'$P_{ss}^{(2)}$')
axarr[1].plot(fps,Pcalcs_total,'^-',color=colors[3],
              label=r'$P_{ss}^{(2)}+P_{ss}^{(3)}$')

axarr[1].set_ylabel(r'$P$')
axarr[1].text(3,1.4,rf'$\rho={rho_num}$')
axarr[1].legend(frameon=False)
axarr[1].set_xlabel('$\mathrm{Pe}$')

fig.subplots_adjust(left=0.2,hspace=0.05,top=0.95)
fig.savefig('results/2020_05_19_threebody_correction'
            +'_P_vs_fp_rho_is_0.05_and_0.5.pdf')
plt.show()

