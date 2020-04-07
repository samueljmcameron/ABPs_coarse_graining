"""
Plot steady state pressure vs activity and vs density

"""



import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from single_pressure import D0, swim

sys.path.append('../../plotting_scripts')
from jupyterplots import JupyterPlots



figsize = JupyterPlots()

colors = sns.color_palette()


prefix1 = 'data/'

fps = np.array([0,1,5,10,20,40,60,80,100],int)


# density (in lj units)
rho = '0.4'
tcut = 200000

fig,axarr = plt.subplots(2,2,sharex='col',
                         figsize=[1.5*figsize[0],2*figsize[1]])

Pavs = np.empty([len(fps)],float)
Perrs = np.empty([len(fps)],float)
Dlasts = np.empty([len(fps)],float)


for i,fp in enumerate(fps):

    data = np.loadtxt(prefix1 + f'pressure_{fp}_{rho}.txt.fixprint')


    ts = data[:,0]
    Ts = data[:,1]
    Ds = data[ts>tcut,2]
    Ps = data[ts>tcut,3]

    #ll = LogLoader(prefix1 + f'log_{fp}_{rho}')

    #ts = ll.data['Step']
    #RMSs = ll.data['c_mymsdd[4]']
    #Ps = ll.data['c_press']

    
    
    Pavs[i] = np.mean(Ps)
    Perrs[i] = np.std(Ps)/np.sqrt(len(Ps))
    Dlasts[i] = Ds[-1]




axarr[0][0].errorbar(fps,Pavs,yerr=Perrs,fmt='.',color=colors[0],
                     label=rf'$\rho={rho}$')
axarr[0][0].plot(fps,swim(fps,float(rho))+fps*0,'k:',
                 label=rf'$\rho f_P^2/6$')
axarr[0][0].set_ylabel(r'$<P>$')

axarr[1][0].plot(fps,Dlasts,'.',color=colors[0],
                 label=rf'$\rho={rho}$')
axarr[1][0].plot(fps,D0(fps)+fps*0,'k--',
                 label=rf'$1+f_P^2/6$')
axarr[1][0].set_ylabel(r'$D(t_{\mathrm{final}})$')

axarr[1][0].set_xlabel(r'$f_P$')


fp = 100

rhos = ['0.10','0.15','0.20','0.25',
        '0.30','0.35','0.4','0.5','0.6','0.7']

Pavs = np.empty([len(rhos)],float)
Perrs = np.empty([len(rhos)],float)
Dlasts = np.empty([len(rhos)],float)

for i,rho in enumerate(rhos):

    if rho == '0.7':
        prefix1 = 'data2/'
    else:
        prefix1 = 'data/'
        
    data = np.loadtxt(prefix1 + f'pressure_{fp}_{rho}.txt.fixprint')


    ts = data[:,0]
    Ts = data[:,1]
    Ds = data[ts>tcut,2]
    Ps = data[ts>tcut,3]

    Pavs[i] = np.mean(Ps)
    Perrs[i] = np.std(Ps)/np.sqrt(len(Ps))
    Dlasts[i] = Ds[-1]


rhonum = np.array(rhos,float)

axarr[0][1].errorbar(rhonum,Pavs,yerr=Perrs,fmt='.',color=colors[1],
             label=rf'$f_P={fp}$')
axarr[0][1].plot(rhonum,swim(fp,rhonum)+rhonum*0,'k:',
         label=rf'$\rho f_P^2/6$')

axarr[1][1].plot(rhonum,Dlasts,'.',color=colors[1],
                 label=rf'$f_P={fp}$')
axarr[1][1].plot(rhonum,D0(fp)+rhonum*0,'k--',
         label=rf'$1+f_P^2/6$')
#axarr[1][1].set_ylabel(r'$D(t_{\mathrm{final}})$')

axarr[1][1].set_xlabel(r'$\rho$')

axarr[0][0].legend(frameon=False)
axarr[0][1].legend(frameon=False)
axarr[1][0].legend(frameon=False)
axarr[1][1].legend(frameon=False)
fig.subplots_adjust(hspace=0.05,wspace=0.05)
fig.savefig('results/pressure_no_forces.pdf')

plt.show()


