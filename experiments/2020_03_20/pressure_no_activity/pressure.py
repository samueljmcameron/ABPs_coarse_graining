import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../scripts')

from virialcoeffs import VirialCoeffs

prefix = '../../2020_03_19/raw_data_processing/raw_data/'


tcut = 1000000

fp = 0

fig,axarr = plt.subplots(2,sharex=True)

rhos = ['0.05','0.1','0.2','0.4']

rhonums = np.array(rhos,float)

Davs = np.empty([len(rhos)],float)
Pavs = np.empty([len(rhos)],float)
Perrs = np.empty([len(rhos)],float)

for i,rho in enumerate(rhos):


    vc = VirialCoeffs(l0=fp,epsilon=1,Pi=3./2.)
    
    data = np.loadtxt(prefix + f'pressure_{fp}_{rho}.txt.fixprint')


    ts = data[:,0]
    Ts = data[:,1]
    Ds = data[ts>tcut,2]
    Ps = data[ts>tcut,3]

    Davs[i] = np.mean(Ds)
    Pavs[i] = np.mean(Ps)
    Perrs[i] = np.std(Ps)/np.sqrt(len(Ps))

    dens = float(rho)


from scipy.optimize import curve_fit

def f(x,a,b,c):

    return b*x*x

popt,pcov = curve_fit(f,rhonums[:-1],Pavs[:-1])

rhops = np.linspace(0,0.4,num=100)
pre = np.cbrt(2)
P = (rhops)**2*pre*vc.B_2()

print(popt)
print(vc.B_2()*np.sqrt(2))

axarr[0].plot(rhonums,Davs,'o-')
axarr[1].errorbar(rhonums,Pavs,yerr=Perrs,fmt='o-')
#axarr[1].plot(rhops,f(rhops,*popt),'k-')
axarr[1].plot(rhops,P,'k-')


axarr[0].set_ylabel(r'$D$')

axarr[1].set_ylabel(r'$P$')
axarr[1].set_ylim(0,0.1)
axarr[1].set_xlim(0.02,0.22)
#axarr[1].set_xscale('log')
axarr[1].set_xlabel(r'$\rho$')

axarr[0].legend(frameon=False)


fig.savefig(f'results/pressure.pdf')

plt.show()
