import numpy as np
import matplotlib.pyplot as plt

prefix = '../raw_data_processing/raw_data/'

def diffusion(t,fp):

    return (1+2./3.*fp**2*(1-np.exp(-6*t*0.00001)))


fps = np.array([0,1,5,10,20,40,60,80,100])

for rho in ['0.05','0.1','0.2','0.4']:

    fig,axarr = plt.subplots(2,sharex=True)

    for i,fp in enumerate(fps):


        data = np.loadtxt(prefix + f'pressure_{fp}_{rho}.txt.fixprint')


        ts = data[:,0]
        Ts = data[:,1]
        Ds = data[:,2]
        Ps = data[:,3]

        axarr[0].plot(ts,Ds,label=rf'$f_p={fp}$')
        axarr[0].plot(ts,diffusion(ts,fp),'k--')
        axarr[1].plot(ts,Ps,label=rf'$f_p={fp}$')


    axarr[0].set_ylabel(r'$D$')
    axarr[1].set_ylabel(r'$P$')
    axarr[1].set_ylim(0.04,20)
    axarr[1].set_yscale('log')
    axarr[1].legend(frameon=False,ncol=3)


    fig.savefig(f'results/pressure_{rho}.pdf')

plt.show()
