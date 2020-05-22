import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append("../../scripts")

from effectivepotentials import W2potential

sys.path.append('../../plotting_scripts')

from jupyterplots import JupyterPlots

if __name__=="__main__":

    figsize = JupyterPlots()
    
    colors = sns.color_palette()

    # start by plotting at different Omegas
    
    Pi = 3./2.

    l0s = np.array([0.1,10.0],float)
    fig,axarr = plt.subplots(2,sharex=True,figsize=[figsize[0],figsize[1]*2])

    axin = axarr[0].inset_axes([2,1.0,2,2],
                               transform=axarr[0].transData)

    r0 = 0.9
    for i,l0 in enumerate(l0s):
    
        ss = W2potential(l0=l0,Pi=Pi)


        r0new = ss.r_0()
        if r0new < r0:
            r0 = r0new
        rs = np.linspace(r0/1.1,10,num=10000,endpoint=True)

        axarr[0].plot(rs,ss.evaluate(rs),color=colors[i],
                      label=rf'$\mathrm{{Pe}}={l0}$')
        axin.plot(rs,ss.evaluate(rs),color=colors[i],
                  label=rf'$\mathrm{{Pe}}={l0}$')
        axin.plot(rs,rs*0,'k--',lw=0.5)

        axarr[1].plot(rs,ss.prime(rs),color=colors[i],
                      label=rf'$\mathrm{{Pe}}={l0}$')

        if i == len(l0s)-1:
            r0label=r'$r_0$'
        else:
            r0label=None
            
        axarr[0].plot(r0new,ss.evaluate(r0new),'ko',ms=2,label=r0label)
        axin.plot(r0new,ss.evaluate(r0new),'ko',ms=2)

    #ax.set_yscale('log')

    #ax.set_ylim(bottom=1e-3,top=max(ss.w(rs))*1000)
    #ax.set_ylim(bottom=-max(ss.w(rs)),top=max(ss.w(rs))*2)

    axarr[0].set_ylabel(r'$W_2(r)$')
    axarr[0].legend(frameon=False)

    axin.set_xlim(left=0.9,right=1.5)
    axin.set_ylim(bottom=-0.0005,top=0.001)
    axin.set_xticks([])
    axin.set_yticks([])
    axarr[0].legend(frameon=False)

    axarr[1].set_xlim(left=r0/1.4)
    axarr[1].set_ylim(bottom=-20,top=20)
    axarr[1].set_ylabel(r'$W_2^{\prime}(r)$')
    axarr[1].set_xlabel(r'$r$')


    fig.subplots_adjust(left=0.15,right=0.85,bottom=0.2)
    fig.savefig('results/2020_05_15_W_2_and_prime.pdf')
    plt.show()
