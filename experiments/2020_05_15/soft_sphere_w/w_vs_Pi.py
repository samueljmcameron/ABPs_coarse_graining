import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append("../../scripts")

from softsphere import SoftSphere

sys.path.append('../../plotting_scripts')

from jupyterplots import JupyterPlots

if __name__=="__main__":

    figsize = JupyterPlots()
    
    colors = sns.color_palette()

    # start by plotting at different Omegas
    
    Pi = 3./2.

    l0s = np.array([0.1,50.0,100.0],float)
    fig,ax = plt.subplots(figsize=[figsize[0],figsize[1]])

    
    for i,l0 in enumerate(l0s):
    
        ss = SoftSphere(l0=l0,Pi=Pi)

        r0 = ss.r_0()
        rs = np.linspace(r0/2,5,num=1000,endpoint=True)

        ax.plot(rs,ss.w(rs),color=colors[i],
                      label=rf'$\mathrm{{Pe}}={l0}$')
        if i == len(l0s)-1:
            r0label=r'$r_0$'
        else:
            r0label=None
            
        ax.plot(r0,ss.w(r0),'ko',ms=2,label=r0label)

    #ax.set_yscale('log')
    #ax.set_ylim(bottom=1e-3,top=max(ss.w(rs))*1000)
    ax.set_ylim(bottom=-max(ss.w(rs)),top=max(ss.w(rs))*2)
    ax.set_xlim(left=r0/1.4)
    ax.set_ylabel(r'$w(r)$')
    ax.set_xlabel(r'$r$')
    ax.legend(frameon=False)


    fig.subplots_adjust(left=0.2,right=0.9)
    fig.savefig('results/2020_05_15_soft_sphere.pdf')
    plt.show()
