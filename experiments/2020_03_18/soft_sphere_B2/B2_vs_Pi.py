import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append("../../scripts")
import scipy.integrate as integrate

from softsphere import SoftSphere

if __name__=="__main__":

    colors = sns.color_palette()

    # start by plotting at different Omegas
    
    Pi = 3./2.
    Omegas = [0.1,10,1000]


    fig,axarr = plt.subplots(2,sharex=True)
    fig.set_size_inches(3.6,4)
    
    for i,Omega in enumerate(Omegas):
    
        ss = SoftSphere(Omega=Omega,Pi=Pi)

        r0 = ss.r_0()
        rs = np.linspace(r0/2,5,num=1000,endpoint=True)

        axarr[0].plot(rs,ss.w(rs),color=colors[i],
                      label=rf'$\Omega={Omega}$')
        if i == len(Omegas)-1:
            r0label=r'$r_0$'
        else:
            r0label=None
            
        axarr[0].plot(r0,ss.w(r0),'ko',ms=2,label=r0label)

    axarr[0].set_yscale('log')
    axarr[0].set_ylim(bottom=1e-3,top=max(ss.w(rs))*1000)
    axarr[0].set_ylabel(r'$w(r)$')
    axarr[0].legend(frameon=False)

    # now plot at different Pis
    
    Omega = 1
    Pis = [0.1,10,1000]
    

    for i,Pi in enumerate(Pis):
    
        ss = SoftSphere(Omega=Omega,Pi=Pi)

        r0 = ss.r_0()
        rs = np.linspace(r0/100,5,num=1000,endpoint=True)

        axarr[1].plot(rs,ss.w(rs),color=colors[i+len(Omegas)],
                      label=rf'$\Pi={Pi}$')
        if i == len(Pis)-1:
            r0label=r'$r_0$'
        else:
            r0label=None
            
        axarr[1].plot(r0,ss.w(r0),'ko',ms=2,label=r0label)

    axarr[1].set_yscale('log')
    axarr[1].set_ylim(bottom=1e-3,top=max(ss.w(rs))*1000)
    axarr[1].legend(frameon=False)
    axarr[1].set_xlabel(r'$r$')
    axarr[1].set_ylabel(r'$w(r)$')



    fig.subplots_adjust(left=0.2,right=0.9)
    fig.savefig('results/fig3.pdf')
    plt.show()
