import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append("../../scripts")

from softsphere import SoftSphere

if __name__=="__main__":

    colors = sns.color_palette()

    # start by plotting at different Omegas
    
    Pi = 1./2**(1./3.)
    Omegas = np.array([0.1,10,1000],float)*2**(7./6.)


    fig,axarr = plt.subplots(2,sharex=True)
    fig.set_size_inches(3.6,4)
    
    for i,Omega in enumerate(Omegas):
    
        ss = SoftSphere(Omega=Omega,Pi=Pi,miniter=50,maxiter=100,
                        int_err=1e-5)

        r0 = ss.r_0()
        if r0 > 2**(1./6.):
            print (f'r0 = {r0}')
        rs = np.linspace(r0/2,5,num=2000,endpoint=True)/2**(1./6.)
        
        axarr[0].plot(rs,ss.w(rs),color=colors[i],
                      label=rf'$\Omega={Omega:.2e}$')
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
    
    Omega = 1*2**(7./6.)
    Pis = np.array([0.1,10,1000],float)/2**(1./3.)
    

    for i,Pi in enumerate(Pis):
    
        ss = SoftSphere(Omega=Omega,Pi=Pi)

        r0 = ss.r_0()
        rs = np.linspace(r0/100,5,num=1000,endpoint=True)/2**(1./6.)

        axarr[1].plot(rs,ss.w(rs),color=colors[i+len(Omegas)],
                      label=rf'$\Pi={Pi:.2e}$')
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
