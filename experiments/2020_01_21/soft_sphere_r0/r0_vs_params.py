import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append("../../scripts")

from softsphere import SoftSphereBase

if __name__=="__main__":

    colors = sns.color_palette()

    epsilon = 1

    l0s = np.logspace(-3,5,num=10001,endpoint=True)/2**(1./6.)

    r0s = np.copy(l0s)

    r0smalls = np.copy(l0s)
    r0larges = np.copy(l0s)

    ss = SoftSphereBase(l0s[0],epsilon,np.nan)

    for i,l0 in enumerate(l0s):

        ss.Pe = l0
        
        r0s[i] = ss.r_0()

        #if ss.a_0() <= 1:
        r0smalls[i] = ss.r_0_small_a_0()
        #else:
        #    r0smalls[i] = np.nan

        if 0.5<= ss.a_0()*(2**(7./6.)):
            r0larges[i] = ss.r_0_large_a_0()
        else:
            r0larges[i] = np.nan


    fig,axarr = plt.subplots(2,sharex=True)
    fig.set_size_inches(3.6,4)

    a0s = l0s/(24*np.sqrt(3)*epsilon)
    
    axarr[0].plot(a0s,r0s,'-',color=colors[3])
    axarr[0].set_ylabel(r'$r_0$')

    axarr[1].plot(a0s,np.abs(r0smalls-r0s)/r0s,'--',
                  label=r'$\Omega\ll1$',color=colors[1])
    axarr[1].plot(a0s,np.abs(r0larges-r0s)/r0s,'-.',
                  label=r'$\Omega\gg1$',color=colors[2])
    axarr[1].set_xlabel(r'$\Omega$')
    axarr[1].set_ylabel(r'fractional error')
    axarr[1].set_yscale('log')
    axarr[1].set_xscale('log')
    axarr[1].set_ylim(bottom=1e-11)
    axarr[1].legend(frameon=False,loc='upper left')

    fig.subplots_adjust(left=0.2,right=0.9)
    fig.savefig('results/fig2.pdf')
    plt.show()
