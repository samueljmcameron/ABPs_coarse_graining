import numpy as np
import matplotlib.pyplot as plt
from numerical_solve_homogeneous_ODE import evaluate_ode, f, jac
from analytical_solve_homogeneous_ODE import W, Wp


if __name__=="__main__":

    l0 = 0.1
    rf = 40
    yp0 = 1.86e-9
    ts,ys = evaluate_ode(yp0,l0,rf,False,f,jac)

    Ws = W(ts,l0)

    Wps = Wp(ts,l0)


    fig,axarr = plt.subplots(2,sharex=True)

    fig.set_size_inches(3.6,5)

    axarr[0].plot(ts,ys[:,0],'r.',label='numeric')
    axarr[0].plot(ts,Ws,'k-',label='analytic')
    axarr[0].set_xlabel(r'$r$')
    axarr[0].set_ylabel(r'$W_h(r)$')
    axarr[0].text(20,-0.2,r'$\Pi=0.1$')

    axarr[0].text(-0.05, 1.1, '(a)', transform=axarr[0].transAxes,
                  fontsize=10, fontweight='bold', va='top', ha='right')
    
    axarr[1].plot(ts,ys[:,1],'r.',label='numeric')
    axarr[1].plot(ts,Wps,'k-',label='analytic')
    axarr[1].set_xlabel(r'$r$',fontsize=10)
    axarr[1].set_ylabel(r'$W_h^{\prime}(r)$',fontsize=10)
    axarr[1].legend(frameon=False,fontsize=10)
    axarr[1].text(-0.05, 1.1, '(b)', transform=axarr[1].transAxes,
                  fontsize=10, fontweight='bold', va='top', ha='right')


    fig.subplots_adjust(left=0.2,right=0.9)
    fig.savefig('Figures/fig1.pdf')
    
    plt.show()
