import numpy as np
import scipy.special as special
import scipy.optimize as optimize
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import time
from vectorquad import quadrature_vec

class WCApotential():

    def __init__(self,epsilon=1):
        
        self.epsilon = epsilon

        return

    def V_WCA(self,r):
        
        a = 4*self.epsilon*(r**(-12)-r**(-6))+self.epsilon

        return np.where(r<2**(1./6.),a,0*a)
    
    def limited_domain_WCA_deriv(self,r):

        # this function is only valid when r <= 1 !!!
        
        return -24*self.epsilon*(2*r**(-13)-r**(-7))

    def V_WCA_prime_scalar(self,r):

        if r <= 2**(1./6.):
            a = self.limited_domain_WCA_deriv(r)

        else:
            a = 0

        return a
    
    def V_WCA_prime(self,r):

        a = self.limited_domain_WCA_deriv(r)

        return np.where(r<=2**(1./6.),a,a*0)
        

if __name__ == "__main__":

    rs = np.logspace(-0.2,0.2,num=10000,endpoint=True)

    cp = WCApotential()

    fig,axarr = plt.subplots(3,sharex=True)

    fig.set_size_inches(4,4*2)

    axarr[0].plot(rs,cp.V_WCA(rs))

    true_p = cp.V_WCA_prime(rs)
    num_p = np.gradient(cp.V_WCA(rs),rs)
    axarr[1].plot(rs,true_p,'.')
    axarr[1].plot(rs,num_p,'k-')

    axarr[2].plot(rs[1:],np.abs(true_p-num_p)[1:],'o')
    axarr[2].set_yscale('log')

    plt.show()
