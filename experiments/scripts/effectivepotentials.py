import numpy as np

import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.special as spec
import matplotlib.pyplot as plt

from softsphere import SoftSphere

class W2potential(SoftSphere):

    def __init__(self,Omega=None,l0=1,epsilon=1,Pi=1):

        super().__init__(Omega=Omega,l0=l0,epsilon=epsilon,Pi=Pi)
        
        return


    def evaluate(self,r):
        
        return (self.V_WCA(r)
                -0.5*self.l0**2*self.w(r)**2*r**2)    


    def prime(self,r):

        w = self.w(r)

        return (self.V_WCA_prime(r)-self.l0**2*r*w**2
                -self.l0**2*r**2*w*self.w_prime(r))

class W3potential(SoftSphere):

    def __init__(self,Omega=None,l0=1,epsilon=1,Pi=1):

        super().__init__(Omega=Omega,l0=l0,epsilon=epsilon,Pi=Pi)
        
        return

    
    def evaluate(self,u,v,theta):
        
        return (-3./2.*self.l0**2*self.w(u)*self.w(v)
               *u*v*np.cos(theta))

    def prime_u(self,u,v,theta):

        return (-3./2.*self.l0**2*self.w_prime(u)*self.w(v)
                *u*v*np.cos(theta)
                -3./2.*self.l0**2*self.w(u)*self.w(v)
                *v*np.cos(theta))

    def prime_v(self,u,v,theta):

        return self.prime_u(v,u,theta)
        

class Chi2potential(SoftSphere):

    def __init__(self,Omega=None,l0=1,epsilon=1,Pi=1):
            
        super().__init__(Omega=Omega,l0=l0,epsilon=epsilon,Pi=Pi)
        
        return

    def evaluate(self,r):

        w = self.w(r)

        return self.V_WCA(r)-2*np.log(spec.i0(self.l0*w*r))
    

    def prime(self,r):

        w = self.w(r)

        return (self.V_WCA_prime(r)
                -2*self.l0*spec.ivp(0,self.l0*w*r)/spec.i0(self.l0*w*r)
                *(w+r*self.w_prime(r)))    


if __name__ == "__main__":


    fps = np.array([0,1,5,10],float)
    
    Pi = 3./2.*2**(1./3.)

    rxs = np.linspace(0.6,2,num=1000,endpoint=True)

    fig,axarr = plt.subplots(2,2,sharex=True)
    
    for i,fp in enumerate(fps):
        print(fp)
        l0 = fp*2**(1./6.)

        chi_2 = Chi2potential(l0 = l0, Pi=Pi)
        W_2 = W2potential(l0 = l0, Pi = Pi)

        W2s = W_2.evaluate(rxs)
        W2primes = W_2.prime(rxs)
        
        chi2s = chi_2.evaluate(rxs)
        chi2primes = chi_2.prime(rxs)

        label = rf'$f_p={fp}$'
        
        axarr[0][0].plot(rxs,W2s,label=label)
        axarr[0][1].plot(rxs,chi2s,label=label)

        axarr[1][0].plot(rxs,W2primes,label=label)
        axarr[1][1].plot(rxs,chi2primes,label=label)
        axarr[1][0].plot(rxs,np.gradient(W2s,rxs),'k--')
        axarr[1][1].plot(rxs,np.gradient(chi2s,rxs),'k--')
    axarr[0][0].set_ylim(-10,10)
    axarr[0][1].set_ylim(-10,10)
    axarr[1][0].set_ylim(-10,10)
    axarr[1][1].set_ylim(-10,10)
    plt.show()
