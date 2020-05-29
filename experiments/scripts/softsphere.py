import numpy as np
import scipy.special as special
import scipy.optimize as optimize
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import time
from vectorquad import quadrature_vec
from conservativepotential import WCApotential

from homemade_ufuncs.src import nufuncs

class SoftSphereBase(WCApotential):

    def __init__(self,Pe,epsilon,Pi):

        super().__init__(epsilon=epsilon)
        
        self.Pe = Pe
        self.Pi = Pi

        return

    def a_0(self):

        return self.Pe/(24*np.sqrt(3)*self.epsilon)

    def r_0_large_a_0(self):

        q = (2**(7./6.)*self.a_0())**(-1/13)
        p = (2**(7./6.)*self.a_0())**(-6/13)
        f7 = 13-7*p
        g4 = 13-4*p
        
        return 2**(1./6.)*q*(1+1/14*f7/g4*(1-np.sqrt(1+28*p*g4/f7**2)))

    def r_0_small_a_0(self):

        return 2**(1./6.)*(1+1/21*(1-np.sqrt(1+7*(2**(7./6.)*self.a_0()))))

    
    def r_0_approx(self):

        if self.epsilon == 0 or self.Pe == 0:
            return 2**(1./6.)

        a0 = self.a_0()

        if a0 < 2**(-7./6.):

            r0 = self.r_0_small_a_0()

        else:

            r0 = self.r_0_large_a_0()

        return r0

    
    def forcebalance(self,r):
        
        return self.a_0() - 2*r**(-13) + r**(-7)

    def r_0(self):

        r0 = self.r_0_approx()

        if abs(r0-2**(1./6.))>1e-14:

            r0 = optimize.newton(self.forcebalance,r0)

        return r0

    def nu_1_integrand(self,s):

        return (special.k1(np.sqrt(self.Pi)*s)
                *self.limited_domain_WCA_deriv(s)*s)

    def nu_2_integrand(self,s):

        return (special.i1(np.sqrt(self.Pi)*s)
                *self.limited_domain_WCA_deriv(s)*s)


    def k_1(self):

        return -0.5*integrate.quad(self.nu_1_integrand,
                                   self.r0,2**(1./6.))[0]

    def k_2(self,dumk1):

        temp = -1/(3*np.sqrt(self.Pi))

        temp -= dumk1*special.ivp(1,np.sqrt(self.Pi)*self.r0,1)

        return temp/special.kvp(1,np.sqrt(self.Pi)*self.r0,1)

    def c_2(self,dumk2):

        return dumk2 - 0.5*integrate.quad(self.nu_2_integrand,
                                          self.r0,2**(1./6.))[0]

    def constants(self):

        dumk1 = self.k_1()
        dumk2 = self.k_2(dumk1)
        dumc2 = self.c_2(dumk2)

        return dumk1, dumk2, dumc2

    

class SoftSphere(SoftSphereBase):

    def __init__(self,Omega=None,Pe=1,epsilon=1,Pi=1,
                 miniter=1,maxiter=50,int_err=1e-6):

        if Omega != None:
            epsilon = 1
            Pe = Omega*24*np.sqrt(3)*epsilon

        super().__init__(Pe,epsilon,Pi)

        self.r0 = self.r_0()

        self.k1,self.k2,self.c2 = self.constants()

        print(self.k1,self.k2,self.c2)

        self.miniter = miniter
        self.maxiter=maxiter
        self.int_err = int_err

        return

    def nu_1_prime(self,r):

        return 0.5*self.nu_1_integrand(r)

    def nu_2_prime(self,r):

        return -0.5*self.nu_2_integrand(r)


    def nu_1_old(self,r):

        return self.k1+0.5*quadrature_vec(self.nu_1_integrand,
                                          self.r0,r,
                                          miniter=self.miniter,
                                          maxiter=self.maxiter,
                                          tol=self.int_err,
                                          rtol=self.int_err)[0]
    
    def nu_2_old(self,r):

        return self.k2-0.5*quadrature_vec(self.nu_2_integrand,
                                          self.r0,r,
                                          miniter=self.miniter,
                                          maxiter=self.maxiter,
                                          tol=self.int_err,
                                          rtol=self.int_err)[0]
    
    def nu_1_new(self,r):
        
        return self.k1 + 0.5*self.epsilon*nufuncs.k1_V_int(r,
                                                           self.Pi,
                                                           self.r0)

    def nu_2_new(self,r):
        return self.k2 - 0.5*self.epsilon*nufuncs.i1_V_int(r,
                                                           self.Pi,
                                                           self.r0)


    def test_bounds(self,r):
        
        r = np.atleast_1d(r)
        
        flag = False
        if (np.sqrt(self.Pi)*self.r0 < 0.5 or
            np.sqrt(self.Pi)*self.r0 > 30.0):
            print(f'product sqrt(Pi)*r0= {np.sqrt(self.Pi)*self.r0} '
                  ' is out of bounds (should be between 0.5 and 30.0).')
            flag = True
        elif (np.sqrt(self.Pi)*r[0] < 0.5 or
              np.sqrt(self.Pi)*r[0] > 30.0):

            print(f'product sqrt(Pi)*r[0] = {np.sqrt(self.Pi)*r[0]} '
                  ' is out of bounds (should be between 0.5 and 30.0).')
            flag = True
        elif (np.sqrt(self.Pi)*r[-1] < 0.5 or
              np.sqrt(self.Pi)*r[-1] > 30.0):
            
            print(f'product sqrt(Pi)*r[-1] = {np.sqrt(self.Pi)*r[-1]} '
                  ' is out of bounds (should be between 0.5 and 30.0).')
            flag = True

        return flag
    
    def nu_1(self,r,flag=False):
        
        if flag == True:
            print('Reverting to slow calculation of nu_1(r)')
            return self.nu_1_old(r)
        else:
            return self.nu_1_new(r)

    def nu_2(self,r,flag=False):
                
        if flag == True:
            print('Reverting to slow calculation of nu_1(r)')
            return self.nu_1_old(r)
        else:
            return self.nu_1_new(r)
    
    def w_minus(self,r,flag=None):

        if abs(self.r0-2**(1./6.)) < 1e-15:

            return r*0

        if flag == None:

            flag = self.test_bounds(r)
        
        return (special.i1(np.sqrt(self.Pi)*r)*self.nu_1(r,flag=flag)/r
                +special.k1(np.sqrt(self.Pi)*r)*self.nu_2(r,flag=flag)/r)

    def w_minus_prime(self,r,flag = None):

        if abs(self.r0-2**(1./6.)) < 1e-15:

            return r*0

        if flag == None:
            flag = self.test_bounds(r)
        
        sPi = np.sqrt(self.Pi)

        return (sPi*special.ivp(1,sPi*r,1)
                *self.nu_1(r,flag=flag)/r
                +special.i1(sPi*r)*self.nu_1_prime(r)/r
                +sPi*special.kvp(1,sPi*r,1)
                *self.nu_2(r,flag=flag)/r
                +special.k1(sPi*r)*self.nu_2_prime(r)/r
                -self.w_minus(r,flag=flag)/r)

    
    def w_plus(self,r):
        
        return self.c2*special.k1(np.sqrt(self.Pi)*r)/r

    def w_plus_prime(self,r):

        sPi = np.sqrt(self.Pi)

        return (self.c2*sPi*special.kvp(1,sPi*r,1)/r
                - self.w_plus(r)/r)
                

    def w(self,r):
        
        return np.where(r>=2**(1./6.),self.w_plus(r),self.w_minus(r))

    
    def w_prime(self,r):

        a = np.where(r>=2**(1./6.),self.w_plus_prime(r),
                     self.w_minus_prime(r))

        return a
        

if __name__ == "__main__":



    rs = np.linspace(0.8,10,num=5000,endpoint=True)*2**(1./6.)

    Pe = 1.0/2**(1./6.)
    epsilon = 1.0
    Pi = 1.0/2**(1./3.)

    ss = SoftSphere(Pe,epsilon,Pi)

    print(ss.r0)
    fig,axarr = plt.subplots(3,sharex=True)

    fig.set_size_inches(4,4*2)

    axarr[0].plot(rs,ss.w(rs))
    axarr[0].plot(2**(1./6.),ss.w(2**(1./6.)),'ko')
    axarr[0].plot(2**(1./6.),ss.w_plus(2**(1./6.)),'r.')
    axarr[0].plot(2**(1./6.),ss.w_minus(2**(1./6.)),'b.')
    
    true_p = ss.w_prime(rs)
    num_p = np.gradient(ss.w(rs),rs)
    axarr[1].plot(rs,true_p,'.')
    axarr[1].plot(rs,num_p,'k-')

    axarr[2].plot(rs,np.abs(true_p-num_p))
    axarr[2].set_yscale('log')

    plt.show()
