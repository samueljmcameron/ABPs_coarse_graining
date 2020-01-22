import numpy as np
import scipy.special as special
import scipy.optimize as optimize
import scipy.integrate as integrate
import matplotlib.pyplot as plt


class SoftSphereBase:

    def __init__(self,l0,epsilon,Pi):

        self.l0 = l0
        self.epsilon = epsilon
        self.Pi = Pi

        return

    def a_0(self):

        return self.l0/(12*np.sqrt(3)*self.epsilon)

    def r_0_large_a_0(self):

        q = self.a_0()**(-1/13)
        p = self.a_0()**(-6/13)
        f7 = 13-7*p
        g4 = 13-4*p
        
        return q*(1+1/14*f7/g4*(1-np.sqrt(1+28*p*g4/f7**2)))

    def r_0_small_a_0(self):

        return 1+1/21*(1-np.sqrt(1+7*self.a_0()))

    
    def r_0_approx(self):

        if self.epsilon == 0 or self.l0 == 0:
            return 1

        a0 = self.a_0()

        if a0 < 1:

            r0 = self.r_0_small_a_0()

        else:

            r0 = self.r_0_large_a_0()

        return r0

    
    def forcebalance(self,r):
        
        return self.a_0() - r**(-13) + r**(-7)

    def r_0(self):
        
        r0 = self.r_0_approx()

        if r0 != 1:

            r0 = optimize.newton(self.forcebalance,r0)

        return r0

    def limited_domain_WCA_deriv(self,r):

        # this function is only valid when r <= 1 !!!
        
        return -12*self.epsilon*(r**(-13)-r**(-7))

    def nu_1_integrand(self,s):

        return (special.k1(np.sqrt(self.Pi)*s)
                *self.limited_domain_WCA_deriv(s)*s)

    def nu_2_integrand(self,s):

        return (special.i1(np.sqrt(self.Pi)*s)
                *self.limited_domain_WCA_deriv(s)*s)


    def k_1(self):

        return -0.5*integrate.quad(self.nu_1_integrand,
                                   self.r0,1)[0]

    def k_2(self,dumk1):

        temp = -1/(3*np.sqrt(self.Pi))

        temp -= dumk1*special.ivp(1,np.sqrt(self.Pi)*self.r0,1)

        return temp/special.kvp(1,np.sqrt(self.Pi)*self.r0,1)

    def c_2(self,dumk2):

        return dumk2 - 0.5*integrate.quad(self.nu_2_integrand,self.r0,1)[0]

    def constants(self):

        dumk1 = self.k_1()
        dumk2 = self.k_2(dumk1)
        dumc2 = self.c_2(dumk2)

        return dumk1, dumk2, dumc2

    

class SoftSphere(SoftSphereBase):

    def __init__(self,Omega=None,l0=1,epsilon=1,Pi=1):

        if Omega != None:
            epsilon = 1
            l0 = Omega*12*np.sqrt(3)*epsilon

        super().__init__(l0,epsilon,Pi)

        self.r0 = self.r_0()

        self.k1,self.k2,self.c2 = self.constants()

        return

    def nonvec_nu_1(self,r):

        # can only take scalars as input

        return self.k1+0.5*integrate.quad(self.nu_1_integrand,
                                          self.r0,r)[0]

    def nonvec_nu_2(self,r):

        # can only take scalars as input

        return self.k2-0.5*integrate.quad(self.nu_2_integrand,
                                          self.r0,r)[0]

    def nu_1(self,r):

        # only valid for r <= 1 !!!!!

        return np.vectorize(self.nonvec_nu_1)(r)


    def nu_2(self,r):

        # only valid for r <= 1 !!!!!

        return np.vectorize(self.nonvec_nu_2)(r)


    def w_minus(self,r):

        if self.r0 == 1:

            return r*0
        
        return (special.i1(np.sqrt(self.Pi)*r)*self.nu_1(r)/r
                +special.k1(np.sqrt(self.Pi)*r)*self.nu_2(r)/r)

    def w_plus(self,r):

        return self.c2*special.k1(np.sqrt(self.Pi)*r)/r

    def w(self,r):

        return np.where(r>=1,self.w_plus(r),self.w_minus(r))
                        

if __name__ == "__main__":

    rs = np.linspace(1,40,num=500,endpoint=True)

    l0 = 1.0
    epsilon = 10.0
    Pi = 1.0

    ss = SoftSphere(l0,epsilon,Pi)

    fig,ax = plt.subplots()

    fig.set_size_inches(4,4)

    ax.plot(rs,ss.w(rs))

    print(min(ss.w(rs)))

    plt.show()
