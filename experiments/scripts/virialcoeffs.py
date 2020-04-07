import numpy as np

import scipy.integrate as integrate
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt

from softsphere import SoftSphere

class VirialCoeffs():

    def __init__(self,Omega=None,l0=1,epsilon=1,Pi=1):

        super().__init__(Omega=Omega,l0=l0,epsilon=epsilon,Pi=Pi)
        
        return
    
    
    def W_3(self,r1,r2,theta):
        
        return (-3./2.*self.l0**2*self.w(r1)*self.w(r2)
               *r1*r2*np.cos(theta))

    def f_12_scalar(self,r):

        return np.exp(-self.W_2(r))-1
    
    def f_12(self,r,cutoff=0.8):


        r0 = r[r<=cutoff]
        r1 = r[r>cutoff]
        
        return np.concatenate((r0*0-1,self.f_12_scalar(r1)))
    
    def f_123(self,r1,r2,theta):
        
        return np.exp(-self.W_3(r1,r2,theta))-1
    
    
    def b_2_integrand(self,r,cutoff=0.9):

        if r<=cutoff:
            ans = -r
        else:
            ans = r*self.f_12_scalar(r)
            
        return ans
    
    def B_2(self):
        
        return -np.pi*integrate.quad(self.b_2_integrand,
                                     0,20)[0]


if __name__ == "__main__":

    rs = np.linspace(0,2,num=500,endpoint=True)

    l0 = 5.0
    epsilon = 10.0
    Pi = 1.0
    
    vc = VirialCoeffs(l0=l0,epsilon=epsilon,Pi=Pi)

    fig,ax = plt.subplots()

    fig.set_size_inches(4,4)

    ax.plot(rs,vc.f_12(rs))

    print(vc.B_2())

    plt.show()
