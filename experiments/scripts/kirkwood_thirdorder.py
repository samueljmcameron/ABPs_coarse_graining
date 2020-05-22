import numpy as np

import scipy.integrate as integrate
import scipy.interpolate as interpolate
import scipy.special as spec
import matplotlib.pyplot as plt

from softsphere import SoftSphere

import integrate_3d

from effectivepotentials import W3potential

class Kirkwood_3(SoftSphere):

    def __init__(self,potential='W_2',Omega=None,l0=1,epsilon=1,
                 Pi=1,convert_units=True,N_distbins = 200,N_theta = 200,
                 dist_lower = 0.5,dist_upper = 8.0,
                 theta_lower = 0.0, theta_upper = np.pi):

        super().__init__(Omega=Omega,l0=l0,epsilon=epsilon,Pi=Pi)

        if convert_units:
            self.prefactor = 2**(1./6.)
        else:
            self.prefactor = 1.

        self.potential = W3potential(Omega=Omega,l0=l0,
                                     epsilon=epsilon,Pi=Pi)

        self.N_u = N_distbins
        self.N_v = N_distbins
        self.N_theta = N_theta

        self.u_a = dist_lower
        self.u_b = dist_upper
        self.v_a = dist_lower
        self.v_b = dist_upper
        self.theta_a = theta_lower
        self.theta_b = theta_upper

        self.cutoff=dist_upper
        
        return

    def partial_integrand(self,u,v,theta):

        # although integrand should have u*v*(u*W_3'(u)+v*W_3'(v)),
        # since W_3(u,v,theta) = W_3(v,u,theta) and
        # g3(u,v,theta) = g3(v,u,theta), we can just integrate one
        # of the above terms twice and then multiply the result by two,
        # hence the factor of 2 in the returned value below.
        
        return 2*u*u*v*self.potential.prime_u(u,v,theta)
    
    def g3_integral(self,g3data):

        # The three body correlation functions are measured
        # in LAMMPS using r* = r/sigma_LJ where sigma_LJ is the
        # value at which the lennard jones potential crosses zero.
        
        # Our analytic results using \tilde{r} = r/sigma, where
        # sigma corresponds to the minimum in the WCA potential.

        # so need to change units in evaluating the integral

        # will return the g3 integral in theory units (with
        # length scaled by sigma, not sigma_LJ).

        
        ans = integrate_3d.simps(g3data,self.partial_integrand,
                                 self.N_u,self.N_v,self.N_theta,
                                 self.prefactor,self.u_a,self.u_b,
                                 self.v_a,self.v_b,self.theta_a,
                                 self.theta_b,self.cutoff)

        return ans

    def third_order_correction(self,rho,g3data):


        # here rho must be in units of rmin (i.e. distance of WCA
        # potential's minimum and cutoff), NOT lj sigma units

        a = -rho**3*self.g3_integral(g3data)*np.pi/3.0

        return a


if __name__ == "__main__":

    pass
