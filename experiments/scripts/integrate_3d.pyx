#!python
#cython: language_level=3


import numpy as np

from libc.math cimport cos, sqrt, floor, fabs
from libc.stdio cimport printf


cimport numpy as np


# data type of the g3 correlation function
DTYPE = np.float64

# I don't understand why, but this next line is necessary
# according to cython docs
ctypedef np.float64_t DTYPE_t




cdef double square(double x):
    """
    compute square of x
    """
    return x*x


cdef double compute_rjksquare(double u, double v, double theta):
    """
    compute distance |u-v| = u^2 + v^2 - 2*u*v*cos(theta)
    """
    return square(u) + square(v) - 2*u*v*cos(theta)


cdef double max_v(double u, double theta, double cutoff):
    """
    compute maximum value of v, given u and theta, which
    still satisfies the cutoff inequality
    
        cutoff^2 > u^2 + v^2 -2*u*v*cos(theta).

    The above inequality is satisfied for v\in(v_{-},v_{+})
    where

    v_{\pm} = u cos(theta) \pm sqrt(u^2(cos^2(theta)-1)+cutoff^2).

    So v_max = v_{+}.
    """
    
    return u*cos(theta) + sqrt(square(u)*(square(cos(theta))-1)
                               +square(cutoff))


cdef int max_N(double u, double theta, double v_a,
               double dv, double cutoff,int N_v):
    """
    Calculate the number of discretised v values (with spacing
    dv), max_N, which will still satisfy the inequality

        cutoff^2 > u^2 + v^2 -2*u*v*cos(theta),  (1)

    while also not exceeding the maximum number of bins in
    allowed, so max_N< N_v.

    Note v = v_a + (j+0.5*dv), j=0..max_N-1. Then, determine
    the maximum v which satisfies eqn (1), and call it v_max.
    Then,

    max_N = min(floor((v_max-v_a)/dv -0.5) + 1.,N_v)

    """

    cdef int N_max

    N_max = <int>(floor((max_v(u,theta,cutoff)-v_a)/dv-0.5))+1

    if N_max > N_v:
        N_max = N_v
    
    return N_max



def simps_test(int N_u,int N_v,int N_theta,
          double u_a, double u_b, double v_a, double v_b,
          double theta_a, double theta_b,double cutoff):
    """
    Checks to see if two different methods of
    indexing are equal for a flattened three dimensional array
    with a condensed form
    
    The array (which is not actually needed for this test
    function) represents a three dimensional function,
    f(u,v,theta).

    Parameters
    ----------
    N_u :
        number of unique u values in f(u,v,theta)
        discretisation.
    N_v :
        number of unique v values in f(u,v,theta)
        discretisation.
    N_theta :
        number of unique theta values in f(u,v,theta)
        discretisation.
    u_a, v_a, theta_a :
        lower limits of u, v, and theta, respectively
    u_b, v_b, theta_b :
        upper limits of u, v, and theta, respectively
    cutoff :
        cutoff in distance for u,v and |u-v|.
    
    However, it has a more complicated indexing due to the
    presence of a cutoff value, and so this function aims to
    ensure that this indexing is accessed in the correct manner.

    A naive treatment of storing values for f(u,v,theta) would be
    to take f[k + N_v * (j + N_v *i)], where 
     
            theta = theta_a + (i+0.5)*(theta_b-theta_a)/N_theta,
                u = u_a + (j+0.5)*(u_b-u_a)/N_u,
                v = v_a + (k+0.5)*(v_b-v_a)/N_v,
     
    and i=0..N_theta-1,  j=0..N_u-1, and k=0..N_v-1. With this naive
    treatment, it would be easy to index this array in an efficient
    manner, using
     
             count = 0
             for theta in thetas:
                 for u in us:
                     f[count:count+N_v] = f(u,vs,theta)
                     count += N_v.

    However, for the array that this function is looking to prototype,
    the indexing is
     
             count = 0
             for theta in thetas:
                 for u in us:
                     for v in vs:
                         if (u^2+v^2-2*u*v<cutoff^2):
                             f[count] = f(u,v,theta)
                             count += 1,
     
    and so the way to index in an efficient manner (i.e. using
    numpy slicing wherever possible) is less obvious. This
    function exists to check that I am doing the numpy slicing
    correctly.

    """
    
    cdef int i, j, k,low_ind, high_ind, count, N_max
    cdef double s, du,dv,dtheta
    cdef double u, theta
    cdef v_max

    cdef np.ndarray[DTYPE_t, ndim=2] udum = np.empty([N_u*N_v*N_theta,3],
                                                     dtype=np.float)
    cdef int i_dum
    
    du = (u_b-u_a)/N_u
    dv = (v_b-v_a)/N_v
    dtheta = (theta_b-theta_a)/N_theta

    count = 0
    for i in range(N_theta):
        
        theta = theta_a + (i+0.5)* dtheta
        
        for j in range(N_u):

            u = u_a + (j + 0.5) * du
            
            for k in range(N_v):

                v = v_a + (k + 0.5) *dv

                if (compute_rjksquare(u,v,theta)
                    < square(cutoff)):

                    udum[count,0] = theta
                    udum[count,1] = u
                    udum[count,2] = v
                    count += 1
    
    s = 0
    count = 0

    for i in range(N_theta):
        theta = theta_a + (i+0.5) * dtheta

        for j in range(N_u):

            u = u_a + (j+0.5)*du

            if (compute_rjksquare(u,v_a + 0.5*dv,theta)
                < square(cutoff)):

                N_max = max_N(u,theta,v_a,dv,cutoff,N_v)

                low_ind = count
                high_ind = low_ind + N_max

                i_dum = count
                for k in range(N_max):

                    v = v_a + (i_dum-count+ 0.5)*dv

                    if (fabs(udum[i_dum,0] - theta) > 1e-12):
                        printf("incorrect theta! theta_true = %f, "
                               "theta_calc = %f\n",udum[i_dum,0],theta)
                    if (fabs(udum[i_dum,1] - u) > 1e-12):
                        printf("incorrect u! u_true = %f, "
                               "u_calc = %f\n",udum[i_dum,1],u)
                    if (fabs(udum[i_dum,2] - v) > 1e-12):
                        printf("incorrect v! v_true = %f, "
                               "v_calc = %f\n",udum[i_dum,2],v)
                    i_dum += 1
                                
                count += N_max

    print(count)
        
    return 


def simps(np.ndarray[DTYPE_t, ndim=2] f,func,
          int N_u,int N_v,int N_theta,double prefactor,
          double u_a, double u_b, double v_a, double v_b,
          double theta_a, double theta_b,double cutoff):
    """
    Calculate the three dimensional integral with
    integrand f(u,v,theta)*func(u,v,theta).


    Parameters
    ----------
    f :
        discretised, flattened 3D array of a hypothetical
        function of three variables, u, v, and theta.
        Note that it is the fourth row which actually
        gives the function values, since the array is
        indexed as

             count = 0
             for theta in thetas:
                 for u in us:
                     for v in vs:
                         if (u^2+v^2-2*u*v<cutoff^2):
                             f[count,0] = count
                             f[count,1] = theta
                             f[count,2] = u
                             f[count,3] = v
                             f[count,4] = f(u,v,theta)
                             count += 1,

    func :
        pythonic, continuous valued function of three
        variables u, v, and theta. Note that this
        function slows everything down, since it is
        relies on python code which cannot be converted
        to C.
    N_u :
        number of unique u values in f(u,v,theta)
        discretisation.
    N_v :
        number of unique v values in f(u,v,theta)
        discretisation.
    N_theta :
        number of unique theta values in f(u,v,theta)
        discretisation.
    prefactor :
        conversion of regular units (u,v,theta) to
        the units func uses. Necessary because distances
        in f are measured in simulation units, u = u_t/sigma,
        v = v_t/sigma, theta = theta_t, but func distances
        are in analytical units of u* = u_t/r_min,
        v* = v_t/r_min, and theta* = theta_t. In general,
        r_min == sigma*prefactor which is why this conversion
        is necessary.
    u_a, v_a, theta_a :
        lower limits of u, v, and theta, respectively (all
        in simulation units)
    u_b, v_b, theta_b :
        upper limits of u, v, and theta, respectively (all
        in simulation units)
    cutoff :
        cutoff in distance for u,v and |u-v| (in simulation
        units)
    
    The three dimensional integral I care about evaluating is
    with f discretising the three-body correlation function
    g^{(3)}(u,v,theta) (with distances in LJ units), and func
    is W_3(u*,v*,theta*) (with distances in WCA units).
    """
    
    # make sure f is the correct data type
    assert f.dtype == DTYPE
    
    cdef int i, j, k,low_ind, high_ind, count, N_max
    cdef double s, du,dv,dtheta
    cdef double u, theta
    cdef v_max
    # intentionally not defining vs, as this gives a small speedup
    # to the calculation
    
    du = (u_b-u_a)/N_u
    dv = (v_b-v_a)/N_v
    dtheta = (theta_b-theta_a)/N_theta
    
    s = 0
    count = 0

    for i in range(N_theta):
        theta = theta_a + (i+0.5) * dtheta

        for j in range(N_u):

            u = u_a + (j+0.5)*du

            if (compute_rjksquare(u,v_a + 0.5*dv,theta)
                < square(cutoff)):

                N_max = max_N(u,theta,v_a,dv,cutoff,N_v)
       
                low_ind = count
                high_ind = low_ind + N_max
                
                count += N_max

                vs = f[low_ind:high_ind,3] /prefactor

                #for k in range(N_max):
                #    vs[k] = (v_a + (k + 0.5)*dv)/prefactor
                u = u /prefactor

                
                s += np.sum(f[low_ind:high_ind,4]*func(u,vs,theta))

    print(count)

    du = du/prefactor
    dv = dv/prefactor
    
    return s*du*dv*dtheta

