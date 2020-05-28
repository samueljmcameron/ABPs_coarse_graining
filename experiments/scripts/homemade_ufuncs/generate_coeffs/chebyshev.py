import numpy as np
from scipy.special import i1
import matplotlib.pyplot as plt

def chbevl(x,array,n):
    """
    Compute chebeyshev approximation to a function. Stolen from
    cephes C library (chbevl.c). Evaluates the series

    y = \sum_{i=0}^{n-1} array[i] T_i (x/2)

    of Chebyshev polynomials T_i at argument x/2

    Parameters
    ==========
    x - argument that the function takes, x\in[-2,2].

    array - coefficients of the chebyshev series, which will be
        uniquely defined for each function.

    n - number of coefficients (length of array).


    Returns
    =======
    f(x) - the value of the function at x

    Notes
    =====
    The coefficients are computed as

    array[j] = 2/n*\sum_{k=0}^{n-1} f(x_k) T_j(x_k)

    where x_k = cos( pi*(k+0.5)/n).

    """

    index = 0
    
    b0 = array[index]
    index += 1
    b1 = 0.0

    i = n-1
    
    b2 = b1
    b1 = b0
    b0 = x*b1 - b2 + array[index]
    index += 1
    i = i -1
    while (i > 0):
        b2 = b1
        b1 = b0
        b0 = x*b1 - b2 + array[index]
        index += 1
        i = i -1

    return 0.5*(b0-b2)

def c_j_one(j,N,a,b,func,mapping):
    """
    Compute coefficient j of the chebyshev approximation, i.e.
    array in 

    f(y) = \sum_{i=0}^{n-1} array[i] T_i (y/2).

    where y = 2*(2*x - b - a)/(b-a) \in[-2,2] is the mapping from
    an arbitrary range of values x for which the function is defined
    over.

    Parameters
    ==========
    j - index of the coefficient being computed

    N - total number of coefficients to be computed

    a - lower bound of valid x in f(x).

    b - upper bound of valid x in f(x).

    func - f(x)


    Returns
    =======
    array[j] - coefficient j of the chebyshev polynomials

    Notes
    =====
    The coefficients are computed as

    array[j] = 2/n*\sum_{k=0}^{n-1} f(x_k) T_j(x_k)

    where x_k = cos( pi*(k+0.5)/n).

    """

    # sum up terms in s
    s = 0
    
    for k in range(N):

        # transform from [-1,1] to the x range
        x_k =np.cos(np.pi*(k+0.5)/N)
        
        y = mapping(a,b,x_k)

        
        s += func(y)*np.cos(np.pi*j*(k+0.5)/N)

    return s*2.0/N

def map_back_direct(a,b,y_k):
    
    return 0.5*(b-a)*y_k+0.5*(b+a)

def map_back_inf(a,b,y_k):
    
    return 2*b/(y_k+1)

def map_direct(a,b,x):

    return (x-0.5*(b+a))/(0.5*(b-a))

def map_inf(a,b,x):

    return 2*b/x -1


def c_j_direct(a,b,N,func):

    return np.array([c_j_one(j,N,a,b,func,map_back_direct)
                     for j in range(N-1,-1,-1)],float)

def c_j_to_inf(b,N,func):
    
    return np.array([c_j_one(j,N,0,b,func,map_back_inf)
                     for j in range(N-1,-1,-1)],float)

def f_approx(x,func,low_bounds,upp_bounds,coeff_arrays,
             Ns):

    if x < low_bounds[0]:
        print('error! x is too small')
        return np.nan

    else:
        i = 0

        while (x > upp_bounds[i] and i < len(Ns)):

            i += 1

        if i < len(Ns):
            low = low_bounds[i]
            upp = upp_bounds[i]
            coeff = coeff_arrays[i]
            y = 2*map_direct(low,upp,x)
            f =chbevl(y,coeff,Ns[i])
        else:
            print(f'error!! x must be less than {b_C}')
            f = np.nan
    return f


if __name__ == "__main__":


    def i1_exp_inv(x):

        return np.exp(-x)*i1(x)/x

    def a_k(k):
        
        if k == 0:
            return 1
    
        s = 1
        for i in range(1,k+1):
            s *= (4-(2*i-1)**2)/(i*8)
            
        return s
    
    def i_asymp_exp_sqrt(x,nmax):

        s = 1

        for k in range(1,nmax):
            a = a_k(k)/x**k
            
            s += (-1)**k*a

        return s/np.sqrt(2*np.pi)


    def i1_exp_sqrt(x):

        if x > 800:
            return i_asymp_exp_sqrt(x,10)
        elif x > 400:
            return i_asymp_exp_sqrt(x,12)
        elif x > 200:
            return i_asymp_exp_sqrt(x,14)
        elif x > 100:
            return i_asymp_exp_sqrt(x,28)
        else:
            return i1(x)*np.exp(-x)*np.sqrt(x)


    a = 0.0
    b = 8.0

    
    N = 28

    A = c_j_direct(a,b,N,i1_exp_inv)

    B = c_j_to_inf(b,N,i1_exp_sqrt)

    print(A)
    print(B)
