import scipy.integrate as integrate

from scipy.special import roots_legendre
import numpy as np

import warnings

class AccuracyWarning(Warning):
    pass

def _cached_roots_legendre(n):
    """
    Cache roots_legendre results to speed up calls of the fixed_quad
    function.
    """
    if n in _cached_roots_legendre.cache:
        return _cached_roots_legendre.cache[n]

    _cached_roots_legendre.cache[n] = roots_legendre(n)
    return _cached_roots_legendre.cache[n]


_cached_roots_legendre.cache = dict()

def fixed_quad_vec(func, a, bs, args=(), n=5):
    """
    Compute a definite integral using fixed-order Gaussian quadrature.
    Integrate `func` from `a` to `bs` using Gaussian quadrature of
    order `n`. Note that bs is an array.
    Parameters
    ----------
    func : callable
        A Python function or method to integrate (must accept vector inputs).
        If integrating a vector-valued function, the returned array must have
        shape ``(..., len(x))``.
    a : float
        Lower limit of integration.
    bs : array
        Upper limits of integration.
    args : tuple, optional
        Extra arguments to pass to function, if any.
    n : int, optional
        Order of quadrature integration. Default is 5.
    Returns
    -------
    val : array (of size len(bs))
        Gaussian quadrature approximation to the integral
    none : None
        Statically returned value of None
    See Also
    --------
    quad : adaptive quadrature using QUADPACK
    dblquad : double integrals
    tplquad : triple integrals
    romberg : adaptive Romberg quadrature
    quadrature : adaptive Gaussian quadrature
    romb : integrators for sampled data
    simps : integrators for sampled data
    cumtrapz : cumulative integration for sampled data
    ode : ODE integrator
    odeint : ODE integrator
    Examples
    --------
    >>> from scipy import integrate
    >>> f = lambda x: x**8
    >>> integrate.fixed_quad(f, 0.0, 1.0, n=4)
    (0.1110884353741496, None)
    >>> integrate.fixed_quad(f, 0.0, 1.0, n=5)
    (0.11111111111111102, None)
    >>> print(1/9.0)  # analytical result
    0.1111111111111111
    >>> integrate.fixed_quad(np.cos, 0.0, np.pi/2, n=4)
    (0.9999999771971152, None)
    >>> integrate.fixed_quad(np.cos, 0.0, np.pi/2, n=5)
    (1.000000000039565, None)
    >>> np.sin(np.pi/2)-np.sin(0)  # analytical result
    1.0
    """
    x, w = _cached_roots_legendre(n)
    x = np.real(x)
    y = np.multiply.outer((bs-a),(x+1)/2.0)+np.multiply.outer(a,0*x+1)
    return (bs-a)/2.0 * np.sum(w*func(y, *args), axis=1), None


def quadrature_vec(func, a, bs, args=(), tol=1.49e-8, rtol=1.49e-8, maxiter=50,
                   miniter=1):
    """
    Compute a definite integral using fixed-tolerance Gaussian quadrature.
    Integrate `func` from `a` to `b` using Gaussian quadrature
    with absolute tolerance `tol`.
    Parameters
    ----------
    func : function
        A Python function or method to integrate.
    a : float
        Lower limit of integration.
    b : float
        Upper limit of integration.
    args : tuple, optional
        Extra arguments to pass to function.
    tol, rtol : float, optional
        Iteration stops when error between last two iterates is less than
        `tol` OR the relative change is less than `rtol`.
    maxiter : int, optional
        Maximum order of Gaussian quadrature.
    vec_func : bool, optional
        True or False if func handles arrays as arguments (is
        a "vector" function). Default is True.
    miniter : int, optional
        Minimum order of Gaussian quadrature.
    Returns
    -------
    val : float
        Gaussian quadrature approximation (within tolerance) to integral.
    err : float
        Difference between last two estimates of the integral.
    See also
    --------
    """
    if not isinstance(args, tuple):
        args = (args,)

    val = np.inf
    err = np.inf
    maxiter = max(miniter+1, maxiter)
    for n in range(miniter, maxiter+1):
        newval = fixed_quad_vec(func, a, bs, (), n)[0]
        err = abs(newval-val)
        oldval = val
        val = newval

        if np.max(err) < tol or np.any(err < rtol*abs(val)):
            break
    else:
        warnings.warn(
            f"maxiter ({maxiter}) exceeded."
            f"Latest difference = {err} at bs = {bs[abs(newval-oldval)>tol]}",
            AccuracyWarning)
    return val, err


if __name__ == "__main__":

    import time
    
    def f(x):

        return x**2

    n = 4
    
    x,w = _cached_roots_legendre(n)

    bs = np.linspace(0,5,num=1000,endpoint=True)
    #b = 4
    a = 1
    #y = (b-a)*(x+1)/2.0 + a
    y = np.multiply.outer((bs-1),(x+1)/2.)+np.multiply(a,x*0+1)


    t_start = time.time()
    
    intresults = quadrature_vec(f,a,bs)[0]

    t_fast = time.time()-t_start
    
    intslow= np.empty([len(bs)],float)

    t_start = time.time()
    
    for i,b in enumerate(bs):

        intslow[i] = integrate.quadrature(f,a,b)[0]

    t_slow = time.time()-t_start
    print(t_fast,t_slow)


    import matplotlib.pyplot as plt

    plt.plot(bs,np.abs(intresults-intslow),'.')
    plt.show()
