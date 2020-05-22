import scipy.stats as stats
import numpy as np
import numpy.fft as fft
import scipy.optimize as optim


def standard_error_block_method(x):


    sigma2s = []
    sigma2errs = []
    Ts = []
    
    n = len(x)

    # find highest power of two which is less than n

    p = int(np.log2(n))
    n_p = int(2**p)
    print(np)
    
    x_p = x[n-n_p:]
    print(len(x_p))

    x_pav = np.mean(x_p)
    
    c0 = np.dot(x_p-x_pav,x_p-x_pav)/n_p

    sigma2s.append(c0/(n_p-1))
    sigma2errs.append(np.sqrt(2/(n_p-1))*c0/(n_p-1))

    count = 0

    Ts.append(0.5*(2**count-1))
    
    while n_p > 2:
        count += 1
        x_p = 0.5*(x_p[::2] + x_p[1::2])
        x_pav = np.mean(x_p)
        n_p = n_p//2
        c0 = np.dot(x_p-x_pav,x_p-x_pav)/n_p
        sigma2s.append(c0/(n_p-1))
        sigma2errs.append(np.sqrt(2/(n_p-1))*c0/(n_p-1))
        Ts.append(0.5*(2**count-1))

    return np.array(sigma2s),np.array(sigma2errs),np.array(Ts)
    
    

def standard_error_naive_estimate(T,x,corrs = None):
    """
    Estimate the standard error of a correlated time series.
    See

    Flyvberg and Peterson. "Error Estimates on Averages of
    Correlated Data". J Chem Phys, 91, 461 (1989).

    for details. This implementation is their equation 14.
    """
    cts = corrs
    try:
        len(cts)
    except:
        cts = autocorrelation(x)


    n = len(x)

    dum = np.sum([(1-i/n)*cts[i] for i in range(1,T+1)])
    dum = np.real(dum)

    ans = (cts[0] + 2*dum)/(n-2*T-1+T*(T+1)/n)

    return ans,cts

def exp_fit(ys,tlim=50,p0=30.):

    t0s = np.linspace(0,tlim,num=tlim,endpoint=False)
    ys_fit = ys[:tlim]
    def f(x,a):
        return np.exp(-x/a)
    popt,pcov = optim.curve_fit(f,t0s,ys_fit,
                                p0=[p0])
    return popt,pcov,t0s,f(t0s,*popt)

def autocorrelation(x):

    """
    Compute (biased) estimator of autocorrelation,

    c_t = 1/(n-t) sum_{k=1}^{n-t}(x_k-x_av)(x_{k+t}-x_av)
    
    thanks to Paul Panzer:
    taken from https://stackoverflow.com/questions/47850760/
    using-scipy-fft-to-calculate-autocorrelation-of-a-signal-
    gives-different-answer

    """

    # first, shift input data to change the sum of the
    # discrete fourier transform to go from -n//2 to n//2-1
    
    xp = fft.ifftshift(x-np.mean(x))
    n, = xp.shape
    # zero pad data to ensure that there is no artificial
    # correlation due to endpoints
    xp = np.r_[xp[:n//2],np.zeros_like(xp),xp[n//2:]]
    f = fft.fft(xp)
    p = np.abs(f)**2
    pi = fft.ifft(p)

    # return function, where n-t is given by
    # np.arange(n//2)[::-1]+n//2
    denom = np.linspace(0,n//2,num=n//2,endpoint=False)[::-1]+1+n//2
    if n % 2 != 0:
        # the calculation is not exact for odd n, but this shift
        # in the denominator gets it a bit closer.
        denom += 1
    return pi[:n//2]/denom
    

def compute_Deff(ts,MSD,dt=1e-5,d=2,tcut = 1000000):

    if tcut != None:
        MSD = MSD[ts>tcut]
        ts = ts[ts>tcut]
    
    Deff,yint,rval,pval,stderr = stats.linregress(ts*dt,MSD)

    Deff = Deff/(2.*d)
    return Deff,stderr,yint

def ideal_Deff(fp,d=2,Dr = 3.):

    return 1 + fp**2/(d*(d-1)*Dr)

def ideal_swim_press(fp,rho,d=2,Dr = 3.):

    return rho*fp**2/(d*(d-1)*Dr)


def ideal_pressure(fp,rho,d=2,Dr=3.):

    return rho+ideal_swim_press(fp,rho,d=d,Dr=Dr)


if __name__ == "__main__":

    from logloader import LogLoader
    import matplotlib.pyplot as plt
    
    #inpath = '../2020_04_07/rdotf_naive_pressure/data/'
    inpath = '../2020_04_08/winkler_pressure/data/'

    fname = inpath + 'log_100_0.2.lammps.log'

    ll = LogLoader(fname,remove_chunk=0,merge_data=True)

    ts = ll.data['Step']
    press = ll.data['c_press']

    tcut = 1000000
    press = press[ts>=tcut]

    Ts = np.linspace(0,5000,num=5001,endpoint=True,dtype=int)
    sigma2s = np.empty([len(Ts)],float)

    cts = None
    for i,T in enumerate(Ts):
        
        sigma2s[i],cts = standard_error_naive_estimate(T,press,corrs = cts)

    print(sigma2s[-20:])
    plt.plot(Ts,sigma2s)

    s2s,s2errs,T2s = standard_error_block_method(press)

    
    plt.errorbar(T2s,s2s,yerr=s2errs,fmt='k-')
    
    plt.show()

    

    """
    Below is code to test the autocorrelation function:



    def c(t,x):

        n = len(x)
        print(n)
        xav = np.mean(x)

        x_k = (x[:n-t]-xav)

        x_kpt = (x[t:]-xav)

        return 1/(n-t)*np.dot(x_k,x_kpt)

    from logloader import LogLoader
    import matplotlib.pyplot as plt
    
    inpath = '../2020_04_08/winkler_pressure/data/'

    fname = inpath + 'log_1_0.4.lammps.log'

    ll = LogLoader(fname,remove_chunk=0,merge_data=True)

    ts = ll.data['Step']
    press = ll.data['c_press']

    tcut = 1000000
    press = press[ts>=tcut]
    ts = ts[ts>=tcut]


    tplots = ts[:len(ts)//2]
    corrs = np.empty([len(tplots)],float)
    for i in range(len(tplots)):
        corrs[i] = c(i,press)

        
    #plt.plot(corrs,'r.')

    plt.plot(autocorrelation(press)[:len(corrs)]-corrs,'k-')
    
    plt.show()
    """
