import numpy as np
import matplotlib.pyplot as plt
        

def f(t,y,l0,cutoff):

    if t < cutoff:
        v0 = 4*(1/t**12-1/t**6)+1
    else:
        v0 = 0

    return [y[1],-3*y[1]/t+l0*y[0]-v0/t]

def jac(t,y,l0):

    return [[0,1],[l0,-3/t]]
  

from scipy.integrate import ode



t0 = 0.5
tf = 1000
dt = 0.5

l0 = 1.0
cutoff = 2**(1/6)


r = ode(f,jac).set_integrator('vode',method='bdf',with_jacobian=True)

def evaluate_ode(y0,l0,cutoff):

    r.set_initial_value(y0,t0).set_f_params(l0,cutoff).set_jac_params(l0)

    while r.successful() and r.t > tf:
        
        r.integrate(r.t+dt)

    return r.y

if __name__=="__main__":

    
    
    fig,axarr = plt.subplots(2)
    
    fig.set_size_inches(4,4)
    
    axarr[0].plot(ts,y0s)
    
    axarr[1].plot(ts,y1s)
    
    plt.show()
