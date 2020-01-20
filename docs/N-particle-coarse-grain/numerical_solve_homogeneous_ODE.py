import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode        

def f(t,y,l0,usepotential):

    cutoff = 2**(1/6)
    
    if t < cutoff and usepotential:
        v0 = 4*(1/t**12-1/t**6)+1
    else:
        v0 = 0

    return [y[1],-3*y[1]/t+l0*y[0]+v0/t]

def jac(t,y,l0):

    return [[0,1],[l0,-3/t]]
  



def evaluate_ode(yp0,l0,tf,usepotential,f,jac):

    r = ode(f,jac).set_integrator('vode',method='bdf',with_jacobian=True)
    t0 = 1
    num = 500

    yf = [0,yp0]
    
    r.set_initial_value(yf,tf).set_f_params(l0,usepotential).set_jac_params(l0)

    dt = (tf-t0)/(num-1)
    
    ys = np.empty([num,2],float)

    ts = np.empty([num],float)

    count = num-1
    
    ys[count,:] = yf

    ts[count] = tf
    
    
    while r.successful() and count > 0:

        count -= 1
        
        r.integrate(r.t-dt)

        ys[count,:] = r.y
        ts[count] = r.t
        
    return ts,ys

def Wboundary(W,Wp):

    return -1-3*Wp-3*W

if __name__=="__main__":


    l0 = 0.1
    yp0 = -1.86e-9
    tf = 40
    ts,ys = evaluate_ode(yp0,l0,tf,False,f,jac)

    print(ys[0,0],ys[0,1],Wboundary(ys[0,0],ys[0,1]))
    
    fig,axarr = plt.subplots(2)
    
    fig.set_size_inches(4,4)
    
    axarr[0].plot(ts,ys[:,0])
    
    axarr[1].plot(ts,ys[:,1])

    
    
    plt.show()
