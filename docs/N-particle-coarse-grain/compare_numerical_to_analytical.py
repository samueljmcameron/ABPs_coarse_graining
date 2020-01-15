import numpy as np
import matplotlib.pyplot as plt
from numerical_solve_homogeneous_ODE import evaluate_ode, f, jac
from analytical_solve_homogeneous_ODE import W, Wp


if __name__=="__main__":

    l0 = 0.1
    rf = 40
    yp0 = 1.86e-9
    ts,ys = evaluate_ode(yp0,l0,rf,False,f,jac)

    Ws = W(ts,l0)

    Wps = Wp(ts,l0)


    fig,axarr = plt.subplots(2)

    axarr[0].plot(ts,ys[:,0],'r.',label='numeric')
    axarr[0].plot(ts,Ws,'k-',label='analytic')

    axarr[1].plot(ts,ys[:,1],'r.',label='numeric')
    axarr[1].plot(ts,Wps,'k-',label='analytic')

    plt.show()
