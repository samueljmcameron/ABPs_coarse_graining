import numpy as np
import scipy.special as special
import matplotlib.pyplot as plt

def c_2(l0):

    K1prime = (-special.kn(0,np.sqrt(l0))
               -1/np.sqrt(l0)*special.kn(1,np.sqrt(l0)))

    return 1/(3*np.sqrt(l0)*K1prime)

def W0(x,l0):

    return special.kn(1,np.sqrt(l0)*x)/x

def W(x,l0):

    return c_2(l0)*W0(x,l0)


def Wp(x,l0):

    c2 = c_2(l0)

    return -c2*(W(x,l0)/(x*x) + np.sqrt(l0)/x
                *(special.kn(0,np.sqrt(l0)*x)
                  +1/(np.sqrt(l0)*x)*special.kn(1,np.sqrt(l0)*x)))
                

if __name__ == "__main__":

    xs = np.linspace(1,20,num=100,endpoint=True)

    l0 = 0.1
    ys = W(xs,l0)
    yps = Wp(xs,l0)

    fig,axarr = plt.subplots(2)

    fig.set_size_inches(4,4)

    axarr[0].plot(xs,ys)

    axarr[1].plot(xs,yps)

    plt.show()
