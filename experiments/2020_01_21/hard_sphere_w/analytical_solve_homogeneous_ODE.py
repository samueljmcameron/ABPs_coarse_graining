import numpy as np
import scipy.special as special
import matplotlib.pyplot as plt

def c_2(Pi):

    K1prime = (-special.kn(0,np.sqrt(Pi))
               -1/np.sqrt(Pi)*special.kn(1,np.sqrt(Pi)))

    return -1/(3*np.sqrt(Pi)*K1prime)

def W0(x,Pi):

    return special.kn(1,np.sqrt(Pi)*x)/x

def W(x,Pi):

    return c_2(Pi)*W0(x,Pi)


def Wp(x,Pi):

    c2 = c_2(Pi)

    return -c2*(W(x,Pi)/(x*x) + np.sqrt(Pi)/x
                *(special.kn(0,np.sqrt(Pi)*x)
                  +1/(np.sqrt(Pi)*x)*special.kn(1,np.sqrt(Pi)*x)))
                

if __name__ == "__main__":

    xs = np.linspace(1,20,num=100,endpoint=True)

    Pi = 0.1
    ys = W(xs,Pi)
    yps = Wp(xs,Pi)

    fig,axarr = plt.subplots(2)

    fig.set_size_inches(4,4)

    axarr[0].plot(xs,ys)

    axarr[1].plot(xs,yps)

    plt.show()
