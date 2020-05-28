from scipy.special import i1, k1
from scipy.integrate import quad


def i1_over12(x):

    return i1(x)/x**12

def s_1(x):

    return quad(i1_over12,0.5,x)[0]

def i1_over6(x):

    return i1(x)/x**6

def s_2(x):

    return quad(i1_over6,0.5,x)[0]

def k1_over12(x):

    return k1(x)/x**12

def q_1(x):

    return quad(k1_over12,0.5,x)[0]

def k1_over6(x):

    return k1(x)/x**6

def q_2(x):

    return quad(k1_over6,0.5,x)[0]
