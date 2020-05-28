import numpy as np
from scipy.special import i1
import matplotlib.pyplot as plt
import chebyshev 
from scipy.integrate import quad
from file_io import write_to_coeffs
from true_functions import s_1, s_2, q_1, q_2

"""
generate coefficients for s_1(x) where 0.5 < x < 30.0 with
absolute accuracy of 1e-12 ish 
"""


import sys

try:

    test = sys.argv[1]
except:
    test = 'no'


func_names = ['s_1','s_2','q_1','q_2']

func_defs = {'s_1' : 'int_{0.5}^x I_1(u)/u^12 du',
             's_2' : 'int_{0.5}^x I_1(u)/u^6 du',
             'q_1' : 'int_{0.5}^x K_1(u)/u^12 du',
             'q_2' : 'int_{0.5}^x K_1(u)/u^6 du' }

func_list = {'s_1' : s_1,
             's_2' : s_2,
             'q_1' : q_1,
             'q_2' : q_2}

N_alls = {'s_1' : [40,44,44],
          's_2' : [28,40,30,28,28],
          'q_1' : [40,44,44],
          'q_2' : [30,40,30]}

low_bound_alls = {'s_1' : [0.5,1.0,5.0],
                  's_2' : [0.5,1.0,5.0,20.0,25.0],
                  'q_1' : [0.5,1.0,5.0],
                  'q_2' : [0.5,1.0,5.0]}

upp_bound_alls = {'s_1' : [1.0,5.0,30.0],
                  's_2' : [1.0,5.0,20.0,25.0,30.0],
                  'q_1' : [1.0,5.0,30.0],
                  'q_2' : [1.0,5.0,30.0]}




for func_name in func_names:

    Ns = N_alls[func_name]
    
    low_bounds = low_bound_alls[func_name]
    upp_bounds = upp_bound_alls[func_name]
    
    coeff_arrays = []

    num_intervals = len(Ns)

    func = func_list[func_name]

    func_def = func_defs[func_name]
    
    for i in range(num_intervals):
        N = Ns[i]
        a = low_bounds[i]
        b = upp_bounds[i]

        coeff_arrays.append(chebyshev.c_j_direct(a,b,N,func))


    write_to_coeffs(func_name,func_def,coeff_arrays,
                    low_bounds,upp_bounds)





    if test == 'test':

        xs = np.linspace(0.7,30,num=200,endpoint=True)
        fs = xs*0

        ftrue = fs*0

        for i,x in enumerate(xs):

            fs[i] = chebyshev.f_approx(x,func,low_bounds,
                                       upp_bounds,coeff_arrays,Ns)
            ftrue[i] = func(x)


        plt.plot(xs,(fs-ftrue))
        plt.show()

        plt.plot(xs,ftrue)
        plt.show()
