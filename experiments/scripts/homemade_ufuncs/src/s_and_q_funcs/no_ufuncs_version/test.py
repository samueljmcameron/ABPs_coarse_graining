import s_and_q_funcs as funcs
import sys
import numpy as np

import matplotlib.pyplot as plt

sys.path.append('../../../generate_coeffs')

from true_functions import s_1, s_2, q_1, q_2

xs = np.linspace(0.5,30.0,num=100001,endpoint=True)

fnames = ['s_1','s_2','q_1','q_2']

chb_funcs = {'s_1' : funcs.s1_py, 's_2' : funcs.s2_py,
             'q_1' : funcs.q1_py, 'q_2' : funcs.q2_py}

true_funcs = {'s_1' : s_1, 's_2' : s_2,
             'q_1' : q_1, 'q_2' : q_2}
import time
tstart = time.time()
for f in fnames:
    
    ys = chb_funcs[f](xs)

tend = time.time()

print(f'total_time = {tend-tstart}')




#ytrues = xs*0    
#for i,x in enumerate(xs):
#    ytrues[i] = true_funcs[f](x)



    #plt.plot(xs,ys-ytrues,label=rf'${f}(x)$')


#plt.show()
