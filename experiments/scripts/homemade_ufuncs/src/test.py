
import sys
import numpy as np

import nufuncs

sys.path.append('../../')

from softsphere import SoftSphere

Pe = 4.0
Pi = 3.0
epsilon= 5.0

ss = SoftSphere(Pe,epsilon,Pi)

rs = np.linspace(0.5,30,num=1000,endpoint=True)

import time as time
tstart = time.time()
old = ss.nu_1_old(rs)
tend = time.time()

nu_1_old = tend-tstart

print(f'time for nu_1_old = {nu_1_old}')

tstart = time.time()
new = ss.nu_1_new(rs)
tend = time.time()

nu_1_new = tend-tstart

print(f'time for nu_1_new = {nu_1_new}')

print(old[np.logical_not(np.isclose(old,new))])


tstart = time.time()
old = ss.nu_2_old(rs)
tend = time.time()

nu_2_old = tend-tstart

print(f'time for nu_2_old = {nu_2_old}')


tstart = time.time()
new = ss.nu_2_new(rs)
tend = time.time()

nu_2_new = tend-tstart

print(f'time for nu_2_new = {nu_2_new}')

print(old[np.logical_not(np.isclose(old,new))])
