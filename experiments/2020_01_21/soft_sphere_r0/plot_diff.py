import numpy as np
import matplotlib.pyplot as plt

data_old = np.loadtxt('data_oldscaling.txt')
data_new = np.loadtxt('data_newscaling.txt')

Omega_old = data_old[0,:]
r0_old = data_old[1,:]
rsmall_old = data_old[2,:]
rlarge_old = data_old[3,:]

Omega_new = data_new[0,:]
r0_new = data_new[1,:]
rsmall_new = data_new[2,:]
rlarge_new = data_new[3,:]

print(Omega_old[np.logical_not(np.isclose(Omega_old,2*Omega_new*2**(1./6.)))])

print(r0_old[np.logical_not(np.isclose(r0_old,r0_new/2**(1./6.)))])

print(rsmall_old[np.logical_not(np.isclose(rsmall_old,rsmall_new/2**(1./6.)))])

print(rlarge_old[np.logical_not(np.isclose(rlarge_old,rlarge_new/2**(1./6.),
                                           equal_nan=True))])
