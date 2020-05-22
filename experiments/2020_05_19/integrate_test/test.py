import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../scripts')
from integrate_3d import simps_test



if __name__=="__main__":
    
    g3data_shape = [5983390,5]
    print(g3data_shape)

    N_u = N_v = N_theta = 200

    u_a = v_a = 0.5
    u_b = v_b = cutoff = 8.0

    theta_a = 0.0
    theta_b = np.pi
    
    simps_test(N_u,N_v,N_theta,u_a,u_b,v_a,v_b,
               theta_a,theta_b,cutoff)

