import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.append('../../analysis_scripts')

from dumpfile import DumpFile

from pickle_dump import save_obj, load_obj
from spatialcorrelations import calculate_items

if __name__ == "__main__":

    rho = sys.argv[1]
    
    fps = np.array([0])#,1,5,10,20,40,60,80,100])

    load_prefix = '../raw_data_processing/pickled_data/'

    for fp in fps:



        dc_name = load_prefix + f'ret_o_{fp}_{rho}'

        ret_o = load_obj(dc_name)

        gmatrix = ret_o['sum_g']
        Nsamples = ret_o['g_cnt']

        rs = gmatrix[:,0]/Nsamples

        gs = gmatrix[:,1]/Nsamples

        plt.plot(rs,gs)
        plt.show()
