import numpy as np
import sys

sys.path.append('../../analysis_scripts')

from dumpfile import DumpFile

from pickle_dump import save_obj, load_obj
from spatialcorrelations import calculate_items

if __name__ == "__main__":

    rho = sys.argv[1]
    
    fps = np.array([0,1,5,10,20,40,60,80,100])

    load_prefix = 'raw_data/'

    save_prefix = 'pickled_data/'

    for fp in fps:

        fname = f'dump_{fp}_{rho}.lammpstrj'

        dfile = DumpFile(load_prefix + fname)

        data = dfile.read_dump(pstatus=True,min_step=1000000)
        save_obj(data,save_prefix+fname)

        ret_o = {}
        calculate_items(ret_o,data, min_neigh=4, cutoff=1.5, MAXnb=100,
                        nbins=2000, nbinsq=50, Pe=10, rho_0=0.60)

        dc_name = save_prefix + f'ret_o_{fp}_{rho}'

        save_obj(ret_o,dc_name)
