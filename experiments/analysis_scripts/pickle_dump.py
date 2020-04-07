import numpy as np


from scipy.special import jv, yv, i0e,i1e

import matplotlib.pyplot as plt

import sys, os
import pickle
import fortrantools as ft
from dumpfile import DumpFile




def save_obj(obj, name,pstatus=True):
    """
    This function pickles an object and saves it to the disk.
    """
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        print ('saved data to \'%s\'.' % (name + '.pkl'))
        
        
def load_obj(name,pstatus=False):
    """
    This function loads a pickled object.
    """
    with open(name + '.pkl', 'rb') as f:
        if pstatus:
            print ('loading from \'%s\' ...' % (name + '.pkl'))
        return pickle.load(f)



if __name__ == "__main__":

    fnames = ['test/dump_80_0.6.lammpstrj',
              'test/dump_10_0.04.lammpstrj']

    fname = 'data/ABP/dump_60_0.4.lammpstrj'

    for fname in fnames:
    
        dfile = DumpFile(fname)
    
        data = dfile.read_dump(pstatus=True,min_step=1000000)
        
        save_obj(data,fname)
