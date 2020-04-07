import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.signal import savgol_filter
sys.path.append('../../analysis_scripts')
from pickle_dump import load_obj

Nsamples = 20
Nbins = 2000
gs = np.empty([Nbins,Nsamples],float)
fp = 0

for rho in ['0.05','0.1','0.2','0.4']:
    print(f'on rho = {rho}')
    
    for i in range(Nsamples):
        fname = 'corr_files/'
        fname += f'g_{fp}_{rho}_{i}.rdf'
        data = np.loadtxt(fname)
        print(i)
        rbins = data[:,1]
        gs[:,i] = data[:,2]


    gcorrs = np.mean(gs,axis=1)
    gstd = np.std(gs,axis=1)/np.sqrt(Nsamples)

    prefix = '../../2020_03_19/raw_data_processing/pickled_data/'

    dc_name = prefix + f'ret_o_{fp}_{rho}'

    ret_o = load_obj(dc_name)

    gmatrix = ret_o['sum_g']
    N_ft_samples = ret_o['g_cnt']

    rs = np.array(gmatrix[:,0]/N_ft_samples).flatten()
    g_fts = np.array(gmatrix[:,1]/N_ft_samples).flatten()

    fig,axarr = plt.subplots(2)

    axarr[0].plot(rbins,gcorrs,'o-')#yerr=gstd,fmt='o-')
    axarr[0].plot(rs,g_fts,'ko-')
    axarr[0].plot(rbins,savgol_filter(gcorrs,9,2),'r-')

    print(rbins[110])
    axarr[1].hist(gs[130,:],bins = 5)

    axarr[0].set_xlim(0.5,5)
    plt.show()
