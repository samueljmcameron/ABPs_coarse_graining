import numpy as np
import matplotlib.pyplot as plt
import sys

import time

sys.path.append('../../scripts')
from kirkwood import Kirkwood
from kirkwood_thirdorder import Kirkwood_3
sys.path.append('../../analysis_scripts')
from rdfloader import RdfLoader
from logloader import LogLoader
import thermodata_analysis as thermo

sys.path.append('../../plotting_scripts')
from jupyterplots import JupyterPlots

def gexps(rexps,g2spline):
    return (np.max(g2spline)-1)*np.exp(-(rexps-1)*3.6)+1

figsize = JupyterPlots()
Nsamples = 20

prefix_pressure = '../../2020_05_10/long_time_data_pressure/data/'
prefix_rdf = '../../2020_05_18/twobody_running_avs/data/'
prefix_g3 = '../../2020_05_18/threebody_running_avs/data/'

fps = np.array(['0','0.25','0.5','0.75','1.0','2.0','4.0'],str)

#fps = np.array(['2.0'],str)#,'4'],str)
               
try:
    rho = sys.argv[1]
except:
    rho = '0.05'

    

Pcalcs_rdf = np.empty([len(fps)],float)
Pcalcs_total = np.empty([len(fps)],float)
Pno_activity = np.empty([len(fps)],float)
Ploads = np.empty([len(fps)],float)

tcut = 1000000

prefactor = 2**(1./6.)

for i, fp in enumerate(fps):

    # first load in pressures calculated via lammps:

    flog = prefix_pressure + f'log_{fp}_{rho}.lammps.log'
    ll = LogLoader(flog,remove_chunk=0,merge_data=True)

    ts = ll.data['Step']
    Ps = ll.data['c_press'][ts>tcut]
         # + thermo.ideal_swim_press(float(fp),float(rho))

    Ploads[i] = prefactor**2*np.mean(Ps)

    print(f"LAMMPS pressure = {Ploads[i]}")

    # next, rescale parameters to be in units of r_min instead of
    # sigma
    dens = float(rho)*prefactor**2
    l0 = float(fp)*prefactor

    Pi = 3./2.*prefactor**2

    # second order correction to the pressure

    # first, load in the local information required from the
    # two body correlation function

    rdf = np.loadtxt(prefix_rdf + f'g2_av_{fp}_{rho}.rdf',skiprows=4)

    rs = rdf[:,1]

    g2s = rdf[:,2]

    # then calculate
    
    kw = Kirkwood(l0 = l0, Pi=Pi)

    Pcalcs_rdf[i] = kw.pressure_correct(dens,g2s,rs)

    print(f"second order correction = {Pcalcs_rdf[i]}")
    
    # third order correction to the pressure

    # first, load in the local three body data
    tstart = time.time()
    g3data = np.loadtxt(prefix_g3 + f'g3bod_av_{fp}_{rho}.g3bod',skiprows=4)
    tend = time.time()
    print(f"Time taken to load g3data when fp = {fp},"
          f" rho = {rho}: {tend-tstart}")


    kw3 = Kirkwood_3(l0=l0, Pi=Pi)

    tstart = time.time()    
    delP = kw3.third_order_correction(dens,g3data)
    tend = time.time()
    print(f"Time taken to compute three body term when fp = {fp},"
          f" rho = {rho}: {tend-tstart}")
    
    print(f"change from second to third order = {delP}")
    
    Pcalcs_total[i] = (Pcalcs_rdf[i]
                       + delP)
    print(f"third order correction = {Pcalcs_total[i]}")

    kwbad = Kirkwood(l0 = 0,Pi = Pi)

    Pno_activity[i] = kwbad.pressure_correct(dens,g2s,rs)



fig,axarr = plt.subplots(2,sharex=True,figsize=[figsize[0],2*figsize[1]])

fpnums = np.array(fps,float)
print(fpnums)
axarr[0].plot(fpnums,Ploads,'o',label='LAMMPS')
axarr[0].plot(fpnums,Pcalcs_rdf,'.',label=r'$\rho^2$')
axarr[0].plot(fpnums,Pcalcs_total,'.',label=r'$\rho^3$')
axarr[0].plot(fpnums,Pno_activity,'.',label=r'$f_P=0$')

axarr[0].legend(frameon=False)

axarr[1].plot(fpnums,(Pcalcs_rdf-Ploads)/Ploads,'ro')
axarr[1].plot(fpnums,(Pcalcs_total-Ploads)/Ploads,'ko')

fig.savefig(f'results/errors_when_rho_is_{rho}.pdf')
plt.show()

