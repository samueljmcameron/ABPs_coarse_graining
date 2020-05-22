import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../scripts')
from kirkwood import Kirkwood

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

prefix1 = '../../2020_04_08/winkler_pressure/data/'
prefix2 = '../../2020_04_08/winkler_pressure/rdf_splits/'

fps = np.array(['0','0.25','0.5','0.75','1','2','4','8'],
               str)
try:
    rho = sys.argv[1]
except:
    rho = '0.05'

    
# load in radial distribution functions for data set
rdf = RdfLoader(prefix2,Nsamples)

Pcalcs = np.empty([len(fps)],float)
Pno_activity = np.empty([len(fps)],float)
Ploads = np.empty([len(fps)],float)

tcut = 1000000

prefactor = 2**(1./6.)

for i, fp in enumerate(fps):

    # first load in pressures calculated via lammps:

    flog = prefix1 + f'log_{fp}_{rho}.lammps.log'
    ll = LogLoader(flog,remove_chunk=0,merge_data=True)

    ts = ll.data['Step']
    Ps = ll.data['c_press'][ts>tcut]
         # + thermo.ideal_swim_press(float(fp),float(rho))

    Ploads[i] = prefactor**2*np.mean(Ps)

    # next, load in the local information required from the
    # two body correlation function

    rdf = RdfLoader(prefix2,Nsamples)
    rs,g_means = rdf.read_gmeans(fp,rho,savgol=True)

    # next, rescale parameters to be in units of r_min instead of
    # sigma
    dens = float(rho)*prefactor**2
    l0 = float(fp)*prefactor

    Pi = 3./2.*prefactor**2

    kw = Kirkwood(l0 = l0, Pi=Pi)

    Pcalcs[i] = kw.pressure_correct(dens,g_means,rs)

    kwbad = Kirkwood(l0 = 0,Pi = Pi)

    Pno_activity[i] = kwbad.pressure_correct(dens,g_means,rs)



fig,axarr = plt.subplots(2,sharex=True,figsize=[figsize[0],2*figsize[1]])

fpnums = np.array(fps,float)
print(fpnums)
axarr[0].plot(fpnums,Ploads,'o',label='LAMMPS')
axarr[0].plot(fpnums,Pcalcs,'-.',label='Kirkwood')
axarr[0].plot(fpnums,Pno_activity,'-.',label=r'$f_P=0$')

axarr[0].legend(frameon=False)

axarr[1].plot(fpnums,(Pcalcs-Ploads)/Ploads,'o')

fig.savefig(f'results/errors_when_rho_is_{rho}.pdf')
plt.show()

