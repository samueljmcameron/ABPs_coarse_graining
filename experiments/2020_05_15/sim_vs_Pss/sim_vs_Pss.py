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



fig,ax = plt.subplots(figsize=[figsize[0],figsize[1]])

fpnums = np.array(fps,float)
print(fpnums)
ax.plot(fpnums,Ploads,'o',label='simulation')
ax.plot(fpnums,Pcalcs,'-.',label=r'$P_{ss}$')
#ax.plot(fpnums,Pno_activity,'-.',label=r'$f_P=0$')
ax.set_xlabel(r'$\mathrm{Pe}$')
ax.set_ylabel(r'$P-P_{\mathrm{id}}$')

ax.legend(frameon=False)

fig.subplots_adjust(left=0.25,right=0.9,bottom=0.2,top=0.9)
fig.savefig(f'results/2020_05_15sim_vs_Pss_{rho}.pdf')
plt.show()

