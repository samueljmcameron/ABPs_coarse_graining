import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append('../../scripts')
from effectivepotentials import W2potential,Chi2potential

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

fp = sys.argv[1]
rho = sys.argv[2]

# load in radial distribution functions for data set
rdf = RdfLoader(prefix2,Nsamples)

rbins,gs = rdf.read_gs(fp,rho)

fig,axarr = plt.subplots(2,4,sharex=True,figsize=[figsize[0]*4,figsize[1]*2])

axarr[0,0].plot(rbins,gs[:,0])
axarr[0,1].plot(rbins,gs[:,1])
axarr[0,2].plot(rbins,gs[:,2])
axarr[0,3].plot(rbins,gs[:,3])
axarr[1,0].plot(rbins,gs[:,4])
axarr[1,1].plot(rbins,gs[:,5])
axarr[1,2].plot(rbins,gs[:,6])
axarr[1,3].plot(rbins,gs[:,7])
axarr[1,3].set_xlim(left=0.5)


# average rdf
fig,axarr = plt.subplots(1,2,figsize=[figsize[0]*2,figsize[1]],
                         sharex=True)

rbins,gcorrs,gstd = rdf.read_gmeans(fp,rho)

axarr[0].plot(rbins,gcorrs)
rbins,gsmooth = rdf.read_gmeans(fp,rho,savgol=True)
axarr[1].plot(rbins,gsmooth)


rs = np.linspace(0.75,10,num=10000,endpoint=True)

g2spline = rdf.g2_spline(rs,gsmooth,rbins)

axarr[0].set_xlim(0.5,10)
axarr[0].plot(rs,g2spline)
axarr[1].plot(rs,g2spline)

# average rdf shifted to be in units of r_min
fig,axarr = plt.subplots(1,2,figsize=[figsize[0]*2,figsize[1]],
                         sharex=True)

rbins,gcorrs,gstd = rdf.read_gmeans(fp,rho)

axarr[0].plot(rbins*2**(-1./6.),gcorrs)
rbins,gsmooth = rdf.read_gmeans(fp,rho,savgol=True)
axarr[1].plot(rbins*2**(-1./6.),gsmooth)

g2spline = rdf.g2_spline(rs,gsmooth,rbins*2**(-1./6.))

axarr[0].set_xlim(0.5,10)
axarr[0].plot(rs,g2spline)
axarr[1].plot(rs,g2spline)
rexps = rs[rs>1]



axarr[1].plot(rexps,gexps(rexps,g2spline))

fig,axarr = plt.subplots(1,2,figsize=[figsize[0]*2,figsize[1]])


l0 = float(fp)*2**(1./6.)
Pi = 3./2.*2**(1./3.)
w2 = W2potential(l0=l0,Pi=Pi)

axarr[0].plot(rs,w2.evaluate(rs))
axarr[1].plot(rs,w2.prime(rs))

fig,axarr = plt.subplots(1,2,figsize=[figsize[0]*2,figsize[1]])


chi2 = Chi2potential(l0=l0,Pi=Pi)

axarr[0].plot(rs,chi2.evaluate(rs))
axarr[1].plot(rs,chi2.prime(rs))

fig,axarr = plt.subplots(1,2,figsize=[figsize[0]*2,figsize[1]])

stitches = np.where(rs<1,g2spline,gexps(rs,g2spline))

axarr[0].plot(rs,rs*rs*w2.prime(rs)*stitches)
axarr[0].plot(rs,rs*rs*w2.prime(rs)*g2spline)
axarr[1].plot(rs,rs*rs*chi2.prime(rs)*stitches)
axarr[0].plot(rs,rs*rs*chi2.prime(rs)*g2spline)


plt.show()
