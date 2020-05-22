import numpy as np
import matplotlib.pyplot as plt
import sys
import seaborn as sns
sys.path.append("../../scripts")

from effectivepotentials import W3potential

sys.path.append('../../plotting_scripts')

from jupyterplots import JupyterPlots

if __name__=="__main__":

    figsize = JupyterPlots()
    
    colors = sns.color_palette()

    # start by plotting at different Omegas
    
    Pi = 3./2.

    l0s = np.array([0.1],float)
    fig,ax = plt.subplots(figsize=[figsize[0],figsize[1]])

    r0 = 0.9
    Npoints=400
    for i,l0 in enumerate(l0s):
    
        ss = W3potential(l0=l0,Pi=Pi)


        r0new = ss.r_0()
        if r0new < r0:
            r0 = r0new
        rs = np.linspace(r0/1.1,2,num=Npoints,endpoint=True)

        R1,R2 = np.meshgrid(rs,rs)
        Z = np.empty([Npoints,Npoints],float)

        for j,r in enumerate(rs):
            Z[j,:] = ss.evaluate(r,rs,np.pi/3.0)

        CS = ax.contourf(R1,R2,Z,cmap=plt.cm.bone,origin='lower')

    #ax.set_yscale('log')

    #ax.set_ylim(bottom=1e-3,top=max(ss.w(rs))*1000)
    #ax.set_ylim(bottom=-max(ss.w(rs)),top=max(ss.w(rs))*2)

    cbar = fig.colorbar(CS)
    ax.set_ylabel(r'$r_{12}$')
    ax.set_xlabel(r'$r_{13}$')
    ax.legend(frameon=False)
    ax.set_title(rf'$W_3(r_{{12}},r_{{13}},\theta=\pi/3;\mathrm{{Pe}}={l0})$',
                 fontsize=10)


    fig.subplots_adjust(left=0.2,right=0.9,bottom=0.2)
    fig.savefig(f'results/2020_05_15_W_3_Pe_{l0}.pdf')
    plt.show()
