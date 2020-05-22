if __name__=="__main__":

    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import sys

    sys.path.append('../../analysis_scripts')

    from dumpfile import DumpFile
    from pickle_dump import save_obj,load_obj
    from logloader import LogLoader
    import fortrantools as ft
    import force_calculations as fc
    import seaborn as sns
    
    sys.path.append('../../plotting_scripts')

    from jupyterplots import JupyterPlots

    colors = sns.color_palette()
    
    fig_x,fig_y = JupyterPlots()
    
    prefix = 'data/'


    rho = sys.argv[1]

    fps = np.array([0,1,10,100,1000],int)

    tcut = 1000000

    Ps = np.empty([len(fps)],float)
    Perrs = np.empty([len(fps)],float)
    
    for i,fp in enumerate(fps):

        flog = prefix + f'log_{fp}_{rho}.lammps.log'
        ll = LogLoader(flog,remove_chunk=0,merge_data=True)


        ts = ll.data['Step']
        Press = ll.data['c_press']

        P_cuts = Press[ts>tcut]
        Ps[i] = np.mean(P_cuts)
        Perrs[i] = np.std(P_cuts)/np.sqrt(len(P_cuts))

    fig,ax = plt.subplots(figsize=[fig_x,2*fig_y])
    ax.errorbar(fps,Ps,yerr=Perrs,fmt='-o',color=colors[0],
                label=rf'$\rho={rho}$')

    ax.set_ylabel(r'$P$')
    ax.set_xlabel(r'$f^P$')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()
    

    fig.subplots_adjust(left=0.25)
    fig.savefig(f'results/vs_fp_press_{rho}.pdf')

    plt.show()
