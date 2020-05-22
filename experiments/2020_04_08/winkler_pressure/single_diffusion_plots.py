if __name__=="__main__":

    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import sys

    sys.path.append('../../analysis_scripts')

    from pickle_dump import save_obj,load_obj
    from logloader import LogLoader
    import seaborn as sns

    import thermodata_analysis as thermo
    
    sys.path.append('../../plotting_scripts')

    from jupyterplots import JupyterPlots

    colors = sns.color_palette()
    
    fig_x,fig_y = JupyterPlots()
    
    prefix = 'data/'

    fp = int(sys.argv[1])
    rho = sys.argv[2]

    label = rf'$f^P={fp},\rho={rho}$'
    
    flog = prefix + f'log_{fp}_{rho}.lammps.log'
    ll = LogLoader(flog,remove_chunk=0,merge_data=True)


    fig,ax = plt.subplots(figsize=[fig_x,fig_y])

    ts = ll.data['Step']
    MSD = ll.data['c_mymsdd[4]']

    l1 = 'LAMMPS'
    l2 = 'Python'
    tcut = 1000000

    dt = 0.00001
    ax.plot(ts[ts>tcut]*dt,MSD[ts>tcut],'o',color=colors[0],
            label=label)
    ax.set_ylabel(r'$<R^2>$')

    Deff,Derr,yint = thermo.compute_Deff(ts,MSD)

    ax.plot(ts*dt,4*Deff*(ts*dt)+yint,'k-')
    ax.legend()

    fig.subplots_adjust(left=0.25)
    fig.savefig(f'results/diffusion_{fp}_{rho}.pdf')

    plt.show()
