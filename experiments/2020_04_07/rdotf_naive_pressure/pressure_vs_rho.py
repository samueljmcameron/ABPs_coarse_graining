def p_id(fp,rho,d):

    return rho*(1+fp**2/(d*(d-1)*3))

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

    fps = np.array([0,1,100],int)
    rhos = np.array([0.05,0.1,0.2,0.3,0.4,0.5,0.6],float)

    tcut = 1000000

    fig,ax = plt.subplots(figsize=[fig_x,fig_y])


    
    for i,fp in enumerate(fps):
    

        pkl_name = prefix + f'P_stuff_for_fp_equals_{fp}'

        if os.path.isfile(pkl_name+'.pkl'):
            Ps,Perrs = load_obj(pkl_name,pstatus=True)
        else:

            Ps = np.empty([len(rhos)],float)
            Perrs = np.empty([len(rhos)],float)
            
            for j,rho in enumerate(rhos):

                flog = prefix + f'log_{fp}_{rho}.lammps.log'
                ll = LogLoader(flog,remove_chunk=0,merge_data=True)


                ts = ll.data['Step']
                Press = ll.data['c_press']

                P_cuts = Press[ts>tcut]
                Ps[j] = np.mean(P_cuts)
                Perrs[j] = np.std(P_cuts)/np.sqrt(len(P_cuts))

            save_obj([Ps,Perrs],pkl_name,pstatus=True)


        ax.errorbar(rhos,Ps,
                    yerr=Perrs,fmt='-o',
                    color=colors[i],label=rf'$f^P={fp}$')

    ax.set_ylabel(r'$P$')
    ax.set_xlabel(r'$\rho$')
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.legend()
    

    fig.subplots_adjust(left=0.25)
    fig.savefig(f'results/vs_rho_press.pdf')

    plt.show()
