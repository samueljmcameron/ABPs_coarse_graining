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

    fps = np.array(['0','0.25','0.5','0.75','1','2','4','8','100'])


    tcut = 1000000

    fig,ax = plt.subplots(figsize=[fig_x,fig_y])

    dt = 0.00001
    
    for i,fpstr in enumerate(fps):

        fp = float(fpstr)
        
        if fpstr == '100':
            rhos = np.array([0.05,0.1,0.2,0.6],float)
        else:
            rhos = np.array([0.05,0.1,0.2,0.3,0.4,0.5,0.6],float)
            
        pkl_name = prefix + f'Deff_stuff_for_fp_equals_{fpstr}'

        if os.path.isfile(pkl_name+'.pkl'):
            Ds,Derrs = load_obj(pkl_name,pstatus=True)
        else:

            Ds = np.empty([len(rhos)],float)
            Derrs = np.empty([len(rhos)],float)
            
            for j,rho in enumerate(rhos):

                flog = prefix + f'log_{fpstr}_{rho}.lammps.log'
                ll = LogLoader(flog,remove_chunk=0,merge_data=True)

                ts = ll.data['Step']
                MSD = ll.data['c_mymsdd[4]']
                Ds[j],Derrs[j],yint = thermo.compute_Deff(ts,MSD)

            save_obj([Ds,Derrs],pkl_name,pstatus=True)


        ax.errorbar(rhos,Ds/thermo.ideal_Deff(fp),
                    yerr=Derrs/thermo.ideal_Deff(fp),fmt='-o',
                    color=colors[i+1],label=rf'$f^P={fpstr}$')

    ax.set_ylabel(r'$D_{\mathrm{eff}}$')
    ax.set_xlabel(r'$\rho$')
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.legend()
    

    fig.subplots_adjust(left=0.25)
    fig.savefig(f'results/vs_rho_diffusion.pdf')

    plt.show()
