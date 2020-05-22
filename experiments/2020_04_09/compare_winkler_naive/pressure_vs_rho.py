def p_id(fp,rho,d):

    return rho*(1+fp**2/(d*(d-1)*3))

def swim_add(fp,rho,d):
    return rho*fp**2/(d*(d-1)*3)

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
    
    prefix1 = '../../2020_04_07/rdotf_naive_pressure/data/'
    prefix2 = '../../2020_04_08/winkler_pressure/data/'

    fps = np.array([1,100],int)
    rhos_winkler = np.array([0.05,0.1,0.2,0.6],float)

    tcut = 1000000

    fig,axarr = plt.subplots(2,sharex=True,figsize=[fig_x,2*fig_y])

    rhos_naive = np.array([0.05,0.1,0.2,0.3,0.4,0.5,0.6],float)

    mask_naive = (np.isin(rhos_naive,rhos_winkler))
    print(mask_naive)    
    for i,fp in enumerate(fps):
    

        pkl_press_naive = prefix1 + f'P_stuff_for_fp_equals_{fp}'
        pkl_diff_naive = prefix1 + f'Deff_stuff_for_fp_equals_{fp}'
        pkl_press_winkler = prefix2 + f'P_stuff_for_fp_equals_{fp}'
        pkl_diff_winkler = prefix2 + f'Deff_stuff_for_fp_equals_{fp}'

        Ps_naive,Perrs_naive = load_obj(pkl_press_naive,pstatus=True)
        Ps_winkler,Perrs_winkler = load_obj(pkl_press_winkler,pstatus=True)
        Deffs_naive,Derrs_naive = load_obj(pkl_diff_naive,pstatus=True)
        Deffs_winkler,Derrs_winkler = load_obj(pkl_diff_winkler,pstatus=True)


        Ps_winkler += swim_add(fp,rhos_winkler,2)

        Ps_naive = Ps_naive[mask_naive]

        Deffs_naive = Deffs_naive[mask_naive]
        Derrs_naive = Derrs_naive[mask_naive]
        
        print(Deffs_naive)
        print(Deffs_winkler)

        print(Derrs_naive,Derrs_winkler)

        press_c = np.abs(Ps_naive-Ps_winkler)/Ps_winkler

        diff_c = np.abs(Deffs_naive-Deffs_winkler)/Deffs_winkler

        axarr[0].plot(rhos_winkler,press_c,'-o',
                      color=colors[i],label=rf'$f^P={fp}$')
        axarr[1].plot(rhos_winkler,diff_c,'-o',
                      color=colors[i],label=rf'$f^P={fp}$')

    axarr[0].set_ylabel(r'$|P_{wink}-P_{naive}|/P_{wink}$')
    axarr[1].set_ylabel(r'$|D_{wink}-D_{naive}|/D_{wink}$')
    axarr[1].set_xlabel(r'$\rho$')
    axarr[0].set_yscale('log')
    axarr[1].set_yscale('log')
    #ax.set_xscale('log')
    axarr[0].legend()
    

    fig.subplots_adjust(left=0.25)
    fig.savefig(f'results/vs_rho_comparison.pdf')

    plt.show()
