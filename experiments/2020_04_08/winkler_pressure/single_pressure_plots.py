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


    fig,axarr = plt.subplots(4,sharex=False,figsize=[fig_x,4.5*fig_y])

    ts = ll.data['Step']
    Press = ll.data['c_press']+thermo.ideal_swim_press(fp,float(rho))

    l1 = 'LAMMPS'
    l2 = 'Python'
    if float(rho)>0.5 and fp >90:
        tcut = 15000000
    else:
        tcut = 1000000
    axarr[0].plot(ts[ts>tcut],Press[ts>tcut],'o',color=colors[0],
                  label=label)

    axarr[0].set_ylabel(r'$P$')
    axarr[0].set_xlabel('time step')
    axarr[0].legend()


    autocorr = thermo.autocorrelation(Press[ts>tcut])
    axarr[1].plot(autocorr[:50])

    popt,pcov,t0s,fs = thermo.exp_fit(autocorr,tlim=50,
                                      p0=5.)
    print(popt,pcov)

    tau_r = popt[0]
    #axarr[1].plot(t0s,fs,'k--',label=rf'$\tau_r={tau_r:.2e}$')
    axarr[1].plot(t0s,t0s**(-1))
    axarr[1].set_ylabel('Autocorrelation')
    axarr[1].set_xlabel(r'$\tau$')
    axarr[1].set_xscale('log')
    axarr[1].set_yscale('log')
    axarr[1].legend()

    axarr[2].hist(Press[ts>tcut],bins=20)
    axarr[2].set_xlabel(r'$P$')
    axarr[2].set_ylabel('count')


    axarr[3].hist(Press[ts>tcut][::10*int(np.ceil(tau_r))],bins=20)
    axarr[3].set_xlabel(r'$P(\Delta t>10\tau_r)$')
    axarr[3].set_ylabel('count')


    fig.subplots_adjust(left=0.25,hspace=0.5)
    fig.savefig(f'results/pressure_{fp}_{rho}.pdf')

    plt.show()
