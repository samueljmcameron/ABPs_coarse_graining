def D0(fp):

    Dt = 1
    taur = 1./(3*Dt)

    return fp**2*taur/2.

def rms(fp,ts):

    Dt = 1
    taur = 1./(3*Dt)
    d = 2

    tts = ts*1e-5
    
    return 4*Dt*tts+fp**2*taur**2/(d*(d-1))*(2*d*tts/taur+np.exp(-2*d*tts/taur)-1)

def swim(fp,rho):

    Dt = 1
    taur = 1./(3*Dt)

    return rho*fp*fp*taur/2.0

if __name__=="__main__":

    import numpy as np
    import matplotlib.pyplot as plt
    import sys

    sys.path.append('../../plotting_scripts')
    from jupyterplots import JupyterPlots
    sys.path.append('../../analysis_scripts')
    from logloader import LogLoader


    figsize = JupyterPlots()

    prefix1 = 'data2/'
    tcut = -1
    fps = np.array([1,5,10,20,40,60,80,100],int)

    # density (in lj units)
    rho = '0.7'

    fig,axarr = plt.subplots(2,sharex=True,figsize=[figsize[0],figsize[0]*2])

    fp = 100


    fname = f'pressure_{fp}_{rho}'
    ll = LogLoader(prefix1 + f'log_{fp}_{rho}.lammps.log')

    ts = ll.data['Step']
    RMSs = ll.data['c_mymsdd[4]']
    Ps = ll.data['c_press']
    Ds = ll.data['v_Diff']

    
    #data = np.loadtxt(prefix1 + fname + '.txt.fixprint')


    #ts = data[:,0]
    #Ts = data[:,1]
    #Ds = data[:,2]
    #Ps = data[:,3]
    #tcut = 200000
    print(ts)
    #axarr[0].plot(ts[1:],np.gradient(RMSs,ts)[1:]/4e-5,'o',label=rf'$f_p={fp}$')
    axarr[0].plot(ts,RMSs,'o',label=rf'$f_p={fp}$')
    axarr[0].plot(ts,rms(fp,ts),'k-')
    #axarr[0].plot(ts,D0(fp)+0*ts,'k-')
    #axarr[0].plot(ts[1:],Ds[1:],'.',label=rf'$f_p={fp}$')
    axarr[1].plot(ts,Ps,'o',label=rf'$f_p={fp}$')
    axarr[1].plot(ts,swim(fp,float(rho))+ts*0,'k-',label=rf'$f_p={fp}$')

    axarr[0].set_ylabel(r'$<R^2>$')
    axarr[1].set_ylabel(r'$P_{eff}$')
    axarr[1].set_xlabel(r'$t$')
    fig.savefig('results_single_pressure/' + fname + '.pdf')
    plt.show()


