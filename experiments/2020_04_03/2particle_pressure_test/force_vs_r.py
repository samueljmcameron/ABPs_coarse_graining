if __name__=="__main__":

    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import sys

    sys.path.append('../../analysis_scripts')

    from dumpfile import DumpFile
    from pickle_dump import save_obj,load_obj
    import force_calculations as fc

    sys.path.append('../../plotting_scripts')

    from jupyterplots import JupyterPlots

    fig_x,fig_y = JupyterPlots()
    
    prefix = 'data/'

    roundoff = 5e-6

    fp = int(sys.argv[1])

    fdump = prefix + f'dump_{fp}.lammpstrj'

    if os.path.isfile(fdump + '.pkl'):
        dumpdata = load_obj(fdump,pstatus=True)
    else:
        duf = DumpFile(fdump,voronoi_flag=False,cg_flag=False)
        dumpdata = duf.read_dump()
        save_obj(dumpdata,fdump,pstatus=True)

    fig,axarr = plt.subplots(6,sharex=True,figsize=[fig_x,fig_y*6])
    N = len(dumpdata)
    ncount = N

    atom_num = 2
    fxs = np.empty([ncount],float)
    fys = np.empty([ncount],float)
    xs = np.empty([ncount],float)
    ys = np.empty([ncount],float)

    muxs = np.empty([ncount],float)
    muys = np.empty([ncount],float)


    for i,data in enumerate(dumpdata[:ncount+1]):

        coords = data['coords']
        x,y = coords.T
        box = data['box']
        idmask = (data['id'] == atom_num)


        xs[i] = x[0]-x[1]
        ys[i] = y[0]-y[1]
        muxs[i] = data['mux'][idmask][0]
        muys[i] = data['muy'][idmask][0]

        fx = data['fx'][idmask][0]
        fy = data['fy'][idmask][0]
        fz = data['fz'][idmask][0]
        z = data['z'][0]
        fxs[i] = fx
        fys[i] = fy
        print(xs[i],ys[i])

    fs = np.sqrt(fxs**2+fys**2)

    # calculate forces from scratch
    fljx,fljy,flj,fx_err,fy_err,f_err = fc.f_and_df(xs,ys,muxs,muys,
                                                    fp=fp,
                                                    roundoff=roundoff)
    rs = np.sqrt(xs**2+ys**2)

    slabel = 'LAMMPS'  # data output in dump file
    clabel = 'Python'  # calculated in this script

    # plot magnitude of force vs r
    axarr[0].plot(rs,fs,'o',label=slabel)
    axarr[0].plot(rs,flj,'k-',label=clabel)
    axarr[0].set_ylabel(r'$|\bm{F}_i(r)|$')
    axarr[0].legend(loc='best')

    # plot error between calculated and output force
    axarr[1].plot(rs,np.abs(flj-fs),'o',label=slabel)
    axarr[1].plot(rs,np.abs(f_err),'k-',label=clabel)
    axarr[1].set_ylabel(r'$|\delta \bm{F}|$')
    axarr[1].set_yscale('log')

    # plot x-component of force vs r
    axarr[2].plot(rs,fxs,'o')
    axarr[2].plot(rs,fljx,'k-')
    axarr[2].set_ylabel(r'$(\bm{F}_i(r))_x$')

    # plot error for x-component of force
    axarr[3].plot(rs,np.abs(fljx-fxs),'o',label=slabel)
    axarr[3].plot(rs,np.abs(fx_err),'k-',label=clabel)
    axarr[3].set_ylabel(r'$\delta(\bm{F}_i(r))_x$')
    axarr[3].set_yscale('log')

    # plot y-component of force vs r
    axarr[4].plot(rs,fys,'o')
    axarr[4].plot(rs,fljy,'k-')
    axarr[4].set_ylabel(r'$(\bm{F}_i(r))_y$')

    # plot error for y-component of force
    axarr[5].plot(rs,np.abs(fljy-fys),'o',label=slabel)
    axarr[5].plot(rs,np.abs(fy_err),'k-',label=clabel)
    axarr[5].set_ylabel(r'$\delta(\bm{F}_i(r))_y$')
    axarr[5].set_yscale('log')
    axarr[5].set_xlabel(r'$r$')

    
    fig.subplots_adjust(left=0.25,right=0.95,top=0.95,hspace=0.05)
    fig.savefig(f'results/f_vs_r_{fp}.pdf')
    plt.show()

