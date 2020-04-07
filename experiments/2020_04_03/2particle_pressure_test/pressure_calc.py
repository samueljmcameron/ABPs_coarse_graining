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
    
    sys.path.append('../../plotting_scripts')

    from jupyterplots import JupyterPlots

    fig_x,fig_y = JupyterPlots()
    
    prefix = 'data/'

    fp = int(sys.argv[1])

    fdump = prefix + f'dump_{fp}.lammpstrj'
    
    if os.path.isfile(fdump + '.pkl'):
        dumpdata = load_obj(fdump,pstatus=True)
    else:
        duf = DumpFile(fdump,voronoi_flag=False,cg_flag=False)
        dumpdata = duf.read_dump()
        save_obj(dumpdata,fdump,pstatus=True)

    flog = prefix + f'log_{fp}.lammps.log'
    ll = LogLoader(flog,remove_chunk=None,merge_data=False)


    N = len(dumpdata)
    ncount = N

    roundoff=2e-5

    Ps = np.empty([ncount],float)
    Perrs = np.empty([ncount])

    for i,data in enumerate(dumpdata[:ncount]):
            
        Ps[i] = 0
        Perrs[i] = 0
        
        coords = data['coords']
        x,y = coords.T
        box = data['box']
        Area = box[0]*box[1]
        ids = data['id']
        num_atoms = len(ids)
        tmplist = ft.dmatrix(coords,box,1.1224,100,
                             num_atoms,2)

        dmatrix,distance_x,distance_y = tmplist[:3]
        num_nb,list_nb = tmplist[3:]
        
        for atom_num in ids:

            idmask = (ids == atom_num)
            dxs = distance_x[atom_num-1,:]
            dys = distance_y[atom_num-1,:]

            other_atoms = np.logical_not(idmask)
            dxs = dxs[other_atoms]
            dys = dys[other_atoms]
            
            mux = data['mux'][idmask][0]
            muy = data['muy'][idmask][0]

            fx = data['fx'][idmask][0]
            fy = data['fy'][idmask][0]
            x = data['x'][idmask][0]
            y = data['y'][idmask][0]

            # calculate forces from scratch
            c_fx,c_fy,c_f,fx_err,fy_err,f_err = fc.f_and_df(dxs,
                                                            dys,
                                                            mux,muy,fp=fp,
                                                            roundoff=roundoff)
            expect_fxerr = np.abs(fx_err)+0.5*roundoff
            expect_fyerr = np.abs(fy_err)+0.5*roundoff

            total_fxerr = np.abs(c_fx-fx)
            total_fyerr = np.abs(c_fy-fy)

            if total_fxerr>expect_fxerr and i !=0:
                print(f'issue with fx force on atom {atom_num}!'
                      f'at step {i}:')
                print(f'fx = {fx}, but c_fx = {c_fx}')
                print(f'error expected = {expect_fxerr} '
                      f'but measured error = {total_fxerr}')
            if total_fyerr>expect_fyerr and i !=0:
                print(f'issue with fy force on atom {atom_num}!'
                      f'at step {i}:')
                print(f'fy = {fy}, but c_fy = {c_fy}')
                print(f'error expected = {expect_fyerr}, '
                      f'but measured error = {total_fyerr}')

            Ps[i] += x*c_fx+y*c_fy
            xerr = 0.5*x*roundoff
            yerr = 0.6*y*roundoff
            
            Perrs[i] = x*(fx_err)+y*(fy_err)+xerr*c_fx+yerr*c_fy

        Ps[i] = Ps[i]/(2*Area)
    

    fig,axarr = plt.subplots(2,sharex=True,figsize=[fig_x,2*fig_y])

    ts = ll.data[0]['Step']
    Press = ll.data[0]['c_press']

    l1 = 'LAMMPS'
    l2 = 'Python'
    
    axarr[0].plot(ts,Press,'o',label=l1)
    axarr[0].plot(ts,Ps,'k-',label=l2)
    axarr[0].set_ylabel(r'$P$')
    axarr[0].legend()
    
    axarr[1].plot(ts,np.abs(Press-Ps),'o',label=l1)
    axarr[1].plot(ts,np.abs(Perrs),'k-',label=l2)
    axarr[1].set_xlabel('timestep')
    axarr[1].set_ylabel(r'$|\delta P|$')
    axarr[1].set_yscale('log')

    fig.savefig(f'results/pressure_{fp}.pdf')

    plt.show()
