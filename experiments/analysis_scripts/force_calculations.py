import numpy as np


def lj_cut(xij,yij,epsilon=1,cutoff=1.1224):
    """
    evaluate the lennard jones potential to get the x,y and
    total force on a particle.

    xij is an array with N entries, ones for every atom pair
    containing atom i

    yij is an array with N entries, ones for every atom pair
    containing atom i.

    returns an N length array, with the lennard jones
    force exerted on particle i from each of the N
    other atoms.
    """

    rij = np.sqrt(xij**2+yij**2)
    ans = 24*epsilon/rij*(2*(1./rij)**(12)-(1./rij)**6)

    ans_x = -xij/rij*ans
    ans_y = -yij/rij*ans
    
    ans_x = np.where(rij>cutoff,0,ans_x)
    ans_y = np.where(rij>cutoff,0,ans_y)
        
    return ans_x,ans_y

def dfvar_lj(var,xij,yij,xerr,yerr,epsilon=1,
             cutoff=1.1224):
    """
    evaluate the infinitesimal of fvar, where
    var = 'x' (var = 'y') is the x (y) component of
    the force. the infinitesimal of fvar is defined as

    /     dfvar = (partial fvar /partial xij)*xerr
    /            + (partial fvar /partial yij)*yerr

    xij is an array with N entries, ones for every atom pair
    containing atom i OR for a pair of two atoms at N times.

    yij is an array with N entries, ones for every atom pair
    containing atom i OR for a pair of two atoms at N times.

    xerr is an array with N entries, ones for every atom pair
    containing atom i OR for a pair of two atoms at N times.

    yerr is an array with N entries, ones for every atom pair
    containing atom i OR for a pair of two atoms at N times.

    returns an array of the errors for each of the N entries
    of the array. If the xij, yij array contain N atoms, then
    to get the force on atom i you need to sum this output.
    

    """
    
    rij2 = xij**2+yij**2
    
    dfdx = 96*epsilon*xij**2/rij2*(-7/rij2**3+2)/rij2**4
    dfdy = 96*epsilon*yij**2/rij2*(-7/rij2**3+2)/rij2**4
    
    if var == 'x':
        dfdx += 24*epsilon/rij2**4*(2/rij2**3-1)
    elif var == 'y':
        dfdy += 24*epsilon/rij2**4*(2/rij2**3-1)

    ans = dfdx*xerr+dfdy*yerr
    ans = np.where(rij2>cutoff*cutoff,0,ans)

    return ans


def ftotal(xij,yij,mux,muy,fp=0,epsilon=1,
           cutoff=1.1224,single_pair=True):
    """
    Compute the total energy (with activity) of the
    particle pair at several timesteps.

    If single_pair = True, then the xij and yij
    are assumed to be the distances between a single
    particle pair at multiple timesteps. Otherwise,
    it is assumed that xij,yij are data for multiple
    particles at a single time step. This matters
    when applying the active forces, as the two cases
    are no longer equivalent. When single_pair=False,
    mux and muy must be scalars.

    If single_pair=True, returns three arrays:
    
    fx -force on particle i in the x direction
    fy -force on particle i in the y direction
    f - magnitude of force on particle i

    all oflength len(xij).

    If single_pair=False, returns three scalars for
    the forces, fx,fy,f.

    """
    
    fx,fy = lj_cut(xij,yij,epsilon=epsilon,cutoff=cutoff)

    if single_pair:
        fx = fx + mux*fp
        fy = fy + muy*fp

    else:
        fx = np.sum(fx)+mux*fp
        fy = np.sum(fy)+muy*fp
    
    f = np.sqrt(fx**2+fy**2)

    return fx,fy,f


def f_and_df(x12,y12,mux,muy,fp,epsilon=1,roundoff=5e-7,
             cutoff=1.1224,single_pair=True):
    """
    evaluate the infinitesimal of f = np.sqrt(fx**2+fy**2), 
    where the infinitesimal df of f is defined as

    /     df = fx*dfx/f + fy*dfy/f 

    it is assumed that the roundoff error on each of the
    numbers is 1e-15

    x1 is a list of the x-coordinates of the atom which the
    force will be evaluated for (one for each time step)

    x2 is a list of the x-coordinates of the other atom
    (one for each time step).

    returns dfx,dfy, df

    """
    x12err = x12*roundoff
    
    y12err = y12*roundoff

    muxerr = mux*roundoff
    muyerr = muy*roundoff

    fx_err = dfvar_lj('x',x12,y12,x12err,y12err,
                      epsilon=epsilon,cutoff=cutoff)

    fy_err = dfvar_lj('y',x12,y12,x12err,y12err,
                      epsilon=epsilon,cutoff=cutoff)
    
    if not single_pair:

        fx_err = np.sum(fx_err)+muxerr*fp
        fy_err = np.sum(fy_err)+muyerr*fp

    else:

        fx_err = fx_err+muxerr*fp

        fy_err = fy_err+muyerr*fp

    fx,fy,f = ftotal(x12,y12,mux,muy,fp=fp,epsilon=epsilon,
                     cutoff=cutoff,single_pair=single_pair)
    
    f_err = fx/f*fx_err+fy/f*fy_err

    return fx,fy,f,fx_err,fy_err,f_err

if __name__=="__main__":


    # sample Python calculation of forces from lammps data
    
    import matplotlib.pyplot as plt
    import os
    import sys

    from dumpfile import DumpFile
    from pickle_dump import save_obj,load_obj


    sys.path.append('../plotting_scripts')

    from jupyterplots import JupyterPlots

    fig_x,fig_y = JupyterPlots()
    
    prefix = '../2020_04_03/2particle_pressure_test/data/'

    fp = 0

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

    atom_num = 1
    fxs = np.empty([ncount],float)
    fys = np.empty([ncount],float)
    xs = np.empty([ncount],float)
    ys = np.empty([ncount],float)


    muxs = np.empty([ncount],float)
    muys = np.empty([ncount],float)

    for i,data in enumerate(dumpdata[:ncount]):

        coords = data['coords']
        x,y = coords.T
        box = data['box']
        idmask = (data['id'] == atom_num)


        xs[i] = x[1]-x[0]
        ys[i] = y[1]-y[0]

        muxs[i] = data['mux'][idmask][0]
        muys[i] = data['muy'][idmask][0]

        fx = data['fx'][idmask][0]
        fy = data['fy'][idmask][0]
        fz = data['fz'][idmask][0]
        z = data['z'][0]
        fxs[i] = fx
        fys[i] = fy

    fs = np.sqrt(fxs**2+fys**2)

    # calculate forces from scratch
    fljx,fljy,flj,fx_err,fy_err,f_err = f_and_df(xs,ys,muxs,muys,
                                                 fp=fp)
    rs = np.sqrt(xs**2+ys**2)
    
    slabel = 'LAMMPS'  # data output in dump file
    clabel = 'Python'  # calculated in this script

    # plot magnitude of force vs r
    axarr[0].plot(rs[:-1],fs[1:],'o',label=slabel)
    axarr[0].plot(rs,flj,'k-',label=clabel)
    axarr[0].set_ylabel(r'$|\bm{F}_i(r)|$')
    axarr[0].legend(loc='best')

    # plot error between calculated and output force
    axarr[1].plot(rs[:-1],np.abs(flj[:-1]-fs[1:]),'o',label=slabel)
    axarr[1].plot(rs,np.abs(f_err),'k-',label=clabel)
    axarr[1].set_ylabel(r'$|\delta \bm{F}|$')
    axarr[1].set_yscale('log')

    # plot x-component of force vs r
    axarr[2].plot(rs[:-1],fxs[1:],'o')
    axarr[2].plot(rs,fljx,'k-')
    axarr[2].set_ylabel(r'$(\bm{F}_i(r))_x$')

    # plot error for x-component of force
    axarr[3].plot(rs[:-1],np.abs(fljx[:-1]-fxs[1:]),'o',label=slabel)
    axarr[3].plot(rs,np.abs(fx_err),'k-',label=clabel)
    axarr[3].set_ylabel(r'$\delta(\bm{F}_i(r))_x$')
    axarr[3].set_yscale('log')

    # plot y-component of force vs r
    axarr[4].plot(rs[:-1],fys[1:],'o')
    axarr[4].plot(rs,fljy,'k-')
    axarr[4].set_ylabel(r'$(\bm{F}_i(r))_y$')

    # plot error for y-component of force
    axarr[5].plot(rs[:-1],np.abs(fljy[:-1]-fys[1:]),'o',label=slabel)
    axarr[5].plot(rs,np.abs(fy_err),'k-',label=clabel)
    axarr[5].set_ylabel(r'$\delta(\bm{F}_i(r))_y$')
    axarr[5].set_yscale('log')
    axarr[5].set_xlabel(r'$r$')

    
    fig.subplots_adjust(left=0.25,right=0.95,top=0.95,hspace=0.05)
    fig.savefig(f'results/f_vs_r_{fp}.pdf')
    plt.show()

