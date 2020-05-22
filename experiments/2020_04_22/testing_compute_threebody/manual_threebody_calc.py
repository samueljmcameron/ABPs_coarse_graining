
def g3_reader(inpath,fp,rho,binnum=1000):

    lines = []
    filename = inpath + f'g_{fp}_{rho}.rdf'

    with open(filename,'r') as f:
        count = 0
        while True:
            line = f.readline()
            try:
                if line[0] == "#":
                    continue
            except:
                print('end of file reached!')
                break
            if line.count(' ') == 1:
                a = line.split(' ')
                count += 1
                lines.append(np.array([[float(j) for j in next(f).split(' ')]
                                       for i in range(binnum)],float))


            
    return lines


if __name__=="__main__":

    import numpy as np
    import matplotlib.pyplot as plt
    import os
    import sys

    sys.path.append('../../analysis_scripts')

    from dumpfile import DumpFile
    from pickle_dump import save_obj,load_obj
    import force_calculations as fc
    from fortrantools import threebody

    sys.path.append('../../plotting_scripts')

    from jupyterplots import JupyterPlots

    fig_x,fig_y = JupyterPlots()
    
    prefix = 'data/'

    roundoff = 5e-6

    fp = 0
    rho = 0.3

    fdump = prefix + f'dump_{fp}_{rho}.lammpstrj'

    if os.path.isfile(fdump + '.pkl'):
        dumpdata = load_obj(fdump,pstatus=True)
    else:
        duf = DumpFile(fdump,voronoi_flag=False,cg_flag=False)
        dumpdata = duf.read_dump()
        save_obj(dumpdata,fdump,pstatus=True)

                      
                      
    nbins = 20
    nthetabins = 20

    g3data = g3_reader(prefix,fp,rho,binnum=8000)
    
    cutoff = 3.0
    dimension = 2
    iters = len(dumpdata)
    for i,data in enumerate(dumpdata[1:2]):

        coords = data['coords']
        x,y = coords.T
        box = data['box']



        n = len(x)

        g3 = threebody(coords=coords,box=box,cutoff=cutoff,maxnb=100,
                       n=n,d=dimension,nbins=nbins,nthetabins=nthetabins)

        g3flat = g3.flatten('F')
        print('shape of flattened fortran threebody array: '
              f'{g3flat.shape}')
        print('shape of flattened fortran threebody array, masked '
              'to only include\n non-zero entries: '
              f'{g3flat[g3flat!=0].shape}')

        g3_lammps = g3data[1][:,4]
        print('shape of lammps threebody output: '
              f'{g3_lammps.shape}')

        print('shape of lammps threebody output, masked to only '
              'include non-zero\n entries: '
              f'{g3_lammps[g3_lammps!=0].shape}')
        print('array of all entries in lammps threebody output which '
              'are not equal to\n the output of the flattened fortran '
              'threebody array: '
              f'{g3_lammps[np.logical_not(np.isclose(g3_lammps,g3flat))]}')

        nonzerog3s = g3_lammps #[g3_lammps!=0]
        uniques = np.unique(nonzerog3s,return_counts=True,return_index=True)
        print (uniques[0][uniques[2]%2!=0])
        print (uniques[1][uniques[2]%2!=0])
        print (uniques[2][uniques[2]%2!=0])
        #print(f'{uniques}')
