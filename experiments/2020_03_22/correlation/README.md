This folder compares the correlation functions computed directly
via the 'rerun' and 'compute rdf' commands in lammps, with the
correlation functions compute using the python wrapper of fortran
code, in ../../abp_analysis/fortrantools.f95 . The former is much
quicker, and both are similar in output, so I will likely use the
built in lammps one from now on.

To compute the radial distribution function (rdf), append the
following lines too the lammps in.abp file:


compute my_rdf all rdf 2000 cutoff 20.0

fix 1 all ave/time 1 1 1 c_my_rdf[*] file g_${fp}_${rho}.rdf mode vector

comm_modify mode single cutoff 22.0

rerun dump_${fp}_${rho}.lammpstrj first 1000000 dump x y