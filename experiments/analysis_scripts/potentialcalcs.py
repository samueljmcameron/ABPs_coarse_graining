import numpy as np


from scipy.special import jv, yv, i0e,i1e

import matplotlib.pyplot as plt

import sys, os
import pickle_dump as pd
import fortrantools as ft
from dumpfile import DumpFile

sys.path.append('../scripts')

from softsphere import SoftSphere


def w_stuff(fp,ref_distance,ref_dis_x,ref_dis_y):

    if (fp>0) : 
        mask = np.zeros_like(ref_distance, dtype=np.bool)
        mask[np.tril_indices_from(mask, k=-1)] = True
        
        rij = ref_distance[mask]
        rij_x = ref_dis_x[mask]
        rij_y = ref_dis_y[mask]

        Wij = W(fp, rij)
        Wij_rij_x = 2*np.sum(Wij * rij_x)
        Wij_rij_y = 2*np.sum(Wij * rij_y)
        Wij_rij_norm = np.sqrt(Wij_rij_x*Wij_rij_x + Wij_rij_y*Wij_rij_y)
        fac = i1e(Wij_rij_norm)/i0e(Wij_rij_norm)/Wij_rij_norm
        Wr2 = 2.0 * np.sum( Wij * rij*rij )
        print ('factor=%g  Wij_rij_x=%s Wij_rij_y=%s Wij_rij_norm=%s'%(fac,Wij_rij_x,Wij_rij_y,Wij_rij_norm))
        Wr2_full = fac * Wr2 / mask.shape[0] / (mask.shape[1]-1) / 2
        Wr2 = Wr2 / 2 / mask.shape[0] / (mask.shape[1]-1) / 2
    else:
        Wr2 = 0
        Wr2_full = 0



def spatial_correlations(t1,snap1,ret_o, min_neigh=4, cutoff=1.5,
                         MAXnb=100, nbins=2000, nbinsq=50,
                         Pe=10,rho_0=0.60):

    # for each snapshot in the dump file data

    box=snap1['box']
    ref_coords = snap1['ucoords']
    mus  = snap1['mus']

    # the following returns the list of:

    # distance matrix between particle pairs
    # ref_distance, ref_dis_x, ref_dis_y = tmp_list[:3]
    # number of neighbours for all particles
    #ref_num_nb, ref_list_nb = tmp_list[3:5]

    # correlation functions and structure functions
    #g, g6, g6re, g6im, sq = tmp_list[5:10]
    #g_ori, g_dp, g_dp_tr, g_pr, s_pr = tmp_list[10:]
    
    tmp_list = ft.distance_matrix(ref_coords,
                                  snap1['c_psi6[1]'],
                                  snap1['c_psi6[2]'],
                                  snap1['mus'],
                                  snap1['local_density'],
                                  box, cutoff, 1,rho_0,
                                  MAXnb, nbins,nbinsq,
                                  len(ref_coords),2)

    return tmp_list



def calculate_items(ret_o,snaps, min_neigh=4, cutoff=1.5, MAXnb=100,
                    nbins=2000, nbinsq=50, Pe=10, rho_0=0.60,
                    spatial_correlation_flag = True,
                    cluster_flag = False, CG_flag = True):

    """ 


    snaps is an n-dimensional list which holds all the data from a 
    dump file. Each item in the list is a snapshot dict of the per
    atom properties measured in the simulation.
    """

    ts = len(snaps)
    
    for t1 in range(0,ts):

        snap1 = snaps[t1]
        print(t1,ts)
        # for each snapshot in the dump file data

        box=snap1['box']
        ref_coords = snap1['ucoords']
        mus  = snap1['mus']

        # compute (normalized) mean polarisation
        polarisation = np.linalg.norm(np.mean(mus,axis=0))

        p6re =  np.mean(snap1['c_psi6[1]'])
        p6im =  np.mean(snap1['c_psi6[2]'])
        p6   =  np.absolute(np.complex(p6re, p6im))

        mux = np.mean(snap1['mux'])
        mux2 = np.mean(np.array(snap1['mux'])**2)
        muy = np.mean(snap1['muy'])
        muy2 = np.mean(np.array(snap1['muy'])**2)
        
        theta_Ns = np.arctan2(snap1['muy'], snap1['mux'])
        theta = np.mean(theta_Ns)
        theta2 = np.mean(theta_Ns**2)
        
        nematic_Ns    = (2.*np.cos(theta)**2  - 1.)
        nematic = np.mean(nematic_Ns)
        nematic2 = np.mean(nematic_Ns**2)
        
        # compute time averages
        ret_o['g_cnt'] = ret_o.get('g_cnt',0) + 1
        ret_o['sum_psi6'] = ret_o.get('sum_psi6',0) + p6
        ret_o['sum_psi62'] = ret_o.get('sum_psi62',0) + p6*p6
        ret_o['sum_psi6_re'] = ret_o.get('sum_psi6_re',0) + p6re
        ret_o['sum_psi6_im'] = ret_o.get('sum_psi6_im',0) + p6im
        ret_o['sum_mux'] = ret_o.get('sum_mux',0) + mux
        ret_o['sum_mux2'] = ret_o.get('sum_mux2',0) + mux2
        ret_o['sum_muy'] = ret_o.get('sum_muy',0) + muy
        ret_o['sum_muy2'] = ret_o.get('sum_muy2',0) + muy2

        
        ret_o['sum_theta'] = ret_o.get('sum_theta',0) + theta
        ret_o['sum_theta2'] = ret_o.get('sum_theta2',0) + theta2


        ret_o['sum_nematic'] = ret_o.get('sum_nematic',0) + nematic
        ret_o['sum_nematic2'] = ret_o.get('sum_nematic2',0) + nematic2
        ret_o['polarisation'] = ret_o.get('polarisation',0) + polarisation

        
        if spatial_correlation_flag:
        
            tmp_list = spatial_correlations(t1,snap1, ret_o,min_neigh=4,
                                            cutoff=1.5,MAXnb=100,nbins=2000,
                                            nbinsq=50,Pe=10, rho_0=0.60)

            # distance matrix between particle pairs
            ref_distance, ref_dis_x, ref_dis_y = tmp_list[:3]
            # number of neighbours for all particles
            ref_num_nb, ref_list_nb = tmp_list[3:5]
        
            # correlation functions and structure functions
            g, g6, g6re, g6im, sq = tmp_list[5:10]
            g_ori, g_dp, g_dp_tr, g_pr, s_pr = tmp_list[10:]


            # compute time averages

            g_mat = np.matrix(g)
            g6_mat = np.matrix(g6)
            g6re_mat = np.matrix(g6re)
            g6im_mat = np.matrix(g6im)
            sq_mat = np.array(sq)

            ret_o['sum_g'] = ret_o.get('sum_g',0*g_mat)+g_mat
            ret_o['sum_g6'] = ret_o.get('sum_g6',0*g6_mat)+g6_mat
            ret_o['sum_g6re'] = ret_o.get('sum_g6re',0*g6re_mat)+g6re_mat 
            ret_o['sum_g6im'] = ret_o.get('sum_g6im',0*g6im_mat)+g6im_mat

            ret_o['sum_sq'] = ret_o.get('sum_sq',0*sq_mat)+sq_mat

            g_ori_mat = np.array(g_ori)
            g_dp_mat = np.array(g_dp)
            g_dp_tr_mat = np.array(g_dp_tr)
            g_pr_mat = np.array(g_pr)
            pij_rij_mat = s_pr


            ret_o['sum_g_ori'] = (ret_o.get('sum_g_ori',0*g_ori_mat)
                                  + g_ori_mat)
            ret_o['sum_g_dp'] = (ret_o.get('sum_g_dp',0*g_dp_mat)
                                 + g_dp_mat)
            ret_o['sum_g_dp_tr'] = (ret_o.get('sum_g_dp_tr',0*g_dp_tr_mat)
                                    +g_dp_tr_mat)
            ret_o['sum_g_pr'] = (ret_o.get('sum_g_pr',0*g_pr_mat)
                                 +g_pr_mat)
            ret_o['sum_pij_rij'] = (ret_o.get('sum_pij_rij',0*pij_rij_mat)
                                    + pij_rij_mat)



        
        if cluster_flag:

            tmp_list = cluster_momenta(t1,snap1,
                                       min_cluster_size=min_cluster_size,
                                       CG_flag=CG_flag)

            RMS_AngMom,RMS_AngMom2 = tmp_list[:2]
            RMS_LinMom,RMS_LinMom2,cluster_size = tmp_list[2:]


            # beginning of time averages

            ret_o['sum_RMS_AngMom'] = (ret_o.get('sum_RMS_AngMom',0)
                                       + RMS_AngMom)
            ret_o['sum_RMS_AngMom2'] = (ret_o.get('sum_RMS_AngMom2',0)
                                        + RMS_AngMom2)
            ret_o['sum_RMS_LinMom'] = (ret_o.get('sum_RMS_LinMom',0)
                                       + RMS_LinMom)
            ret_o['sum_RMS_LinMom2'] = (ret_o.get('sum_RMS_LinMom2',0)
                                        + RMS_LinMom2)
            ret_o['sum_cluster_size'] = (ret_o.get('sum_cluster_size',0)
                                         +cluster_size)


    return ret_o


if __name__ == "__main__":

    fname = 'data/ABP/dump_60_0.4.lammpstrj'
    
    dfile = DumpFile(fname)

    data = dfile.read_dump(pstatus=True,min_step=1000000)

    print(len(data))
    
    #data = pd.load_obj(fname)

    ret_o = {}
    calculate_items(ret_o,data, min_neigh=4, cutoff=1.5, MAXnb=100,
                    nbins=2000, nbinsq=50, Pe=10, rho_0=0.60)

    gmatrix = ret_o['sum_g']
    Nsamples = ret_o['g_cnt']

    print(Nsamples)
    rs = gmatrix[:,0]/Nsamples

    gs = gmatrix[:,1]/Nsamples

    plt.plot(rs,gs)
    plt.show()

