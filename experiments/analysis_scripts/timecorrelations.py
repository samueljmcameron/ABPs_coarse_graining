import numpy as np
from scipy.spatial.distance import pdist, squareform

from scipy.special import jv, yv, i0e,i1e

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
import sys, os
import pickle

from dumpfile import DumpFile


################################################################################

def save_obj(obj, name):
    """
    This function pickles an object and saves it to the disk.
    """
    with open(name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        print ('saved data to \'%s\'.' % (name + '.pkl'))
        
        
def load_obj(name):
    """
    This function loads a pickled object.
    """
    with open(name + '.pkl', 'rb') as f:
        print ('loading from \'%s\' ...' % (name + '.pkl'))
        return pickle.load(f)

    
def make_sure_path_exists(path):
    import errno
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def calculate_items(snaps, min_neigh=4, cutoff=1.5, MAXnb=100,
                    nbins=2000, nbinsq=50, Pe=10, rho_0=0.60,
                    CG_flag = False):

    """ 
    snaps is an n-dimensional list which holds all the data from a 
    dump file. Each item in the list is a snapshot dict of the per
    atom properties measured in the simulation.
    """
    import fortran_tools as ft


    outer_vec = []
    inner_vec = []
    corr, corr_b, corr_in, count, lindemann = {}, {}, {}, {}, {}
    MSD, Q, Q2, g6t, g6t_re, g6t_im = {}, {}, {}, {}, {}, {}
    for t1,snap1 in enumerate(snaps):

        # for each snapshot in the dump file data

        box=snap1['box']
        ref_coords = snap1['ucoords']
        mus  = snap1['mus']
        if CG_flag:
            # if time averaged velocities are computed
            # at each time step
            # (to reduce noise of velocity per atom
            # at each timestep)
            vs  = snap1['CG_vs']
        else:
            vs = np.column_stack((snapshot['vx'],
                                  snapshot['vy']))
        

        tmp_list = ft.distance_matrix(ref_coords,
                                      snap1['c_psi6[1]'],
                                      snap1['c_psi6[2]'],
                                      snap1['mus'],
                                      snap1['local_density'],
                                      box, cutoff, 1,rho_0,
                                      MAXnb, nbins,nbinsq,
                                      len(ref_coords),2)

        # distance matrix between particle pairs
        ref_distance, ref_dis_x, ref_dis_y = tmp_list[:3]
        # number of neighbours for all particles
        ref_num_nb, ref_list_nb = tmp_list[3:5]
        
        # correlation functions and structure functions
        g, g6, g6re, g6im, sq = tmp_list[5:10]
        g_ori, g_dp, g_dp_tr, g_pr, s_pr = tmp_list[10:]

        # load array which each atom is labelled as being a member
        # of a cluster, with those atoms having label 0 being
        # members of the largest cluster.
        # Size of this array is the number of atoms in the simulation
        cl_i = snap1['new_c_c1']

        # load array of cluster sizes, going from largest to
        # smallest
        # Size of this array is N_clusters
        cl_s = np.array(sorted(snap1['cluster_sizes'].values(),
                               reverse=True))
        Nc = snap1['N_clusters']
        min_cluster_size = 2
        # compute angular momentum and linear momentum of clusters
        cluster_AngM, cluster_LinM = ft.get_cluster_omega(ref_coords,
                                                          vs, box,
                                                          cl_i, cl_s,
                                                          len(ref_coords),
                                                          2, Nc)
        # compute mean squared angular momentum of clusters
        RMS_AngMom2 = np.nanmean(np.where(cl_s>=min_cluster_size,
                                          np.multiply(cluster_AngM,
                                                      cluster_AngM),
                                          np.nan))
        RMS_AngMom = np.sqrt(RMS_AngMom2)
        
        # compute mean squared linear momentum of clusters
        RMS_LinMom2 = np.nanmean(np.where(cl_s>=min_cluster_size,
                                          np.multiply(cluster_LinM,
                                                      cluster_LinM),
                                          np.nan))
        RMS_LinMom = np.sqrt(RMS_LinMom2)

        # compute mean cluster size
        cluster_size = np.nanmean(np.where(cl_s>=min_cluster_size,
                                           cl_s,np.nan))
        # compute (normalized) mean polarisation
        polarisation = np.linalg.norm(np.mean(mus,axis=0))

        print(f"cluster info: RMS_L, RMS_M, size, "
              f"no_clusters(>={min_cluster_size}) "
              "no_clusters")
        print(RMS_AngMom,RMS_LinMom,cluster_size,
              np.count_nonzero( cl_s >= min_cluster_size ),
              len(cl_s))

        if (Pe>0) : 
            mask = np.zeros_like(ref_distance, dtype=np.bool)

            mask[np.tril_indices_from(mask, k=-1)] = True
            rij = ref_distance[mask]
            rij_x = ref_dis_x[mask]
            rij_y = ref_dis_y[mask]

            Wij = W(Pe, rij)
            Wij_rij_x = 2*np.sum(Wij * rij_x)
            Wij_rij_y = 2*np.sum(Wij * rij_y)
            Wij_rij_norm = np.sqrt(Wij_rij_x*Wij_rij_x
                                   + Wij_rij_y*Wij_rij_y)
            fac = i1e(Wij_rij_norm)/i0e(Wij_rij_norm)/Wij_rij_norm

            Wr2 = 2.0 * np.sum( Wij * rij*rij )
            
            print(f"factor={fac}  Wij_rij_x={Wij_rij_x} "
                  f"Wij_rij_y={Wij_rij_y} Wij_rij_norm={Wij_rij_norm}")
            
            Wr2_full = fac * Wr2 / mask.shape[0] / (mask.shape[1]-1) / 2
            Wr2 = Wr2 / 2 / mask.shape[0] / (mask.shape[1]-1) / 2
        else:
            Wr2 = 0
            Wr2_full = 0


#ITEM: ATOMS id type x y xu yu mux muy fx fy tqz v_psi6
# c_psi6[1] c_psi6[2] f_cg[1] f_cg[2] f_cg[3] f_cg[4] f_cg[5]
# f_cg[6] f_cg[7] c_c1 '<psi6_re>', '<psi6_re^2>', '<psi6_im>',
# '<psi6_im^2>', '<mux>', '<mux^2>', '<muy>', '<muy^2>'


        if t1==0 :
            # beginning of time averages
            p6re =  np.mean(snap1['c_psi6[1]'])
            p6im =  np.mean(snap1['c_psi6[2]'])
            p6   =  np.absolute(np.complex(p6re, p6im))
            sum_psi6  = p6
            sum_psi62 = p6*p6
            sum_psi6_cmplx = np.complex(p6re, p6im)
            sum_mux  = np.mean(snap1['mux'])
            sum_mux2 = np.mean(np.array(snap1['mux']) ** 2)
            sum_muy  = np.mean(snap1['muy'])
            sum_muy2 = np.mean(np.array(snap1['muy']) ** 2)
            theta = np.arctan2(snap1['muy'], snap1['mux'])
            sum_theta  = np.mean(theta)
            sum_theta2 = np.mean(theta ** 2)
            nematic    = (2.*np.cos(theta)**2  - 1.)
            sum_nematic  = np.mean( nematic )
            sum_nematic2  = np.mean( nematic**2 )
            
            sum_RMS_AngMom = RMS_AngMom 
            sum_RMS_AngMom2 = RMS_AngMom * RMS_AngMom
            sum_RMS_LinMom = RMS_LinMom
            sum_RMS_LinMom2 = RMS_LinMom * RMS_LinMom
            sum_cluster_size = cluster_size
            sum_polarisation = polarisation
            
            
            sum_g = np.matrix(g)
            sum_g6 = np.matrix(g6)
            
            sum_g6re = np.matrix(g6re)
            sum_g6im = np.matrix(g6im)
            sum_sq = np.array(sq)
            sum_g_ori = np.array(g_ori)
            sum_g_dp = np.array(g_dp)
            sum_g_dp_tr = np.array(g_dp_tr)
            g_pr = np.array(g_pr)
            sum_g_pr = np.array(g_pr)
            sum_Wr2 = Wr2
            sum_Wr2_full = Wr2_full
            sum_pij_rij = s_pr
            
            
            g_cnt = 1
        
        else:
            # add to time averages
            p6re =  np.mean(snap1['c_psi6[1]'])
            p6im =  np.mean(snap1['c_psi6[2]'])
            p6   =  np.absolute(np.complex(p6re, p6im))
            sum_psi6  += p6
            sum_psi62 += p6*p6
            sum_psi6_cmplx += np.complex(p6re, p6im)
            sum_mux  += np.mean(snap1['mux'])
            sum_mux2 += np.mean(np.array(snap1['mux']) ** 2)
            sum_muy  += np.mean(snap1['muy'])
            sum_muy2 += np.mean(np.array(snap1['muy']) ** 2)
            theta = np.arctan2(snap1['muy'], snap1['mux'])
            sum_theta  += np.mean(theta)
            sum_theta2 += np.mean(theta ** 2)
            nematic    = (2.*np.cos(theta)**2  - 1.)
            sum_nematic  += np.mean( nematic )
            sum_nematic2  += np.mean( nematic**2 )
            
            sum_RMS_AngMom += RMS_AngMom
            sum_RMS_AngMom2 += RMS_AngMom * RMS_AngMom
            sum_RMS_LinMom += RMS_LinMom
            sum_RMS_LinMom2 += RMS_LinMom * RMS_LinMom
            sum_cluster_size += cluster_size
            sum_polarisation += polarisation
            
            sum_g += np.matrix(g)
            sum_g6   += np.matrix(g6)
            sum_g6re += np.matrix(g6re)
            sum_g6im += np.matrix(g6im)
            sum_sq += np.array(sq)
            
            
            sum_g_ori += np.array(g_ori)
            sum_g_dp += np.array(g_dp)
            sum_g_dp_tr += np.array(g_dp_tr)
            sum_g_pr += np.array(g_pr)
            sum_Wr2 += Wr2
            sum_Wr2_full += Wr2_full
            sum_pij_rij += s_pr
            
            g_cnt += 1

        # mask for distance matrix
        nb = np.logical_and(ref_distance<cutoff, ref_distance > 0)

        # count number of neighbours each particle has
        num_neighs_ref = ref_num_nb
        # determine whether particles are surrounded by other particles
        inner_particles_ref = num_neighs_ref > min_neigh
        # or not
        outer_particles_ref = np.logical_not(inner_particles_ref)

        # double count number of pairs which have more than four neighs
        norm_in = np.sum(np.multiply(nb[inner_particles_ref]
                                     [:,inner_particles_ref],
                                     nb[inner_particles_ref]
                                     [:,inner_particles_ref]))+0.0
        # double count number of pairs which have less than four neighs
        norm_b =  np.sum(np.multiply(nb[outer_particles_ref]
                                     [:,outer_particles_ref],
                                     nb[outer_particles_ref]
                                     [:,outer_particles_ref]))+0.0

        # double count number of pairs which have neighbours
        norm_all = np.sum(np.multiply(nb,nb))+0.0

        outer_n = outer_particles_ref.sum()
        inner_n = inner_particles_ref.sum()
        print 'boundary/inner = %s/%s' %(outer_n,inner_n )
        inner, outer = get_inner_outer(snap1)
        inner_vec.extend(inner)        
        outer_vec.extend(outer)
        
        for t2 in range(t1+1, ts):
            # for a snapshot at a later time
            snap2 = snaps[t2]
            t = snap2['step'] - snap1['step']
            coords = snap2['ucoords']#[:30]


            tmp_list = ft.distance_matrix(coords,
                                          snap2['c_psi6[1]'],
                                          snap2['c_psi6[2]'],
                                          snap2['mus'],
                                          snap2['local_density'],
                                          box, cutoff, 0,rho_0,
                                          MAXnb, nbins,nbinsq,
                                          len(coords),2)

            # distance matrix between particle pairs
            distance, distance_x, ref_distance_y = tmp_list[:3]
            # number of neighbours for all particles
            num_nb, list_nb = tmp_list[3:5]

            # correlation functions and structure functions,
            # which are all not calculated and so are 0 in
            # what follows
            # g, g6, g6re, g6im, sq = tmp_list[5:10]
            # g_ori, g_dp, g_dp_tr, g_pr, s_pr = tmp_list[10:]

            # mask for distance matrix
            nb1  = np.logical_and(distance<cutoff, distance > 0)


            out = '%s '%t
            count[t] = (count.get(t, 0)) + 1


            # compute number of neighbours which are still neighbours
            # at a later time
            c = np.sum(np.multiply(nb,nb1))
            c /= norm_all
            corr[t] = (corr.get(t, 0)) + c
            out += '%s '%c

            # compute number of neighbours with more than min_neigh
            # neighbours that are still neighbours at a later time
            c = np.sum(np.multiply(nb[inner_particles_ref]
                                   [:,inner_particles_ref],
                                   nb1[inner_particles_ref]
                                   [:,inner_particles_ref]))
            c /= norm_in
            corr_in[t] = (corr_in.get(t, 0)) + c
            out += '%s '%c

            # compute number of neighbours with less than min_neigh
            # neighbours that are still neighbours at a later time
            c = np.sum(np.multiply(nb[outer_particles_ref]
                                   [:,outer_particles_ref],
                                   nb1[outer_particles_ref]
                                   [:,outer_particles_ref]))
            c /= norm_b
            corr_b[t] = (corr_b.get(t, 0)) + c
            out += '%s '%c

            
            d = coords-ref_coords 
            Dsq = [ ((u-d[j-1])**2).sum() for i, u in enumerate(d)
                    for j in ref_list_nb[i,0:ref_num_nb[i]]   ]
            if len(Dsq)> 0 : 
                lindemann[t] = (lindemann.get(t, 0)) + np.mean(Dsq)


            # MSD of clusters
            cl_i = snap1['new_c_c1']
            cl_s = sorted(snap1['cluster_sizes'].values(), reverse=True)
            Nc = snap1['N_clusters'] 
            com_MSD = ft.get_cluster_msd(ref_coords, coords, box,
                                         cl_i, cl_s, len(coords), 2, Nc)
            MSD[t] = (MSD.get(t, 0)) + com_MSD

            # compute overlap autocorrelation function and psi6
            # autocorrelation function
            tmplist = ft.get_overlap(ref_coords, coords, snap1['c_psi6[1]'],
                                     snap1['c_psi6[2]'], snap2['c_psi6[1]'],
                                     snap2['c_psi6[2]'], box, len(coords), 2)

            overlap, g6t_abs, g6_re, g6_im = tmplist

            Q[t] = (Q.get(t, 0)) + overlap
            Q2[t] = (Q2.get(t, 0)) + overlap*overlap
            g6t_re[t] = (g6t_re.get(t, 0)) + g6_re
            g6t_im[t] = (g6t_im.get(t, 0)) + g6_im
            g6t   [t] = (g6t.get(t, 0)) + g6t_abs

            print(f"{out} , {lindemann[t]} , {MSD[t]} , {Q[t]} , "
                  f"{g6t_re[t]} , {g6t_im[t]}")
  
    ret_t=[corr, corr_b, corr_in, count, lindemann,
           MSD, Q, Q2, g6t, g6t_re, g6t_im]
    
    ret_o=[g_cnt,
           sum_psi6, sum_psi62, sum_psi6_cmplx,
           sum_mux, sum_mux2,
           sum_muy, sum_muy2,
           sum_theta, sum_theta2,
           sum_nematic, sum_nematic2,
           sum_g, sum_g6, sum_g6re, sum_g6im, sum_sq,
           sum_g_ori, sum_g_dp , sum_g_dp_tr, sum_g_pr,
           sum_Wr2, sum_Wr2_full, sum_pij_rij,
           sum_RMS_LinMom, sum_RMS_LinMom2, sum_RMS_AngMom,
           sum_RMS_AngMom2, sum_cluster_size, sum_polarisation,
           inner_vec, outer_vec]
    return ret_t, ret_o




if len(sys.argv) == 4 :
   attraction = sys.argv[1]
   phi = sys.argv[2]
   fp = sys.argv[3]
else:
    print 'Usage: %s ttraction phi fp' % (sys.argv[0])
    sys.exit(1)



module_name = 'fortran_tools'
if not os.path.isfile('fortran_tools'+'.so'): 
    status, output, command = compile_fortran(fortran_source, module_name, extra_flags='-fopenmp -lgomp')


trajs=[]
trajs_cnt = 0
avg_ns={}
items_t = ['corr', 'corr_b', 'corr_in', 'count', 'lindemann', 'MSD', 'Q', 'Q2', 'g6t', 'g6t_re', 'g6t_im'] #, '<psi6_re>', '<psi6_re^2>', '<psi6_im>', '<psi6_im^2>', '<mux>', '<mux^2>', '<muy>', '<muy^2>', 'g', 'g6re', 'g6im', 'sq']
items_o = ['count_o','<psi6>', '<psi6^2>', '<psi6_cmplx>', '<mux>', '<mux^2>', '<muy>', '<muy^2>', '<theta>', '<theta^2>', '<nematic_s>', '<nematic_s^2>', 'g', 'g6', 'g6re', 'g6im', '2d_sq', 
           'g_ori', 'g_dp', 'g_dp_tr', 'g_pr', 'Wr2', 'Wr2_full', 'pij_rij', 
           'RMS_M', 'RMS_M^2', 'RMS_L','RMS_L^2','cluster_size', 'polarisation']

sum_dic = {}
for item in items_t + items_o :
  sum_dic[item] = {}
inner, outer  = [], []

paths = ['../RUN_ABP_%s'%i for i in range(1,4)]

for i_traj, path in enumerate(paths) :  
  try:
      #if fp=='0':
      #     filename = '%s/Quincke_0.0_0.0_0.0_10.0/trajs/dump_%s_0.0_%s.lammpstrj'%(path, fp, phi)
      #else:
      filename = '%s/Quincke_0_0_0_%s/trajs/dump_%s_0.0_%s.lammpstrj'%(path, attraction, fp, phi)        

      print (filename)
      import re
      RUN_no=int(re.findall(r'RUN_ABP_[\d]',filename)[0].replace('RUN_ABP_',''))
      traj = read_dump(filename, min_step=0) #, max_step=3000000-1)
      snap_n = len(traj)
      if snap_n < 40 : 
           traj = traj[0::8]
      elif snap_n < 100 :
           traj = traj[0::14]
      elif snap_n < 200 :
           traj = traj[0::30]
      else:
           traj = traj[::6]


      #if RUN_no < 7  and RUN_no > 1 and not RUN_no==3: 
      #   traj = [ snap for snap in traj if snap['step']>1000000]
      print '%s (%s) snapshots are loaded from %s\n'%(len(traj), snap_n, filename)
      if (len(traj)<2): continue
      returned_items_t, returned_items_o = calculate_items(traj, min_neigh=4, cutoff=1.5, MAXnb=50, nbins=8000, nbinsq=100, Pe=np.float(fp), rho_0 = np.float(phi))
      ns = average_cluster_size_dist(traj)
  except:
      print sys.exc_info()[0]
      continue

  trajs_cnt += 1
  try:
    for i, item in enumerate(items_o):
        sum_dic[item] += returned_items_o[i]
  except:
    for i, item in enumerate(items_o):
        sum_dic[item] = returned_items_o[i]

  for i, item in enumerate(items_t): 
      for key in sorted(returned_items_t[0]):
          sum_dic[item][key] = sum_dic[item].get(key,0) + returned_items_t[i][key]

  for key in ns :               
      avg_ns[key] = avg_ns.get(key, 0) + ns[key] 

  inner.extend(returned_items_o[-2])
  outer.extend(returned_items_o[-1])

  save_obj(sum_dic, './tmp/nb_corr_%s_phi%s_fp%s'%(attraction, phi, fp))
  #dummy = load_obj('./tmp/nb_corr_%s_phi%s_fp%s'%(attraction, phi, fp))

if len(sum_dic['corr']) == 0: 
    print "error: failed to analyse trajectorie(s)."
    sys.exit()

for key in avg_ns:
  avg_ns[key] /= trajs_cnt

for item in items_o:
  if item=='count_o' : continue	
  #sum_dic[item] /= trajs_cnt # now count_o has the number of samples
  print '%s total snapshots analyses in %s trajectories.'% (sum_dic['count_o'], trajs_cnt)
  sum_dic[item] /= sum_dic['count_o']




def func(x, tau):
    return np.exp(-x/tau)
def func_stretched(x, tau, beta):
    return np.exp(-(x/tau)**(beta))

out=''
for key in sorted(sum_dic['corr']):
    for item in items_t:
        if item=='count' : continue
        sum_dic[item][key] /= sum_dic['count'][key]
    string =  ' '.join(['%s'%sum_dic[item][key] for item in items_t])
    out += '%s %s\n' % (key, string)




string=''
plt.figure(figsize=(8,6))
for item in ['corr', 'corr_b', 'corr_in', 'Q', 'g6t', 'g6t_re', 'g6t_im']:
    x, y = [], []
    for key in sorted(sum_dic['corr']):
        x.append(key)
        y.append(sum_dic[item][key])
    x, y = np.array(x), np.array(y)
    x    = x*0.000001
    try:
       popt, pcov = curve_fit(func_stretched, x, y) #, bounds=([0], [np.inf]))
       perr = np.sqrt(np.diag(pcov))
       string += '%s %s  %s %s   ' % (popt[0], perr[0], popt[1], perr[1])
       plt.plot(x, y, 'o', label=item)
       plt.plot(x, func_stretched(x, *popt), '-', label=r'$\exp(-(t/%.1f)^{%.1f})$' % (popt[0], popt[1]))
    except:
       popt, pcov = None, None
       perr = None

plt.title('attraction=%s, phi=%s, fp=%s'% (attraction, phi, fp))
plt.xlabel(r'time $/\tau$')
plt.ylabel('$C(t)$')
plt.legend()
with PdfPages('./tmp/nb_corr_%s_phi%s_fp%s.pdf'% (attraction, phi, fp)) as pdf:
    pdf.savefig(transparent=True)

out = '#time %s | %s %s %s   %s\n'%(' '.join(items_t), attraction, phi, fp, string ) + out


make_sure_path_exists('./results')

with open('./results/ns_%s_phi%s_fp%s.txt'%(attraction, phi, fp), 'w') as f:
  print ('# cluster_size freq.')
  for key in sorted(avg_ns):
    print  ('%s %s'%  (key, avg_ns[key] / np.sum(avg_ns.values())))
    f.write('%s %s\n'%(key, avg_ns[key] / np.sum(avg_ns.values())))
    
with open('./results/nb_corr_%s_phi%s_fp%s.txt'%(attraction, phi, fp), 'w') as f:
  f.write('%s'%out)
print '%s\n'%out

for item in ['g', 'g6', 'g6re', 'g6im', 'g_ori', 'g_dp', 'g_dp_tr', 'g_pr']:
   with open('./results/%s_%s_phi%s_fp%s.txt'%(item, attraction, phi, fp), 'w') as f:
       for i in range(len(sum_dic[item][:,0])):
           f.write('%s %s %s %s %s\n'%(sum_dic[item][i,0], sum_dic[item][i,1], sum_dic[item][i,2], sum_dic[item][i,3], sum_dic[item][i,4] ))

item='2d_sq'
print (sum_dic[item])
with open('./results/%s_%s_phi%s_fp%s.txt'%(item, attraction, phi, fp), 'w') as f:
    for a1 in sum_dic[item]:
        for a2 in a1:
            f.write('%s %s %s\n'%(a2[0], a2[1], a2[2]))

items = ['count_o','<psi6>', '<psi6^2>', '<psi6_cmplx>', '<mux>', '<mux^2>', '<muy>', '<muy^2>', '<theta>', '<theta^2>', '<nematic_s>', '<nematic_s^2>', 'RMS_M', 'RMS_M^2', 'RMS_L','RMS_L^2', 'cluster_size', 'polarisation', 'Wr2', 'Wr2_full', 'pij_rij']
with open('./results/fluctuations_%s_phi%s_fp%s.txt'%(attraction, phi, fp), 'w') as f:
  ret = ' '.join(items)+' pij_rij_3 pij_rij_4 pij_rij_5\n'
  for item in items:
    if item=='pij_rij' : 
       ret += '%s %s %s %s' % (sum_dic[item][0], sum_dic[item][1],sum_dic[item][2],sum_dic[item][3] )
    else:
       ret += '%s ' % (sum_dic[item])

  f.write(ret)

item='inner-outer'
with open('./results/%s_%s_phi%s_fp%s.txt'%(item, attraction, phi, fp), 'w') as f:
    for i in range(len(outer)):
        f.write('%s %s\n'%(inner[i], outer[i]))


class StructureFunction(object):

    """ WORK IN PROGRESS I THINK?????? """

    def __init__(self):

        return

    def random_unit_vector3D(self):
        """
        Generates a random 3D unit vector (direction) with a uniform spherical
        distribution.

        Algo from http://stackoverflow.com/questions/5408276/python-uniform
        -spherical-distribution:return:
        """

        phi = np.random.uniform(0,np.pi*2)
        costheta = np.random.uniform(-1,1)

        theta = np.arccos( costheta )
        x = np.sin( theta) * np.cos( phi )
        y = np.sin( theta) * np.sin( phi )
        z = np.cos( theta )
        return (x,y,z)

    def random_unit_vector2D(self):
        """
        Generates a random 2D unit vector (direction) with a uniform spherical
        distribution.

        Algo from http://stackoverflow.com/questions/5408276/python-uniform
        -spherical-distribution:return:
        """

        phi = np.random.uniform(0,np.pi*2)

        x =  np.cos( phi )
        y =  np.sin( phi )
        return (x,y)

    
    def get_sq(x, drs, q): 
        sqs = map(lambda r: 2*np.sin(qv*r)/(qv*r), drs)
        return np.sum(sqs) / len(x)

    def get_average_sq(snap, qmod):
        np.random.seed()
        start = time.time()
        x = snap['coords']
        box = snap['box']
        N = snap['N'] 

        drs = get_pair_distances(x, box)
        Ndir = 30
        sq = 0.0
        for idir in range(Ndir):
            q = np.array(random_unit_vector2D()) * qmod
            sq += get_sq(x, drs, q)        
        end = time.time()
        print(f"averaging s(q) for q = {qmod} over {Ndir} directions "
              f"for {N} atoms. took {end-start} (s). "
              "(process {current_process().pid})")
        return sq/Ndir
