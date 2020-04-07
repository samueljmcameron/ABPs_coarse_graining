import numpy as np
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt

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

class ClusterAnalysis(DumpFile):

    def __init__(self,filer,min_step=0,max_step=1e10,pstatus=False):
        
        super().__init__(filer)

        return
    
    def cluster_size_dist(self,snap):
        """ Count number of clusters of each size.

	Returns a dictionary with keys being the size of the
        cluster, and values being the count of each size occuring.

        E.g. if the dump file has c_c1 column with entries
        [1,1,1,2,3,3,5,5], then

        snap['cluster_sizes'] = {1 :3, 2: 1, 3: 2, 5: 2 }.

        So list(snap['cluster_sizes'].values()) = [3,1,2,2]. The
        output of this function would then be
        {3: 1, 1: 1, 2: 2}, meaning there is one cluster of size 3,
        one cluster of size 1, and 2 clusters of size 2.
        """
        
        cluster_sizes = snap['cluster_sizes']
        cluster_size_hist = self.count_items(list(cluster_sizes.values()))
        return cluster_size_hist

    def average_cluster_size_dist(self,data):
        """ Compute average cluster size, where the average is taken
        over all of the time steps in the dump file.

        Returns a dictionary with keys being the sizes of the clusters,
        and values being the average number of occurence of this
        sized cluster.

        """
        
        avg_hist = {}
        cnt = 0
        for snap in data:
            cnt += 1
            cluster_size_hist = self.cluster_size_dist(snap)
            for key in cluster_size_hist :               
                avg_hist[key] = (avg_hist.get(key, 0)
                                 + cluster_size_hist[key] )
        if cnt>0 :
            for key in sorted(avg_hist):
                avg_hist[key] /= float(cnt)
        return avg_hist

    def PBC(self,d, box):
        """ Apply periodic boundary conditions to the distance between
        two points in 2D space.

        d is the distance in real space (i.e. d can be larger than the
        dimensions of box.

        box is the dimensions of the box where periodic boundary
        conditions are being applied. E.g. if the box was 100 by
        100, box = [100,100].
        """
        
        d[0] = d[0] - np.rint(d[0]/box[0])*box[0]
        d[1] = d[1] - np.rint(d[1]/box[1])*box[1]
        return d

    def get_boundary(self,coords, box, min_neigh=4, cutoff=1.5,
                     dist2bound_flag=False):
        """
        Compute an array mask to determine which particles have more
        than min_neigh neighbours.

        Two particles are neighbours if they are within the cutoff
        distance.

        coords are the x,y coordinates of the atoms being
        considered.
        box is the size of the 'box' (i.e. where periodic boundary
        conditions are applied).

        dist2bound_flag=True calculates the distance each
        atom is from a hole in the cluster, e.g. to another atom
        which does not have min_neigh atoms. If the atom itself has
        less than min_neigh atoms, then this distance is 0.
        """

        # compute distances considering periodic boundary conditions
        distance=squareform(pdist(coords,
                                  lambda u, v:
                                  np.sqrt(((self.PBC((u-v), box))**2).sum()
                                  )))

        distance_dummy = distance.copy()

        # fill diagonals of pair distance matrix with value larger than
        # the cutoff.
        # the '4' is arbitrary, it just needs to be a number > 1.
        np.fill_diagonal(distance, 4*cutoff)

        # mask of all neighbours (within cutoff distance of each other)
        nb = distance < cutoff
        # list of all neighbours (within cutoff distance of each other)
        num_neighs = nb.sum(axis=1)

        # for each atom, true if it has more than min_neigh neighbors
        inner_particles = num_neighs > min_neigh

        # for each atom, true if it has less than min_neigh neighbors
        outer_particles = np.logical_not(inner_particles)

        if dist2bound_flag:

            # compute distance from each atom to the nearest outside atom
            dist2bound = np.array([np.min(distance_dummy[i,outer_particles])
                                   for i in range(len(coords))])
            return outer_particles, inner_particles, nb, dist2bound
        else:
            return outer_particles, inner_particles, nb

        return

    def get_clust_mask(self,i,snap,size=False):
        """
        Compute array mask to choose only those coordinates in
        cluster i. Optionally, also compute the size of cluster
        i.
        """

        # snap['sorted2lammps_id'] is a dictionary of clusters,
        # where key goes from 0 to snap['N_clusters']-1, and
        # value goes from largest cluster label to smallest
        # cluster label 
        mask = (snap['c_c1'] == snap['sorted2lammps_id'][i])
        if size:
            clust_size = snap['cluster_sizes']\
                         [snap['sorted2lammps_id'][i]]
            return mask, clust_size
        else:
            return mask
        
        return
        
    def get_inner_outer(self,snap):
        """
        Count number of particles in dense and less dense regions for 
        each cluster.

        Returns two lists, each of length N, where N is the number of
        clusters in the system. For the first returned list, the i^th
        item (i = 0,...,N-1) is a count of the number of particles which
        are in the dense region of the i^th cluster. The second returned
        list is similar, with each item being a count of the number of
        particles in the less dense region of the cluster.
        """


        # inner will store number of particles in dense region of each
        # cluster
        # outer will do the same, but for particles in less dense region
        inner, outer = [],[]
        #load all x y coords at current time step
        coords=snap['coords']

        for i in range(snap['N_clusters']):

            mask,clust_size = self.get_clust_mask(i,snap,size=True)
            
            # ignore clusters which are smaller than 10 particles
            if clust_size <= 10 : continue
            # load in coordinates of each atom in cluster
            clust_coords = coords[mask,:]
            # compute masks for which particles have less and
            # more than min_neigh = 4 neighbours
            mask_outer, mask_inner,nb = self.get_boundary(clust_coords,
                                                          snap['box'])

            inner.append(mask_inner.sum())
            outer.append(mask_outer.sum())
            
        return inner, outer

    def get_psi6_profile(self,snap, nbins, min_r=0, max_r=10,
                         bin_size=0.1):
        """
        Compute cumulative sum of psi6 order parameter, as a
        function of binned distance.

        Returns two arrays. First is an array of length
        nbins. The first component of this array ([0])
        is the cumulative sum of psi6 for all atoms on the
        boundary (or within min_r distance of the boundary).
        The last component of the array ([nbins-1]) is the
        cumulative sum of psi6 for all atoms which are
        farther than bin_size*nbins distance away from
        a boundary.

        Second returned array is just the count of
        atoms in each bin.
        """

        profile = np.zeros(nbins)
        count   = np.zeros(nbins)

        coords=snap['coords']

        for i in range(snap['N_clusters']):

            mask,clust_size = self.get_clust_mask(i,snap,size=True)

            if clust_size <= 10 : 
                continue
            # load in magnitude of psi6 order parameter for each
            # atom in the cluster
            clust_psi6 = snap['v_psi6'][mask]
            # load in coordinates of each atom in cluster
            clust_coords = coords[mask,:]
            # compute masks for which particles have less and
            # more than min_neigh = 4 neighbours
            # also compute distance to the outside of a cluster
            # for all particles
            outer, inner, nb,dist2bound = self.get_boundary(clust_coords,
                                                            snap['box'],
                                                            dist2bound_flag=
                                                            True)

            # for each atom in cluster, compute an integer
            # which is zero if they are on the boundary,
            # and largest in the middle of the cluster
            js = np.floor((dist2bound - min_r)/bin_size) + 1
            js[js>nbins-1] = nbins-1 
            js[js<0] = 0
            js = js.astype(int)
            for k, j_bin in enumerate(js):
                profile[j_bin] += clust_psi6[k]
                count[j_bin]   += 1

        return profile, count


    def get_average_psi6_profile(self,data, min_r=0, max_r=10,
                                 bin_size=0.1):
        """
        Compute psi6 order parameter, averaged over all clusters,
        as function of binned distance.

        Return average psi6 order parameter array and array with
        count of total number of atoms in each binned distance.
        """

        r = np.arange(min_r - bin_size/2.,
                      max_r + bin_size/2.,bin_size)
        nbins = len(r)
        sum_profile = np.zeros(nbins)
        sum_count   = np.zeros(nbins)
        for snap in data:
            profile, count = self.get_psi6_profile(snap, nbins,
                                                   min_r=min_r,
                                                   max_r=max_r,
                                                   bin_size=bin_size)
            sum_profile += profile
            sum_count += count
        for i, cnt in enumerate(sum_count):
            if cnt>0 : sum_profile[i] /= cnt

        return r, sum_profile, sum_count



    def get_lindemann(traj, min_neigh=4, cutoff=1.5):
        """
        Compute

        q(t) = sum_{t2-t1=t}<(r_j(t1,t2)-r_i(t1,t2))**2>_{neighbours},

        where

        r_i(t1,t2) = r_i(t2)-r_i(t1)

        is the difference in position of atom i at time t2
        compared to time t1, and the angle brackets indicate an
        average over neighbours.
        """
        
        lindemann, count={}, {}
        ts = len(traj)
        for t1 in range(ts):
            snap1 = traj[t1]
            box=snap1['box']
            ref_coords = snap1['ucoords']

            # compute array mask to deterimine which particles have
            # more than min_neigh neighbours

            mask_outer, mask_inner,nb = self.get_boundary(ref_coords,
                                                          box,
                                                          min_neigh=
                                                          min_neigh,
                                                          cutoff=
                                                          cutoff)

            for t2 in range(t1+1, ts):
                snap2 = traj[t2]
                t = snap2['step'] - snap1['step']
                coords = snap2['ucoords']
                # compute distances between atoms at the two times,
                # r_i(t2)-r_i(t1)
                d = np.array(map (lambda x: self.PBC(x,box),
                                  coords-ref_coords))
                # compute (r_j(t2)-r_j(t1)-r_i(t2)+r_i(t1))**2 for
                # for nearest neighbours
                Dsq = [ ((PBC((u-d[j]), box))**2).sum() for i, u
                        in enumerate(d) for j, dummy
                        in enumerate(nb[i]) if dummy  ]

                if len(Dsq)> 0 : 
                  lindemann[t] = (lindemann.get(t, 0))+ np.mean(Dsq)
                  count[t] = (count.get(t, 0)) + 1
        return lindemann, count


    def get_average_lindemann(trajs, min_neigh=4, cutoff=1.5):
        """
        Compute lindemann average over multiple dump files.
        """
        sum_lindemann, sum_count = {}, {}
        out = '#time lindemann count\n'
        for traj in trajs:
          lindemann, count = get_lindemann(traj, min_neigh=4, cutoff=1.5)

          for key in lindemann:
            sum_lindemann[key] = sum_lindemann.get(key, 0) + lindemann[key]
            sum_count[key] = sum_count.get(key, 0) + count[key]

        for key in sorted(sum_lindemann):
            sum_lindemann[key] /= sum_count[key]
            out += f"{key} {sum_lindemann[key]} {sum_count[key]}" + "\n"
        return out


    

if __name__ == "__main__":

    fname = 'test/dump_80_0.6.lammpstrj'

    cs = ClusterAnalysis(fname)

    data = cs.read_dump(pstatus=True)
    
    r,sum_profile,sum_count = cs.get_average_psi6_profile(data)

    fig,axarr = plt.subplots(2,sharex=True)

    axarr[0].plot(r,sum_profile)
    axarr[1].plot(r,sum_count)

    plt.show()
