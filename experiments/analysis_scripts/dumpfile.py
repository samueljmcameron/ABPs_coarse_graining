import numpy as np

class DumpFile(object):

    def __init__(self,filer,voronoi_flag=True,cg_flag=True):

        self.filer = filer
        self.clusters_computed = False
        self.voronoi_flag = voronoi_flag
        self.cg_flag=cg_flag

        return
    
    def read_dump(self,min_step=0, max_step=1e10,pstatus=False):
        """
        Read data from dump file self.filer.

        Example: If you have a simulation with k atoms, and the (per-atom)
        info is dumped into filer n times total, then this function
        returns a list with n items in it. Each item is a dictionary,
        where the 'key' is the name of the measured quantity (e.g. 'x'
        would be the x-coordinate), and value is typically an array (for
        key 'x', the array would be k atoms long). Extra keys/values that
        are manipulations of the original dump data are also included in
        each item within the output list.


        """
        
        snapshot = {} # dictionary to store dump data for current timestep
        data=[] # array which will hold all the snapshots
        snapshot['traj'] = self.filer
        minstep_flag = False
        read_flag = False
        if (max_step<min_step) :
            raise Exception('Max_step is smaller than min_step.')

        cnt = 0

        with open(self.filer,'r') as f:
            while True:
                line=f.readline().strip('\n')
                if not line: break
                items = line.split()
                if items[0] == 'ITEM:':
                    try:
                        if items[1] == 'TIMESTEP':
                            step = int(f.readline().split(' ')[0])
                            if step > max_step :
                                self.max_step_print(max_step,data,pstatus)
                                return data
                            if ( step >= min_step and minstep_flag == False) :
                                self.min_step_print(min_step,pstatus)
                                minstep_flag = True
                            if (step >= min_step and step <= max_step):
                                read_flag = True
                                cnt += 1
                                self.curr_step_print(cnt,step,pstatus)
                            snapshot['step'] =  step
                        if items[1] == 'NUMBER':
                            N = int(f.readline().split()[0])
                            snapshot['N'] =  N
                        if items[1] == 'BOX':
                            line = f.readline().split()
                            box_x = float(line[1]) - float(line[0])
                            line = f.readline().split()
                            box_y = float(line[1]) - float(line[0])
                            line = f.readline().split()
                            box_z = float(line[1]) - float(line[0])
                            snapshot['z_size'] = box_z
                            snapshot['box'] = np.array([box_x, box_y])
                        if items[1] == 'ATOMS':

                            header = items[2:]
                            x = np.zeros( (N, len(header)) )
                            #if pstatus:
                            #    print(header)
                            for i in range(N):
                                line = f.readline().strip('\n').split()
                                for j in range(len(header)):
                                    x[ int(line[0])-1, j] = float(line[j])

                            self.set_snap_atom_vals(header,x,snapshot)

                            self.set_snap_vectors(header,snapshot)
                            if (read_flag): data.append(snapshot.copy())
                            snapshot = {}
                            snapshot['traj'] = self.filer

                    except:
                        print(f'Error detected in reading a snapshot'
                              f'in {self.filer}')

        return data

    def max_step_print(self,max_step,data,pstatus):
        """ Print a message to stdout that the final time-step in the
        dump file has been reached.

        Enabled if pstatus is true.
        """
        
        if pstatus:
            print(f'max_step reached {max_step}')
            print(f"From {self.filer}, last TIMESTEP {data[-1]['step']}")
        return

    def min_step_print(self,min_step,pstatus):
        """ Print a message to stdout that the first time-step in the
        dump file has been reached.

        Enabled if pstatus is true.
        """
        
        if pstatus:
            print(f'From {self.filer}, first TIMESTEP reached {min_step}')
        return

    def curr_step_print(self,cnt,step,pstatus):
        """ Print a message to stdout that the current time-step in the
        dump file is being read.

        Enabled if pstatus is true.
        """
        
        if pstatus:
            print(f'{cnt}: reading TIMESTEP {step}')
        return
    
    def clustering(self,snapshot):
        """ Store sorted cluster information into snapshot dictionary.

        Creates five new keys in snapshot dict:

        'cluster_sizes' value is a dict with keys labeling the
        cluster number and values counting number of atoms in the
        cluster.

        'N_clusters' value is an integer (number of unique clusters).

        'lammps2sorted_id' value is a dictionary, keys are labels of
        original cluster label (as outputted by dump), and values are
        labels of sorted ordering of cluster sizes (i.e. going from
        0 to N_clusters -1).

        'sorted2lammps_id' value is similar to 'lammps2sorted_id'
        above, except keys and values are switched.

        'new_c_c1' is the labels of cluster sizes, maintaining
        the same order as 'c_c1', except that the cluster labels
        have been changed such that the labels of atoms from the
        largest cluster in 'c_c1' (which are completely random)
        are represented as 0s in 'new_c_c1', and the labels of
        atoms from the smallest cluster in 'c_c1' will be labelled
        as N_clusters-1.
        e.g. if snapshot['c_c1'] = [1,2,3,4,4,4,5,1,6], then
        snapshot['new_c_c1'] = [1,2,3,0,0,0,5,6]. Note that
        ordering is random for those clusters with the same
        number of atoms, so the above is equivalent to
        snapshot['new_c_c1'] = [1,6,5,0,0,0,2,3], since there
        are four atoms which are of cluster size one.

        """
        
        cluster_sizes = self.count_items(snapshot['c_c1'])
        snapshot['cluster_sizes'] = cluster_sizes

        snapshot['N_clusters'] = len(cluster_sizes)
        
        # compute list of cluster LABELS, going from largest 
        # cluster to smallest cluster. Note this is NOT a list
        # of cluster sizes.
        sorted_key = sorted(cluster_sizes, key=cluster_sizes.get,
                            reverse=True)

        lammps2sorted_id = { key:c_id  for c_id, key
                             in enumerate(sorted_key)   }
        sorted2lammps_id = { c_id:key  for c_id, key
                             in enumerate(sorted_key)   }
        snapshot['lammps2sorted_id'] = lammps2sorted_id
        snapshot['sorted2lammps_id'] = sorted2lammps_id
        snapshot['new_c_c1'] = [lammps2sorted_id[lmp_id]
                                for lmp_id in snapshot['c_c1']]

        return

    def count_items(self,lst):
        """ Count number of times a label occurs in a list.

        Input lst is the list with labels. 

        Returns a dictionary, with keys being all unique items from
        lst and values being the count of how many times those items
        occur in  lst.

        An example form of this list from a dump file would be the
        cluster (c_c1) row. If the dump entry had 4 atoms, and the
        first and third were in the same cluster, lst = [1,2,1,3]
        and this function would return {'1' : 2, '2' : 1, '3' : 1}.

        """
        
        my_count = {}
        for val in lst:
            my_count[val] = (my_count.get(val, 0)) + 1
        return my_count

    def set_snap_atom_vals(self,header,x,snapshot):
        """ Insert values of the dump file into snapshot dictionary.

        Input header is a string list of names of all the
        quantities measured in the dump file (e.g. 'x', 'mux',
        etc). There are n (where n is dependent on the dump file)
        items in header. E.g. if 'x' and 'y' are the only quantities
        measured in the dump file, then n = 2.

        Input x is an n by N (where N is dependent on the dump file
        and n is defined above) array. E.g. if the dump only outputs
        'x' and 'y' positions of 20 atoms, then n = 2 and N = 20.

        Modifies snapshot dict (since dicts are mutable objects),
        saving each measured quantity (e.g. 'key') for the N atoms
        as snapshot[key] = x[:,key_index].
        """
        
        for j, key in enumerate(header):
            if key in 'id type c_c1 c_nb'.split() : 
                snapshot[key] = np.array( [int(p) for p in x[:,j]]  ) 
                if key == 'c_c1' :
                    self.clustering(snapshot)
                    self.clusters_computed = True
            else:
                snapshot[key] = x[:,j]

        return
    
    def set_snap_vectors(self,header,snapshot):

        """ Store vector quantities in snapshot dict.
        """
        
        snapshot['coords'] = np.column_stack((snapshot['x'],
                                              snapshot['y']))
        
        snapshot['ucoords'] = np.column_stack((snapshot['xu'],
                                               snapshot['yu']))
        
        snapshot['mus'] = np.column_stack((snapshot['mux'],
                                           snapshot['muy']))
        if self.cg_flag:
            # if dump is computing time averaged velocities
            # (averaged over a few time steps)
            snapshot['CG_vs'] = np.column_stack((snapshot['f_cg[1]'],
                                                 snapshot['f_cg[2]']))

        if self.voronoi_flag:
            snapshot['local_density'] = 1/(snapshot['c_voronoi[1]']
                                           /snapshot['z_size'])

        return

if __name__ == "__main__":

    fname = 'test/dump_10_0.04.lammpstrj'
    d = DumpFile(fname)

    data = d.read_dump(pstatus=True)

    snap = data[-1]
    
    print(snap['c_c1'].shape)

    print(type(snap['c_c1']))
    print(type(snap['sorted2lammps_id'][0]))
    print(snap['sorted2lammps_id'][0])

    print(snap['box'])
