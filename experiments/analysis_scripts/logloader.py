# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.

# log tool

oneline = "Read LAMMPS log files and extract thermodynamic data"

docstr = """
l = log("file1")                     read in one or more log files

  if the same time step occurs multiple times, it is assumed that
  all prior data was equilibration, and it is deleted

nvec = l.nvec                        # of vectors of thermo info
nlen = l.nlen                        length of each vectors
names = l.names                      list of vector names
t,pe,... = l.get("Time","KE",...)    return one or more vectors of values
l.write("file.txt")                  write all vectors to a file
l.write("file.txt","Time","PE",...)  write listed vectors to a file

  get and write allow abbreviated (uniquely) vector names
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   nvec = # of vectors
#   nlen = length of each vector
#   names = list of vector names
#   ptr = dictionary, key = name, value = index into data for which column
#   data[i][j] = 2d array of floats, i = 0 to # of entries, j = 0 to nvecs-1
#   firststr = string that begins a thermo section in log file

#   eof = ptr into incremental file for where to start next read

# Imports and external programs

import sys, re, glob
from os import popen

import numpy as np


# Class definition

class LogLoader(object):

    # --------------------------------------------------------------------

    def __init__(self,fname,remove_chunk = 0,chunk_start = 'Step',
                 chunk_end = 'Loop',merge_data = True):


        # fname is name of log file
        self.fname = fname
        self.chunk_start = chunk_start
        self.chunk_end = chunk_end
        
        self.data = self.__load_data()

        if remove_chunk != None:

            self.data = self.remove_chunk(remove_chunk)

        if merge_data:

            self.data = self.merge_data()

        return
    # --------------------------------------------------------------------

    def remove_chunk(self,i):

        data = []
        
        for j,d in enumerate(self.data):
            if j != i:
                data.append(d)
                
        return data


    # --------------------------------------------------------------------

    def merge_data(self,remove_duplicates=True):
        
        dataset = dict(self.data[0])
        
        for j,ddict in enumerate(self.data[1:]):
            
            for key,value in ddict.items():
                try:
                    val = dataset[key]
                except KeyError:

                    raise Exception('Unable to merge thermo chunks: '
                                    f'{key} variable in thermo '
                                    'chunk {j} is not measured in '
                                    'the previous thermo chunk.')
                dataset[key] = np.concatenate((val,
                                               ddict[key]))

        if remove_duplicates:
            steps = np.copy(dataset['Step'])
            
            dmask = steps[1:]-steps[:-1]
            dmask = np.array(np.concatenate(([True],dmask)),bool)

            for key in dataset:
                dataset[key] = dataset[key][dmask]
            
            assert (dataset['Step'] == np.unique(steps)).all()

            
        return dataset
    
    # --------------------------------------------------------------------
    
    def __load_data(self):

        data = []

        lines = self.__load_lines()
        
        master_masks,header_list = self.__chunk_lines(lines)

        for j,chunk in enumerate(header_list):
            
            d_lines = np.array(lines)[master_masks[j]]
            
            dt = np.array([d.split() for d in d_lines])

            for key,value in chunk.items():

                if key == 'Step':
                    chunk[key] = np.array(dt[:,value],int)
                else:
                    chunk[key] = np.array(dt[:,value],float)
            data.append(chunk)
                                    
        return data

    
    # --------------------------------------------------------------------
    
    def __load_lines(self):
        with open(self.fname) as f:
            
            lines = [line.strip('\n') for line in f.readlines()
                     if line.strip()]

        return lines

    
    # --------------------------------------------------------------------

    def __chunk_lines(self,lines,skipchunks=None):
        """
        Build an array mask which is true for all data lines
        (i.e. lines with thermo output) and false for all
        other lines.

        """
        
        chunk_start = self.chunk_start
        chunk_end = self.chunk_end
        header_list = []
        master_masks = []

        cnt = 0
        while cnt < len(lines):
            
            line = lines[cnt]
            items = line.split()
            
            # if chunk of thermo outputs starts here
            
            if items[0] == chunk_start:

                # construct header containing all thermo variables
                # that will be output, with keys being variable name
                # and value being the column of the variable
                headerdict = self.__build_header(items)

                # construct mask which will select only lines of
                # current thermo output chunk
                mask = np.zeros_like(lines,bool)

                # move to next line, since current line has variable
                # names
                cnt += 1

                line = lines[cnt]
                items = line.split()
                
                # while still reading numbers

                while items[0] != chunk_end:
                    # if reading a slurm file, sometimes MPI
                    # produces messages in the middle of data
                    # chunks about processors, so check to
                    # ensure the line in the data chunk is
                    # actually a data line, which is why the
                    # following try is necessary


                    try:
                        int(items[0])
                    except:
                        print("weird interruption in data stating: "
                              f"{' '.join(items)}")
                    else:
                        mask[cnt] = True

                    cnt += 1
                    line = lines[cnt]
                    items = line.split()

                # add header and masks to lists
                header_list.append(headerdict)
                master_masks.append(mask)

            else:
                cnt+= 1
                
        return master_masks,header_list

    def __build_header(self,items):
        newdict = {}

        for j,item in enumerate(items):
            newdict[item] = j

        return newdict
            

if __name__=='__main__':


    prefix = '../2020_03_31/interactions_pressure/data/'
    
    fname = prefix + 'slurm-756108.out'

    ll = LogLoader(fname,remove_chunk=0)

    #print(ll.data[0])
    #print(ll.data[1])

    ll.merge_data()

    
    #ll.read_header()
