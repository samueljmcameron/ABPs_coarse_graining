# cython: language_level = 3



import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

import sys


cimport header

def s1_array(np.ndarray[DTYPE_t,ndim=1] x, int N):

    cdef int i
    cdef np.ndarray[DTYPE_t,ndim=1] out = np.zeros([N],dtype = DTYPE)

    for i in range(N):
        out[i] = header.s1(x[i])

    return out

def s1_py(x):

    if type(x) is np.ndarray:
        N = len(x)
        return s1_array(x,N)
    
    return header.s1(x)


def s2_array(np.ndarray[DTYPE_t,ndim=1] x, int N):

    cdef int i
    cdef np.ndarray[DTYPE_t,ndim=1] out = np.zeros([N],dtype = DTYPE)

    for i in range(N):
        out[i] = header.s2(x[i])

    return out
    
def s2_py(x):

    if type(x) is np.ndarray:
        N = len(x)
        return s2_array(x,N)

    return header.s2(x)


def q1_array(np.ndarray[DTYPE_t,ndim=1] x, int N):

    cdef int i
    cdef np.ndarray[DTYPE_t,ndim=1] out = np.zeros([N],dtype = DTYPE)

    for i in range(N):
        out[i] = header.q1(x[i])

    return out
    
def q1_py(x):

    if type(x) is np.ndarray:
        N = len(x)
        return q1_array(x,N)

    return header.q1(x)


def q2_array(np.ndarray[DTYPE_t,ndim=1] x, int N):

    cdef int i
    cdef np.ndarray[DTYPE_t,ndim=1] out = np.zeros([N],dtype = DTYPE)

    for i in range(N):
        out[i] = header.q2(x[i])

    return out
    
def q2_py(x):

    if type(x) is np.ndarray:
        N = len(x)
        return q2_array(x,N)

    return header.q2(x)
