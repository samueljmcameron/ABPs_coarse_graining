# cython: language_level = 3


from header cimport s1, s2, q1, q2
from cmath import sqrt

cdef void s1_single_point(double x_in,double ans_out) nogil:

    ans_out = s1(x_in)


cdef void s2_single_point(double x_in,double ans_out) nogil:

    ans_out = s2(x_in)

cdef void q1_single_point(double x_in,double ans_out) nogil:


    ans_out = q1(x_in)
    
cdef void q2_single_point(double x_in,double ans_out) nogil:

    ans_out = q2(x_in)


# Boilerplate Cython definitions
#
# Pulls definitions from the Numpy C headers.
#------------------------------------------------------------
from numpy cimport import_array, import_ufunc
from numpy cimport (PyUFunc_FromFuncAndData,
                    PyUFuncGenericFunction)
from numpy cimport NPY_DOUBLE
from numpy cimport PyUFunc_d_d
#Required module initialization
#------------------------------------------------------------
import_array()
import_ufunc()
# ufunc declaration of s1
cdef PyUFuncGenericFunction s1_loop_func[1]
cdef char s1_input_output_types[2]
cdef void *s1_elementwise_funcs[1]
s1_loop_func[0] = PyUFunc_d_d
s1_input_output_types[0] = NPY_DOUBLE
s1_input_output_types[1] = NPY_DOUBLE
s1_elementwise_funcs[0] = <void*>s1_single_point
s1_py = PyUFunc_FromFuncAndData(
    s1_loop_func,
    s1_elementwise_funcs,
    s1_input_output_types,
    1, # number of supported input types
    1, # number of input args
    1, # number of output args
    0, # 'indentity' element, never mind this
    "s1_py", # function name
    "s1_py(x) -> computes int_{0.5}^x I_1(u)/u^12 du", # docstring
    0 # unused
)

# ufunc declaration of s2
cdef PyUFuncGenericFunction s2_loop_func[1]
cdef char s2_input_output_types[2]
cdef void *s2_elementwise_funcs[1]
s2_loop_func[0] = PyUFunc_d_d
s2_input_output_types[0] = NPY_DOUBLE
s2_input_output_types[1] = NPY_DOUBLE
s2_elementwise_funcs[0] = <void*>s2_single_point
s2_py = PyUFunc_FromFuncAndData(
    s2_loop_func,
    s2_elementwise_funcs,
    s2_input_output_types,
    1, # number of supported input types
    1, # number of input args
    1, # number of output args
    0, # 'indentity' element, never mind this
    "s2_py", # function name
    "s2_py(x) -> computes int_{0.5}^x I_1(u)/u^6 du", # docstring
    0 # unused
)

# ufunc declaration of q1
cdef PyUFuncGenericFunction q1_loop_func[1]
cdef char q1_input_output_types[2]
cdef void *q1_elementwise_funcs[1]
q1_loop_func[0] = PyUFunc_d_d
q1_input_output_types[0] = NPY_DOUBLE
q1_input_output_types[1] = NPY_DOUBLE
q1_elementwise_funcs[0] = <void*>q1_single_point
q1_py = PyUFunc_FromFuncAndData(
    q1_loop_func,
    q1_elementwise_funcs,
    q1_input_output_types,
    1, # number of supported input types
    1, # number of input args
    1, # number of output args
    0, # 'indentity' element, never mind this
    "q1_py", # function name
    "q1_py(x) -> computes int_{0.5}^x K_1(u)/u^12 du", # docstring
    0 # unused
)

# ufunc declaration of q2
cdef PyUFuncGenericFunction q2_loop_func[1]
cdef char q2_input_output_types[2]
cdef void *q2_elementwise_funcs[1]
q2_loop_func[0] = PyUFunc_d_d
q2_input_output_types[0] = NPY_DOUBLE
q2_input_output_types[1] = NPY_DOUBLE
q2_elementwise_funcs[0] = <void*>q2_single_point
q2_py = PyUFunc_FromFuncAndData(
    q2_loop_func,
    q2_elementwise_funcs,
    q2_input_output_types,
    1, # number of supported input types
    1, # number of input args
    1, # number of output args
    0, # 'indentity' element, never mind this
    "q2_py", # function name
    "q2_py(x) -> computes int_{0.5}^x K_1(u)/u^6 du", # docstring
    0 # unused
)




"""
cimport header

import numpy as np
cimport numpy as np

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t



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


"""
