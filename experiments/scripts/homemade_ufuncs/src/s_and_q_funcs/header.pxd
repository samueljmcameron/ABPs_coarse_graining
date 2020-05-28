# cython: language_level = 3

cdef extern from "header.h":

  extern double s1(double x) nogil;
  extern double s2(double x) nogil;
  extern double q1(double x) nogil;
  extern double q2(double x) nogil;
  
