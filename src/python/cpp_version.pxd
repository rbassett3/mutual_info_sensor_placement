# distutils: language = c++
# 1 "cython_version.pyx"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "cython_version.pyx"
from libcpp.string cimport string;
cdef extern from 'cpp_version.h':
    cdef int cython_wrapper(double*, int, int, double, int, int, string, int)
