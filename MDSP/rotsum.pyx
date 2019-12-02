# -*- coding: utf-8 -*-
###################################################################################
# rotation and summation cython module of MDSP
# Gael Goret for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr
###################################################################################
#distutils: extra_compile_args = -fopenmp
#distutils: extra_link_args = -fopenmp

import cython
from cython.parallel cimport prange
from cpython cimport bool
cimport numpy
from numpy cimport ndarray
cdef extern from "math.h":
    double  fabs(float)nogil
import numpy

#-----------------------------------------------------------------------------------                       
                    
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def orient_volume(ndarray[numpy.float32_t, ndim = 3] volume not None,
                  ndarray[numpy.float32_t, ndim = 3] mask not None,
                  ndarray[numpy.float32_t, ndim = 2] rotation_array not None):
    cdef int  idim, jdim, kdim,i,j,k,p,q,r
    cdef float m,n,o, intensity, icenter, jcenter, kcenter
    idim = volume.shape[0]
    jdim = volume.shape[1]
    kdim = volume.shape[2]
    cdef ndarray[numpy.float32_t, ndim = 3] out = numpy.zeros_like(volume)
    cdef ndarray[numpy.float32_t, ndim = 3] cpt = numpy.zeros_like(volume)
    icenter = (idim-1)/2.0 
    jcenter = (jdim-1)/2.0 
    kcenter = (kdim-1)/2.0 
    with nogil:
        for i in prange(idim):
            m=i - icenter
            for j in range(jdim):
                n=j - jcenter
                for k in range(kdim):
                    o=k - kcenter
                    p = <int> (rotation_array[0,0]*m +rotation_array[0,1]*n+rotation_array[0,2]*o + icenter) 
                    q = <int> (rotation_array[1,0]*m +rotation_array[1,1]*n+rotation_array[1,2]*o + jcenter)
                    r = <int> (rotation_array[2,0]*m +rotation_array[2,1]*n+rotation_array[2,2]*o + kcenter)                    
                    if (0<=p<idim) and (0<=q<jdim) and (0<=r<kdim) and mask[p,q,r]:
                        intensity = volume[p,q,r]
                        out[i,j,k] += intensity 
                        cpt[i,j,k] += <int>(intensity!=0)
    return out,cpt
    
#-----------------------------------------------------------------------------------           
                    
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def orient_volume_opt(ndarray[numpy.float32_t, ndim = 3] volume not None,
                  ndarray[numpy.float32_t, ndim = 3] mask not None,
                  ndarray[numpy.float32_t, ndim = 2] rotation_array not None):
    cdef int  idim, jdim, kdim,i,j,k,p,q,r
    cdef float m,n,o, intensity, icenter, jcenter, kcenter
    idim = volume.shape[0]
    jdim = volume.shape[1]
    kdim = volume.shape[2]
    cdef ndarray[numpy.float32_t, ndim = 3] out = numpy.zeros_like(volume)
    cdef ndarray[numpy.float32_t, ndim = 3] cpt = numpy.zeros_like(volume)
    icenter = (idim-1)/2.0 
    jcenter = (jdim-1)/2.0 
    kcenter = (kdim-1)/2.0 
    cdef float r00, r01, r02, r10, r11, r12, r20, r21, r22
    r00 = rotation_array[0,0]
    r01 = rotation_array[0,1]
    r02 = rotation_array[0,2]
    r10 = rotation_array[1,0]
    r11 = rotation_array[1,1]
    r12 = rotation_array[1,2]
    r20 = rotation_array[2,0]
    r21 = rotation_array[2,1]
    r22 = rotation_array[2,2]
    with nogil:
        for i in prange(idim):
            m=i - icenter
            for j in range(jdim):
                n=j - jcenter
                for k in range(kdim):
                    o=k - kcenter
                    p = <int> (r00*m +r01*n+r02*o + icenter) 
                    q = <int> (r10*m +r11*n+r12*o + jcenter)
                    r = <int> (r20*m +r21*n+r22*o + kcenter)                    
                    if (0<=p<idim) and (0<=q<jdim) and (0<=r<kdim) and mask[p,q,r]:
                        intensity = volume[p,q,r]
                        out[i,j,k] += intensity 
                        cpt[i,j,k] += <int>(intensity!=0)
    return out,cpt
    
#-----------------------------------------------------------------------------------                       
                    
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def orient_image(ndarray[numpy.float32_t, ndim = 2] img not None,
                  ndarray[numpy.float32_t, ndim = 2] mask not None,
                  ndarray[numpy.float32_t, ndim = 2] rotation_array not None):
    cdef int  idim, jdim, i, j, p, q
    cdef float m, n, intensity, icenter, jcenter
    idim = img.shape[0]
    jdim = img.shape[1]
    cdef ndarray[numpy.float32_t, ndim = 2] out = numpy.zeros_like(img)
    cdef ndarray[numpy.float32_t, ndim = 2] cpt = numpy.zeros_like(img)
    icenter = (idim-1)/2.0 
    jcenter = (jdim-1)/2.0 
    with nogil:
        for i in prange(idim):
            m = i - icenter
            for j in range(jdim):
                n = j - jcenter
                p = <int> (rotation_array[0,0]*m +rotation_array[0,1]*n + icenter) 
                q = <int> (rotation_array[1,0]*m +rotation_array[1,1]*n + jcenter)
                if (0<=p<idim) and (0<=q<jdim) and mask[p,q]:
                    intensity = img[p,q]
                    out[i,j] += intensity 
                    cpt[i,j] += <int>(intensity!=0)
    return out,cpt
