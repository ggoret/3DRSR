# -*- coding: utf-8 -*-
###################################################################################
# Projection cython module of 3D Reciprocal Space Reconstruction
# Gael Goret and Jerome Kieffer for the European Synchrotron Radiation Facility
# gael.goret@esrf.fr
###################################################################################

import cython
from cython.parallel cimport prange
from cpython cimport bool
cimport numpy
from numpy cimport ndarray
cdef extern from "math.h":
    double  fabs(float)nogil
import numpy

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
def func_somme(float q0x, float q0y, float q0z,
           float  dqx, float dqy,float dqz,
           ndarray[numpy.float32_t, ndim = 3] Volume, 
           ndarray[numpy.float32_t, ndim = 3] Mask, 
           ndarray[numpy.float32_t, ndim = 3] Qfin,
           ndarray[numpy.float32_t, ndim = 2] data,
           ndarray[numpy.float32_t, ndim = 2] POL_tmp,
           ndarray[numpy.float32_t, ndim = 2] C3,
           ndarray[numpy.float32_t, ndim = 2] Filter,
           int sym_3_fold_axis,
           ndarray[numpy.float32_t, ndim = 3] sym_3_fold_ops):
    cdef int nz, ny, nx
    cdef int dim1, dim2,dimsym
    cdef int i,j,k,l,m,n
    cdef int nz_2, ny_2, nx_2
    cdef float pv, one, qfx, qfy, qfz,rotqfx, rotqfy, rotqfz
    nx   = Volume.shape[0]
    ny   = Volume.shape[1]
    nz   = Volume.shape[2]
    dim1 = data.shape[0]
    dim2 = data.shape[1]
    dimsym = sym_3_fold_ops.shape[0]
    nx_2 = (nx-1)/ 2
    ny_2 = (ny-1)/ 2 
    nz_2 = (nz-1)/ 2
    one = 1.0 
    with nogil:
        for l in prange(dim1):
            for m in range(dim2):
                if Filter[l,m]:
                    qfx = Qfin[l,m,0]
                    qfy = Qfin[l,m,1]
                    qfz = Qfin[l,m,2]
                    if sym_3_fold_axis:
                        for n in range(dimsym):
                            rotqfx = (sym_3_fold_ops[n,0,0]*qfx +sym_3_fold_ops[n,0,1]*qfy+sym_3_fold_ops[n,0,2]*qfz) 
                            rotqfy = (sym_3_fold_ops[n,1,0]*qfx +sym_3_fold_ops[n,1,1]*qfy+sym_3_fold_ops[n,1,2]*qfz)
                            rotqfz = (sym_3_fold_ops[n,2,0]*qfx +sym_3_fold_ops[n,2,1]*qfy+sym_3_fold_ops[n,2,2]*qfz)
                            i = <int> (nx_2  * (one + (rotqfx-q0x)/dqx ))
                            j = <int> (ny_2  * (one + (rotqfy-q0y)/dqy ))
                            k = <int> (nz_2  * (one + (rotqfz-q0z)/dqz ))
                            if (0<=i<nx) and (0<=j<ny) and (0<=k<nz):
                                pv = data[l,m]/(POL_tmp[l,m] * C3[l,m]) 
                                Volume[i,j,k] += pv    
                                Mask[i,j,k] += <int>(pv!=0)
                    i = <int> (nx_2  * (one + (qfx-q0x)/dqx ))
                    j = <int> (ny_2  * (one + (qfy-q0y)/dqy ))
                    k = <int> (nz_2  * (one + (qfz-q0z)/dqz ))
                    if (i<nx) and (j<ny) and (k<nz ):
                        pv = data[l,m]/(POL_tmp[l,m] *C3[l,m]) 
                        Volume[i,j,k] += pv    
                        Mask[i,j,k] += <int>(pv!=0)
                    
#-----------------------------------------------------------------------------------   
                    
@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)                    
def slice_somme(float dQ0, float dQ1, float dQ2,
           ndarray[numpy.float32_t, ndim = 1] Qoff, 
           ndarray[numpy.float32_t, ndim = 2] Image, 
           ndarray[numpy.float32_t, ndim = 2] Mask, 
           ndarray[numpy.float32_t, ndim = 3] Qn,
           ndarray[numpy.float32_t, ndim = 2] data,
           ndarray[numpy.float32_t, ndim = 2] POL_tmp,
           ndarray[numpy.float32_t, ndim = 2] C3,
           ndarray[numpy.float32_t, ndim = 2] Filter,
           int apply_sym,
           ndarray[numpy.float32_t, ndim = 3] sym_ops,
           ndarray[numpy.float32_t, ndim = 2] G):
    cdef float size, one,half_size
    cdef int dim1, dim2
    cdef int i,j,l,m
    cdef float Qxn_minus_Qxoff,Qyn_minus_Qyoff,Qzn_minus_Qzoff, rotQxn_minus_Qxoff, rotQyn_minus_Qyoff, rotQzn_minus_Qzoff, Qxoff,Qyoff,Qzoff,Qxn,Qyn,Qzn,pv,rotQxn,rotQyn,rotQzn,twodQ1,twodQ2
    size   = Image.shape[0]
    dim1 = data.shape[0]
    dim2 = data.shape[1]
    dimsym = sym_ops.shape[0]
    half_size = size/2.
    twodQ1 = 2.*dQ1
    twodQ2 = 2.*dQ2
    one = 1.0 
    with nogil:
        for l in prange(dim1):
            for m in range(dim2):
                if Filter[l,m]:
                    Qxoff = Qoff[0]
                    Qyoff = Qoff[1]
                    Qzoff = Qoff[2]
                    Qxn = Qn[l,m,0]
                    Qyn = Qn[l,m,1]
                    Qzn = Qn[l,m,2]
                    if apply_sym:
                        for n in range(dimsym):
                            rotQxn = sym_ops[n,0,0]*Qxn + sym_ops[n,0,1]*Qyn + sym_ops[n,0,2]*Qzn
                            rotQyn = sym_ops[n,1,0]*Qxn + sym_ops[n,1,1]*Qyn + sym_ops[n,1,2]*Qzn
                            rotQzn = sym_ops[n,2,0]*Qxn + sym_ops[n,2,1]*Qyn + sym_ops[n,2,2]*Qzn
                            rotQxn_minus_Qxoff = (G[0,0]*rotQxn + G[0,1]*rotQyn + G[0,2]*rotQzn) - Qxoff
                            rotQyn_minus_Qyoff = (G[1,0]*rotQxn + G[1,1]*rotQyn + G[1,2]*rotQzn) - Qyoff
                            rotQzn_minus_Qzoff = (G[2,0]*rotQxn + G[2,1]*rotQyn + G[2,2]*rotQzn) - Qzoff
                            if  fabs(rotQxn_minus_Qxoff)<=dQ1 and fabs(rotQyn_minus_Qyoff)<=dQ2 and fabs(rotQzn_minus_Qzoff)<=dQ0:
                                i = <int> (( size * (one + ( rotQxn_minus_Qxoff )/(twodQ1)) ) - half_size )
                                j = <int> (( size * (one + ( rotQyn_minus_Qyoff )/(twodQ2)) ) - half_size )
                                if (0<=i<size) and (0<=j<size):
                                    pv = data[l,m]/(POL_tmp[l,m] * C3[l,m])
                                    Image[i,j] += pv 
                                    Mask[i,j] +=  <int>(pv!=0)
                    else:
                        Qxn_minus_Qxoff = (G[0,0]*Qxn + G[0,1]*Qyn + G[0,2]*Qzn) - Qxoff
                        Qyn_minus_Qyoff = (G[1,0]*Qxn + G[1,1]*Qyn + G[1,2]*Qzn) - Qyoff
                        Qzn_minus_Qzoff = (G[2,0]*Qxn + G[2,1]*Qyn + G[2,2]*Qzn) - Qzoff
                        if  fabs(Qxn_minus_Qxoff)<=dQ1 and fabs(Qyn_minus_Qyoff)<=dQ2 and fabs(Qzn_minus_Qzoff)<=dQ0:
                            i = <int> (( size * (one + ( Qxn_minus_Qxoff )/(twodQ1)) ) - half_size )
                            j = <int> (( size * (one + ( Qyn_minus_Qyoff )/(twodQ2)) ) - half_size )
                            if (0<=i<size) and (0<=j<size):
                                pv = data[l,m]/(POL_tmp[l,m] * C3[l,m])
                                Image[i,j] += pv 
                                Mask[i,j] +=  <int>(pv!=0)

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
                    
