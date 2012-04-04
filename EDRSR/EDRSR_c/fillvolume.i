%module fillvolume
%{
#include<string.h>
#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include"fillvolume.h"
%}

%{
  ///   #include <numpy/oldnumeric.h>
#include<numpy/arrayobject.h>
%}

%include typemaps.i


%typemap(in) (int  Nz, int Ny, int Nx, float  * VOLUME) (PyArrayObject* tmp=NULL) {

 if(!PyArray_Check($input)) {
    printf(" expected an array as argument \n");
    return NULL;	
 }
 if(!PyArray_ISCONTIGUOUS((PyArrayObject *) $input)) {
    PyErr_SetString(PyExc_ValueError," need contiguous array in typemap   int  Nz, int Ny, int Nx, double * VOLUME    \n ");
    return NULL;
 }
 
 if( ((PyArrayObject *) $input )->descr->type_num != PyArray_FLOAT )  {
   PyErr_SetString(PyExc_ValueError, "VOLUME is not of type float.");
   return NULL;
 }

 // tmp = (PyArrayObject *) $input  ; 
 if(((PyArrayObject *) $input )->nd!=3)  {
   PyErr_SetString(PyExc_ValueError, "VOLUME is not 3D");
   return NULL;
 }
 
 tmp = (PyArrayObject *)PyArray_ContiguousFromObject($input, PyArray_FLOAT, 3, 3);
 // si on voulait  lacher le GIL c' est ici
 $4 = (float   *)tmp->data;
 $1 =tmp->dimensions[0] ;
 $2 =tmp->dimensions[1] ;
 $3 =tmp->dimensions[2] ;
}
%typemap(freearg) (int  Nz, int Ny, int Nx, float  * VOLUME) {Py_DECREF(tmp$argnum);}


void pippo(int  Nz, int Ny, int Nx, float  * VOLUME )  ; 

