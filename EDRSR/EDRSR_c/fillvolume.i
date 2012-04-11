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

  static int ANz,ANy,ANx , Adim1,Adim2 = 0;


%}

%include typemaps.i


%init %{	
  import_array();
%}	



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
 ANz=$1 =tmp->dimensions[0] ;
 ANy=$2 =tmp->dimensions[1] ;
 ANx=$3 =tmp->dimensions[2] ;
}
%typemap(freearg) (int  Nz, int Ny, int Nx, float  * VOLUME) {Py_DECREF(tmp$argnum);}




%typemap(in) (float  * VOLUME_CHECK) (PyArrayObject* tmp=NULL) {

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
 $1 = (float   *)tmp->data;


 if(ANz!=tmp->dimensions[0] ||  ANy!=tmp->dimensions[1]  ||  ANx!=tmp->dimensions[2]  ) {
   PyErr_SetString(PyExc_ValueError, "Volumes  must be same length.");
   return NULL;
 }
 
}
%typemap(freearg) ( float  * VOLUME_CHECK) {Py_DECREF(tmp$argnum);}



%typemap(in) (int  dim1, int dim2,  float  * QFIN) (PyArrayObject* tmp=NULL) {

 if(!PyArray_Check($input)) {
    printf(" expected an array as argument \n");
    return NULL;	
 }
 if(!PyArray_ISCONTIGUOUS((PyArrayObject *) $input)) {
    PyErr_SetString(PyExc_ValueError," need contiguous array in typemap   int   dim1, int dim2,  float  * QFIN   \n ");
    return NULL;
 }
 
 if( ((PyArrayObject *) $input )->descr->type_num != PyArray_FLOAT )  {
   PyErr_SetString(PyExc_ValueError, "QFIN is not of type float.");
   return NULL;
 }

 // tmp = (PyArrayObject *) $input  ; 
 if(((PyArrayObject *) $input )->nd!=3)  {
   PyErr_SetString(PyExc_ValueError, "VOLUME is not 3D");
   return NULL;
 }
 
 tmp = (PyArrayObject *)PyArray_ContiguousFromObject($input, PyArray_FLOAT, 3, 3);
 // si on voulait  lacher le GIL c' est ici
 $3 = (float   *)tmp->data;
 Adim1=$1 =tmp->dimensions[0] ;
 Adim2=$2 =tmp->dimensions[1] ;

 if(3 !=tmp->dimensions[2]) {
   PyErr_SetString(PyExc_ValueError, "QFin  third dimension  must be 3.");
   return NULL;
 }
 
}
%typemap(freearg) (int  dim1, int dim2,  float  * QFIN) {Py_DECREF(tmp$argnum);}



%typemap(in) (float  * IMAGE_check ) (PyArrayObject* tmp=NULL) {

 if(!PyArray_Check($input)) {
    printf(" expected an array as argument \n");
    return NULL;	
 }
 if(!PyArray_ISCONTIGUOUS((PyArrayObject *) $input)) {
    PyErr_SetString(PyExc_ValueError," need contiguous array in typemap   int  Nz, int Ny, int Nx, double * VOLUME    \n ");
    return NULL;
 }
 
 if( ((PyArrayObject *) $input )->descr->type_num != PyArray_FLOAT )  {
   PyErr_SetString(PyExc_ValueError, "IMAGE_check is not of type float.");
   return NULL;
 }

 // tmp = (PyArrayObject *) $input  ; 
 if(((PyArrayObject *) $input )->nd!=2)  {
   PyErr_SetString(PyExc_ValueError, "IMAGE_check is not 2D");
   return NULL;
 }
 
 tmp = (PyArrayObject *)PyArray_ContiguousFromObject($input, PyArray_FLOAT, 2, 2);
 // si on voulait  lacher le GIL c' est ici
 $1 = (float   *)tmp->data;

 if(Adim1!=tmp->dimensions[0] ||  Adim2!=tmp->dimensions[1]    ) {
   PyErr_SetString(PyExc_ValueError, "IMAGE_check and QFIN  must be same length.");
   return NULL;
 }
 
}
%typemap(freearg) ( float  *IMAGE_check ) {Py_DECREF(tmp$argnum);}

void pippo(int  Nz, int Ny, int Nx, float  * VOLUME )  ; 

/*

	comment="""
	# Qfin float (dim1,dim2,3)
	#  Volume,                          float(nz,ny,nz   )
	#  Mask,                            float(nz,ny,nz   ) 
	#  C3,                              float32 (dim1,dim2) 
	#  POL_tmp                          float32 (dim1,dim2)  
	# (data,                            float32 (dim1,dim2) 

*/

void func_somme(	float  q0x ,float q0y ,float q0z,
	       float  dqx ,float dqy ,float dqz ,

		int  Nz, int Ny, int Nx, float  * VOLUME,  
		float *VOLUME_CHECK, 
		int dim1,int dim2, float *QFIN,
		float *IMAGE_check,
		float *IMAGE_check,
		float *IMAGE_check  ) ;
  

/* def func_somme(	float  q0x ,float q0y ,float q0z, */
/* 	       float  dqx ,float dqy ,float dqz , */
/* 	       int nz,int ny,int nx,float *Volume, */
/* 	       float *Mask,  */
/* 	       int dim1,int dim2, float *Qfin, */
/* 	       float *data, */
/* 	       float *POL_tmp, */
/* 	       float *C3  ) { */
