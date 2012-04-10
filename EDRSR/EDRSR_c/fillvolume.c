#include<emmintrin.h>
#define FLOAT_TO_INT(in,out)  \
     out=_mm_cvtss_si32(_mm_load_ss(&(in)));

void pippo(int nz, int ny, int nx,     float *volume ) {
  int i,j,k;
  for(i=0; i<nz; i++) {
    for(j=0; j<ny; j++) {
      for(k=0; k<nx; k++) {
	volume[ (i*ny+j)*nx +k] = 1234.987 ; 
      }
    }
  }
}

/*
	comment="""
	# Qfin float (dim1,dim2,3)
	#  Volume,                          float(nz,ny,nz   )
	#  Mask,                            float(nz,ny,nz   ) 
	#  C3,                              float32 (dim1,dim2) 
	#  POL_tmp                          float32 (dim1,dim2)  
	# (data,                            float32 (dim1,dim2) 

*/

def func_somme(	float  q0x ,float q0y ,float q0z,
	       float  dqx ,float dqy ,float dqz ,
	       int nz,int ny,int nx,float *Volume,
	       float *Mask, 
	       int dim1,int dim2, float *Qfin,
	       float *data,
	       float *POL_tmp,
	       float *C3  ) {

  int i,j,k;
  int nz_2, ny_2, nx_2;
  int npoints, ipoint, icubdim_2;

  npoints = dim1*dim2;

  nx_2 = (nx-1)/ 2; 
  ny_2 = (ny-1)/ 2; 
  nz_2 = (nz-1)/ 2; 

  printf( "RECIPROCAL SPACE CENTER  =%e %e %e", q0x, q0y, q0z);
  printf('-----------------------------\n');
  printf( "Computation of 3D Volume Indices ...\n");
  
  //  time9 = time.time()
    
  for(ipoint=0; ipoint<npoints; ipoint++) {
    FLOAT_TO_INT( (  (nx_2  * (1.0 + ( Qfin[ipoint*3+0]-q0x)/dqx ))),  i  );
    FLOAT_TO_INT( (  (ny_2  * (1.0 + ( Qfin[ipoint*3+1]-q0y)/dqy ))),  j  );
    FLOAT_TO_INT( (  (nz_2  * (1.0 + ( Qfin[ipoint*3+2]-q0z)/dqz ))),  k  );
    if(i<=nx && j<=ny && k<=nz ) {
      Volume[ (k*ny +j)*nx +i ] +=  data[ipoint]/(POL_tmp[ipoint] *C3[ipoint]);  
      Mask[ (k*ny +j)*nx +i ] += 1.0;
    }
}
