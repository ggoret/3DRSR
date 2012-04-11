

void pippo(int nz, int ny, int nx,     float *volume )  ; 

void func_somme(float  q0x ,float q0y ,float q0z,
	       float  dqx ,float dqy ,float dqz ,
	       int nz,int ny,int nx,float *Volume,
	       float *Mask, 
	       int dim1,int dim2, 
	       float *Qfin,
	       float *data,
	       float *POL_tmp,
	       float *C3  );
