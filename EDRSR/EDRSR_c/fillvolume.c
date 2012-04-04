
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
