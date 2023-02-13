#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
/* plain gauss elimination, only not bothering with the zeroes
 *
 *       diag[0]  abovediag[0]             0   .....
 *  belowdiag[0]       diag[1]  abovediag[1]   .....
 *             0  belowdiag[1]       diag[2]
 *             0             0  belowdiag[2]   .....
 */
int solve_tridiag(complex diag[], complex abovediag[], complex belowdiag[], complex rhs[],
                                complex x[], size_t N){

  double *alpha = (double *) malloc (N * sizeof (double));
  double *z = (double *) malloc (N * sizeof (double));

  size_t i, j;

      /* Bidiagonalization (eliminating belowdiag)
         & rhs update
         diag' = alpha
         rhs' = z
       */
      alpha[0] = abovediag[0]/diag[0];
      z[0] = rhs[0]/diag[0];


      for(i = 1; i < N-2; i++){
          alpha[i] = abovediag[i] / (diag[i] - belowdiag[i-1] * alpha[i-1]);
          z[i] = (rhs[i] - belowdiag[i-1]*z[i-1])/(diag[i] - belowdiag[i-1]*alpha[i-1]);
        }
      i = N - 1;
      z[i] = (rhs[i] - belowdiag[i-1]*z[i-1])/(diag[i] - belowdiag[i-1]*alpha[i-1]);
      /* backsubstitution */
      x[N - 1] = z[N - 1];
      if (N >= 2)
        {
          for (i = N - 2, j = 0; j <= N - 2; j++, i--)
            {
              x[i] = z[i] - abovediag[i] * x[i+1];
            }
        }

  if (z != 0)
    free (z);
  if (alpha != 0)
    free (alpha);

  return 0;
}


int main(){    
  // points
  int Nx = 1000;
  int Nz = 1000;
  
  // parameters
  double n0 = 1.5;
  double l0 = 632E-9;
  double w = 4E-6;
  double Lx = 200E-6;

  double zmax = 1.6E-3;

  double dx = Lx/Nx;
  double dz = zmax/Nz;
  double k0 = 2*M_PI/l0;
  double beta = k0 * n0;

  double *n = malloc(Nx * sizeof(double));

  double alpha = 0.5;
  
  complex edx0;
  complex edxN;

  complex val;
  double xi;

  FILE *fp;

  complex *diag = (complex *) malloc(Nx * sizeof(complex));
  complex *abovediag = (complex *) malloc((Nx-1) * sizeof(complex));
  complex *belowdiag = (complex *) malloc((Nx-1) * sizeof(complex));
  complex *b = (complex *) malloc(Nx * sizeof(complex));
  complex *x = (complex *) malloc(Nx * sizeof(complex));
  complex *x_new = (complex *) malloc(Nx * sizeof(complex));

  complex **p, **q, *r;

  for(int i = 0; i < Nx; i++){
    n[i] = n0;
  }

  fp = freopen("output1d.txt", "w", stdout);
  
  // initial field (gaussian)
  for(int i = 0; i < Nx; i++){
    xi = -Lx/2 + i*dx;
    x[i] = exp(-(xi*xi)/(w*w));
    printf("%e %e %e\n", 0.0, -Lx/2 + i*dx, cabs(x[i]));
  }
  printf("\n");
  for(int k = 1; k < Nz; k++){
    int i = 0;
    // b
    val = (2*alpha / (dx*dx)) - alpha*(n[i]*n[i] - n0*n0)*(k0*k0) + 2*I*k0*n0/dz;      
    diag[i] = val;

    // c
    val = -alpha/(dx * dx);
    abovediag[i] = val;

    // r
    val = ((1-alpha)/(dx * dx)) * x[i+1]
        + (
          (1-alpha)*(n[i]*n[i] - n0*n0)*k0*k0
          - 2*(1-alpha)/(dx*dx)
          + 2*I*k0*n0/dz) * x[i];
    b[i] = val;
  
    for(i = 1; i < Nx-1; i++){
      // a
      val = -alpha/(dx * dx);
      belowdiag[i-1] = val;
      
      // b
      val = (2*alpha / (dx*dx)) - alpha*(n[i]*n[i] - n0*n0)*k0*k0 + 2*I*k0*n0/dz;      
      diag[i] = val;

      // c
      val = -alpha/(dx * dx);
      abovediag[i] = val;

      // r
      val = ((1-alpha)/(dx * dx)) * (x[i-1] + x[i+1])
          + (
            (1-alpha)*(n[i]*n[i] - n0*n0)*k0*k0
            - 2*(1-alpha)/(dx*dx)
            + 2*I*k0*n0/dz) * x[i];
      b[i] = val;
    }
    // a
    val = -alpha/(dx * dx);
    belowdiag[i-1] = val;
    
    // b
    val = (2*alpha / (dx*dx)) - alpha*(n[i]*n[i] - n0*n0)*k0*k0 + 2*I*k0*n0/dz;      
    diag[i] = val;

    // r
    val = ((1-alpha)/(dx*dx)) * x[i-1] 
        + (
          (1-alpha)*(n[i]*n[i] - n0*n0)*k0*k0
          - 2*(1-alpha)/(dx*dx)
          + 2*I*k0*n0/dz) * x[i];
    b[i] = val;

    solve_tridiag(diag, abovediag, belowdiag, b, x_new, Nx);

    edx0 = x[0]/x[1];
    x_new[0] = x_new[1]*edx0;

    edxN = x[Nx-1]/x[Nx-2];
    x_new[Nx-1] = x_new[Nx-2]*edxN;

    p = &x;
    q = &x_new;
    r = *p;
    *p = *q;
    *q = r; 
    
    //if((k % (Nz/10)) == 0){
    for(int i = 0; i < Nx; i++){
      printf("%e %e %e\n", k*dz, -Lx/2 + i*dx, cabs(x_new[i]));
    }
      printf("\n");
   //}
    
  }
  fclose(fp);
  free(diag);
  free(abovediag);
  free(belowdiag);
  free(b);
  free(x);
  free(x_new);
}
