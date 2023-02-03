// Saves into a text file the output of a gaussian light beam propagating in a waveguide

/*
Copyright (C) <year>  <name of author>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

int main(){    
    // points
    int Nx = 200;
    int Ny = 200;
    int Nz = 200;
    
    // parameters
    double n0 = 1.48;
    double l0 = 685E-9;
    double alpha = 11E-6;
    double wx = 10E-6;
    double wy = wx * 1.6;
    double Lx = 100E-6;
    double Ly = 100E-6;

    double zmax = 5E-3;//1.6E-3;

    double super_gaussian_power = 2.5;

    double dx = Lx/Nx;
    double dy = Ly/Ny;
    double dz = zmax/Nz;
    double k0 = 2*M_PI/l0;
    double beta = k0 * n0;
    double xi, yj;

    double* dn = malloc(sizeof(double) * Nx * Ny);
    double dn1 = 4E-4;
    double dn2 = dn1;

    double d1x = 0;
    double d1y = 18E-6;
    
    int i, j, k;

    fftw_complex *in, *out;
    fftw_plan p;

    FILE *fp1, *fp2;

    //initialization of FFTW
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    p = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_BACKWARD, FFTW_MEASURE);
    
    // shape of refractive index contrast
    for(i = 0; i < Nx; i++){
        for(j = 0; j < Ny; j++){
            xi = -Lx/2 + i*dx;
            yj = -Ly/2 + j*dy;
            dn[i+Nx*j] = dn1 * exp(-(pow(((xi-d1x)*(xi-d1x))/(wx*wx), super_gaussian_power) + pow(((yj-d1y)*(yj-d1y))/(wy*wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(((xi+d1x)*(xi+d1x))/(wx*wx), super_gaussian_power) + pow(((yj+d1y)*(yj+d1y))/(wy*wy), super_gaussian_power)));
        }
    }
    
    fp1 = freopen("input2d.txt", "w", stdout);
    
    // initial field (gaussian)
    for(i = 0; i < Nx; i++){
        for(j = 0; j < Ny; j++){
            xi = -Lx/2 + i*dx;
            yj = -Ly/2 + j*dy;
            in[i+Nx*j] = exp(-( (xi-d1x)*(xi-d1x) + (yj-d1y)*(yj-d1y) ) / (alpha*alpha));
            printf("%e %e %e\n", -Lx/2 + i*dx, -Ly/2 + j*dy, cabs(in[i+Nx*j]) * cabs(in[i+Nx*j]));
        }
        printf("\n");
    }
    fclose(fp1);


    // frequency indices
    int freqidx[Nx + Ny];

    for(i=0; i < Nx/2; i++){
        freqidx[i] = i;
    }
    for(j=0; j < Ny/2; j++){   
        freqidx[Nx+j] = j;
    }
    for(i=Nx/2; i < Nx; i++){
        freqidx[i] =  i-Nx;
    }
    for(j=Ny/2; j < Ny; j++){
        freqidx[Nx+j] =  j-Ny;
    }    
    

    fftw_complex *phase = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    for(i=0; i < Nx; i++){
        for(j=0; j < Ny; j++){
            phase[i+j*Nx] = cexp(I*dz*( (2*M_PI) * (2*M_PI) * ( (freqidx[i]/Lx) * (freqidx[i]/Lx) + (freqidx[Nx+j]/Ly) * (freqidx[Nx+j]/Ly) )/(4*beta)));
        }
    }
    
    fp2 = freopen("output2d.txt", "w", stdout);


    // main loop
    fftw_complex *aux = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);

    for(k=1; k < Nz; k++){
        fftw_execute(p); // 'out' now points towards the DFT of 'in' 

        p = fftw_plan_dft_2d(Nx, Ny, aux, in, FFTW_FORWARD, FFTW_MEASURE); // setup the inverse DFT
        
        for(i = 0; i < Nx; i++){
            for(j = 0; j < Ny; j++){
                aux[i+j*Nx] = out[i+j*Nx] * phase[i+j*Nx];
            }
        }

        fftw_execute(p); // 'in' now points towards the inverse DFT of 'aux'
        p = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_BACKWARD, FFTW_MEASURE); // setup the second DFT
        for(i = 0; i < Nx; i++){
            for(j = 0; j < Ny; j++){
                in[i+j*Nx] /= (Nx * Ny); // normalization of FFT
                in[i+j*Nx] *= cexp(-I * k0 * dn[i+j*Nx] * dz); // potential operator in real space
            }
        }
        
        fftw_execute(p); // 'out' now points towards the DFT of 'in' 

        p = fftw_plan_dft_2d(Nx, Ny, aux, in, FFTW_FORWARD, FFTW_MEASURE); // setup the inverse DFT

        for(i = 0; i < Nx; i++){
            for(j = 0; j < Ny; j++){
                aux[i+j*Nx] = out[i+j*Nx] * phase[i+j*Nx];
            }
        }

        fftw_execute(p); // 'in' now points towards the inverse DFT of 'aux'
        p = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_BACKWARD, FFTW_MEASURE);
        for(i = 0; i < Nx; i++){
            for(j = 0; j < Ny; j++){
                in[i+j*Nx] /= (Nx * Ny); // normalization of FFT
                if(k == (Nz - 1)){
                    printf("%e %e %e\n", -Lx/2 + i*dx, -Ly/2 + j*dy, cabs(in[i+j*Nx]) * cabs(in[i+j*Nx]));
                }
            }
            if(k == (Nz - 1)){
                printf("\n");
            }
        }
    }
    
    
    fclose(fp2);
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(aux); 
    fftw_free(out);
}
