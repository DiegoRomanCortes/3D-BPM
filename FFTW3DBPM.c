// Saves into a text file the output of a gaussian light beam propagating in a 2D waveguide array

/*
Copyright (C) 2023  Diego Roman-Cortes

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <gsl/gsl_sf_laguerre.h>

void tostring(char str[], int num){
    int i, rem, len = 0, n;

    n = num;
    while(n != 0){
        len++;
        n /= 10;
    }
    for(i = 0; i < len; i++){
        rem = num % 10;
        num = num /10;
        str[len - (i + 1)] = rem + '0';
    }
    str[len] = '\0';
}

int main(){    
    // number of points in grid
    int Nx = 1500;
    int Ny = 1500;
    int Nz = 500;
    
    // parameters
    double n0 = 1.48; // refraction index of borosilicate 
    double l0 = 730E-9; // wavelenght of light
    double alpha = 6E-6; // width parameter in 2D Gaussian beam e^(-(x^2+y^2)/alpha^2)
    double wx = 2.1E-6; // width of the waveguide 
    double wy = 3.0E-6; // height of the waveguide
    double sigma = 2.05E-6; // width of LG-mode
    double l = 1; // azimuthal parameter of LG-mode
    int n = 0; // radial parameter of LG-mode
    double Lx = 250E-6; // width of the grid
    double Ly = 250E-6; // height of the grid

    double zmax = 5E-3; // propagation distance

    double super_gaussian_power = 3; // exponent of super-gaussian waveguide
    
    // auxiliar variables
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    double dz = zmax/Nz;
    double k0 = 2*M_PI/l0;
    double beta = k0 * n0;
    double xi, yj, r, phi;

    double* dn = malloc(sizeof(double) * Nx * Ny);
    
    // 1D array setup
    double dn1 = 30E-4; // contrast of first waveguide
    double dn2 = 4E-4; // contrast of other waveguides

    double d1x = 20E-6; // X separation of waveguides
    double d1y = 20E-6; // Y separation of waveguides
    
    // for animation
    int frames = 60;
    int rem, div;
    char strdiv[10];

    int i, j, k;

    fftw_complex *in, *out;
    fftw_plan p;

    FILE *fp1, *fp2, *fp3;

    //initialization of FFTW
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx * Ny);
    p = fftw_plan_dft_2d(Nx, Ny, in, out, FFTW_BACKWARD, FFTW_MEASURE);
    
    fp1 = freopen("refractive2d.txt", "w", stdout);
    // shape of refractive index contrast
    for(i = 0; i < Nx; i++){
        for(j = 0; j < Ny; j++){
            xi = -0.5*Lx + i*dx;
            yj = -0.5*Ly + j*dy;
            // dn[i+Nx*j] = dn1 * exp(-(pow(fabs((xi)/wx), super_gaussian_power) + pow(fabs((yj)/wy), super_gaussian_power)));
            dn[i+Nx*j] = dn1 * tanh(33 / (exp((xi/wx)*(xi/wx) + (yj/wy)*(yj/wy))));
            // for(int n = 0; n < 11; n++){
            //     for(int m = 0; m < 11; m++){
            //         dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi-n*d1x+10*0.5*d1x)/wx), super_gaussian_power) + pow(fabs((yj-m*d1x+10*0.5*d1x)/wy), super_gaussian_power)));
            //     }
            // }
            printf("%e %e %e\n", xi, yj, dn[i+Nx*j]);
        }
    }
    fclose(fp1);
    
    // save the input (gaussian) in a text file
    fp2 = freopen("0.txt", "w", stdout);
    
    // initial field (gaussian)
    for(i = 0; i < Nx; i++){
        for(j = 0; j < Ny; j++){
            xi = -0.5*Lx + i*dx;
            yj = -0.5*Ly + j*dy;
            r = sqrt(xi*xi + yj*yj);
            phi = atan2(yj, xi);
            in[i+Nx*j] = cexp(-r*r/(sigma*sigma)) * pow(r/sigma, l) *  gsl_sf_laguerre_n(n, l, 2*r*r/(sigma*sigma)) * cexp(I*phi*l); // laguerre-gaussian mode
            printf("%e %e %e\n", xi, yj, cabs(in[i+Nx*j]) * cabs(in[i+Nx*j]));
        }
        printf("\n");
    }
    fclose(fp2);

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
            }
        }

        // save to txt
        rem = k % (Nz/frames);
        if(rem == 0){
            div = k / (Nz/frames);
            tostring(strdiv, div);
            fp3 = freopen(strcat(strdiv, ".txt"), "w", stdout);
            for(i = 0; i < Nx; i++){
                for(j = 0; j < Ny; j++){
                    printf("%e %e %e\n", -0.5*Lx + i*dx, -0.5*Ly + j*dy, cabs(in[i+j*Nx]) * cabs(in[i+j*Nx]));
                }
                printf("\n");
            }
            fclose(fp3);
        }

    }
    
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(aux); 
    fftw_free(out);

}