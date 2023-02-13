// Saves into a text file, read and save with gnuplot the output of a gaussian light beam propagating in a 1D waveguide array

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

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

int main(){    
    // number of points in grid
    int Nx = 500;
    int Ny = 500;
    int Nz = 500;
    
    // parameters
    double n0 = 1.48; // refraction index of borosilicate 
    double l0 = 700E-9; // wavelenght of light
    double alpha = 11E-6; // width parameter in 2D Gaussian beam e^(-(x^2+y^2)/alpha^2)
    double wx = 5E-6; // width of the waveguide 
    double wy = wx * 1.5; // height of the waveguide
    double Lx = 300E-6; // width of the grid
    double Ly = 100E-6; // height of the grid

    double zmax = 20E-3; // propagation distance

    double super_gaussian_power = 3; // exponent of super-gaussian waveguide
    
    // auxiliar variables
    double dx = Lx/Nx;
    double dy = Ly/Ny;
    double dz = zmax/Nz;
    double k0 = 2*M_PI/l0;
    double beta = k0 * n0;
    double xi, yj;

    double* dn = malloc(sizeof(double) * Nx * Ny);
    
    // 1D array setup
    double dn1 = 4E-4; // contrast of first waveguide
    double dn2 = 4E-4; // contrast of first waveguide

    double d1x = 19E-6; // X separation of waveguides
    double d1y = 0; // Y separation of waveguides
    
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
            dn[i+Nx*j] = dn1 * exp(-(pow(fabs((xi-d1x/2)/wx), super_gaussian_power) + pow(fabs((yj-d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi+d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi-3*d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi+3*d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi-5*d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi+5*d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi-7*d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi+7*d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi-9*d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi+9*d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            dn[i+Nx*j] += dn2 * exp(-(pow(fabs((xi-11*d1x/2)/wx), super_gaussian_power) + pow(fabs((yj+d1y/2)/wy), super_gaussian_power)));
            
            printf("%e %e %e\n", xi, yj, dn[i+Nx*j]);
        }
    }
    fclose(fp1);

    // call gnuplot to plot the input data
    FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "set pm3d map\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [-%e:%e]\n", Lx/2, Lx/2);
    fprintf(gnuplotPipe, "set yrange [-%e:%e]\n", Ly/2, Ly/2);
    fprintf(gnuplotPipe, "set xlabel 'X (m)'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y (m)'\n");
    fprintf(gnuplotPipe, "set title 'Refractive Index'\n");
    fprintf(gnuplotPipe, "splot 'refractive2d.txt' using 1:2:3 with image\n");
    fprintf(gnuplotPipe, "set term png\n");
    fprintf(gnuplotPipe, "set output 'refractive.png'\n");
    fprintf(gnuplotPipe, "replot'\n");
    fprintf(gnuplotPipe, "set term x11'\n");
    fflush(gnuplotPipe);
    
    // save the input (gaussian) in a text file
    fp2 = freopen("input2d.txt", "w", stdout);
    
    // initial field (gaussian)
    for(i = 0; i < Nx; i++){
        for(j = 0; j < Ny; j++){
            xi = -0.5*Lx + i*dx;
            yj = -0.5*Ly + j*dy;
            in[i+Nx*j] = exp(-( (xi-d1x/2)*(xi-d1x/2) + (yj-d1y/2)*(yj-d1y/2) ) / (alpha*alpha));
            printf("%e %e %e\n", xi, yj, cabs(in[i+Nx*j]) * cabs(in[i+Nx*j]));
        }
        printf("\n");
    }
    fclose(fp2);

    // call gnuplot to plot the input data
    gnuplotPipe = popen("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "set pm3d map\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [-%e:%e]\n", Lx/2, Lx/2);
    fprintf(gnuplotPipe, "set yrange [-%e:%e]\n", Ly/2, Ly/2);
    fprintf(gnuplotPipe, "set xlabel 'X (m)'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y (m)'\n");
    fprintf(gnuplotPipe, "set title 'Initial Gaussian Field'\n");
    fprintf(gnuplotPipe, "splot 'input2d.txt' using 1:2:3 with image\n");
    fprintf(gnuplotPipe, "set term png\n");
    fprintf(gnuplotPipe, "set output 'input.png'\n");
    fprintf(gnuplotPipe, "replot'\n");
    fprintf(gnuplotPipe, "set term x11'\n");
    fflush(gnuplotPipe);

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
    
    // save the output result of propagation in a text file
    fp3 = freopen("output2d.txt", "w", stdout);


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
    
    
    fclose(fp3);

    // call gnuplot to plot the output data
    gnuplotPipe = popen("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "set pm3d map\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [-%e:%e]\n", Lx/2, Lx/2);
    fprintf(gnuplotPipe, "set yrange [-%e:%e]\n", Ly/2, Ly/2);
    fprintf(gnuplotPipe, "set xlabel 'X (m)'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y (m)'\n");
    fprintf(gnuplotPipe, "set title 'Propagated Gaussian Field'\n");
    fprintf(gnuplotPipe, "splot 'output2d.txt' using 1:2:3 with image\n");
    fprintf(gnuplotPipe, "set term png\n");
    fprintf(gnuplotPipe, "set output 'output.png'\n");
    fprintf(gnuplotPipe, "replot'\n");
    fprintf(gnuplotPipe, "set term x11'\n");
    fflush(gnuplotPipe);
    
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(aux); 
    fftw_free(out);
}
