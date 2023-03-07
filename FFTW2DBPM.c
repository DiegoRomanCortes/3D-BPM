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

#include <math.h>
#include <complex.h>
#include <fftw3.h>

int main(){    
    // points
    int Nx = 10;
    int Nz = 10;
    
    // parameters
    double n0 = 1.5;
    double l0 = 632E-9;
    double alpha = 4E-6;
    double Lx = 200E-6;

    double zmax = 1.6E-3;

    double dx = Lx/Nx;
    double dz = zmax/Nz;
    double k0 = 2*M_PI/l0;
    double beta = k0 * n0;
    double xi;
    
    int i;

    fftw_complex *in, *out;
    fftw_plan p;

    FILE *fp;

    //...
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);
    p = fftw_plan_dft_1d(Nx, in, out, FFTW_BACKWARD, FFTW_MEASURE);
    fp = freopen("output1d.txt", "w", stdout);
    // initial field (gaussian)
    for(i = 0; i < Nx; i++){
        xi = -Lx/2 + i*dx;
        in[i] = exp(-(xi*xi)/(alpha*alpha));
        printf("%e %e %e\n", 0.0, -Lx/2 + i*dx, cabs(in[i]));
    }
    printf("\n");

    // frequency index
    int freqidx[Nx];

    for(i=0; i < Nx/2; i++){
        freqidx[i] = i;
    }

    for(i=Nx/2; i < Nx; i++){
        freqidx[i] =  i-Nx;
    }

    fftw_complex phase[Nx];
    for(i=0; i < Nx; i++){
        phase[i] = cexp(I*dz*(2*M_PI*freqidx[i]/Lx)*(2*M_PI*freqidx[i]/Lx)/(2*beta));
    }

    // main loop
    fftw_complex *aux = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nx);

    for(i=1; i < Nz; i++){
        fftw_execute(p);

        p = fftw_plan_dft_1d(Nx, aux, in, FFTW_FORWARD, FFTW_MEASURE);
        
        for(int j=0; j < Nx; j++){
            aux[j] = out[j] * phase[j];
        }

        fftw_execute(p);
        for(int j=0; j < Nx; j++){
            in[j] /= Nx;
            printf("%e %e %e\n", dz*i, -Lx/2 + j*dx, cabs(in[j]));
        }
        printf("\n");
        p = fftw_plan_dft_1d(Nx, in, out, FFTW_BACKWARD, FFTW_MEASURE);
    }

    fclose(fp);
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(aux); 
    fftw_free(out);
}
