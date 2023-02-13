#make file - this is a comment section
#gcc -O3 fftw1d.c fftw1d -lfftw3 -lm
#gcc -O3 adi1d.c -o adi1d -lgsl -lm
all:    #target name
	 gcc -O3 fftw2d.c -o fftw2d -lfftw3 -lm
	 ./fftw2d
	 gthumb refractive.png &