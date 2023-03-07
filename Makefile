#make file - this is a comment section
#gcc -O3 fftw1d.c fftw1d -lfftw3 -lm
#gcc -O3 adi1d.c -o adi1d -lgsl -lm
all:    #target name
	gcc -O3 FFTW3DBPM.c -o FFTW3DBPM -lgsl -lgslcblas -lfftw3 -lm
	./FFTW3DBPM
	gnuplot plot.gnu
	gthumb 0.png &
clear:
	rm *.png
	rm FFTW3DBPM
debug:
	gcc -g -o FFTW3DBPM FFTW3DBPM.c -lgsl -lgslcblas -lfftw3 -lm
	ddd FFTW3DBPM &