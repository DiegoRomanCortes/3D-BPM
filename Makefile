#make file - this is a comment section
#gcc -O3 fftw1d.c fftw1d -lfftw3 -lm
#gcc -O3 adi1d.c -o adi1d -lgsl -lm
all:    #target name
	 gcc -O3 FFTW3DBPM.c -o FFTW3DBPM -lfftw3 -lm
	 ./FFTW3DBPM
	gthumb refractive.png &
clear:
	rm *.txt
	rm *.png
	rm FFTW3DBPM