# 2D- and 3D-BPM

These C-written codes implement the Beam Propagation Method (BPM) to simulate Micro Photonic Systems.

For 2D systems, the Finite Difference Beam Propagation Method (FDBPM) and the Fast Fourier Transform Propagation Method (FFTBPM) might be helpful.
For 3D systems, the open-source C library [fftw3](https://www.fftw.org/) was used for a 3D-FFTBPM implementation. The result is plotted with [Gnuplot](http://www.gnuplot.info) and saved as a PNG image.

As an example, a photonic 1D homogeneous array and the expected discrete diffraction pattern [[1]](https://doi.org/10.1016/j.physrep.2008.04.004):  

![discrete-diffraction](https://user-images.githubusercontent.com/48524687/218358459-cbeb9713-3dad-415b-b803-cf3f650bcf3c.png)

