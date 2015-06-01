#!/bin/sh

 g++  src/roofline.cc build/libsfft.a ../perfplot/pcm/MeasuringCore.a -I../perfplot/pcm `pkg-config fftw3 --libs --cflags` -fopenmp -lippvm -lipps -o roofline-test -march=native -ffast-math -O3


