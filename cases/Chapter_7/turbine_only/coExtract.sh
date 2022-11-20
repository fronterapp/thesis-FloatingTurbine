#!/bin/bash

file=log.mpirun
printf 'Time \t Delta t \t Co \t ICo \n' > Co.dat
my_output="$(awk '/Courant Number/{Co=$6;getline;ICo=$7;getline;DTIME=$3;getline;TIME=$3;printf("%8.6f\t%8.6f\t%8.6f\t%8.6f\n",TIME,DTIME,Co,ICo);}' $file)"
echo "$my_output" >> Co.dat