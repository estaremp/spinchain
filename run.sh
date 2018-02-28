#!/bin/bash

#Clean
rm -rf output

#Compile files
make

#Run executable
./run

#Plot stuff
echo ' >> Plot'
#Read variables to be used
source 'graphical.data'

mv src/Plotting/*.pyc ./bin

python2 bin/probabilities.pyc $VECTORS
python2 bin/eigenvalues.pyc $VECTORS
python2 bin/dynamics.pyc $TOTALTIME $INITIALVEC
python2 bin/exmap.pyc

mkdir output
mkdir output/Plots
mv ./*.data ./*.out ./output
mv ./*.png ./output/Plots

#Clean everything
make clean
rm -rf ./bin


