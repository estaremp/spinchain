#!/bin/bash

#Clean
rm -rf output

#Compile files
make

#Run executable
./run

#Read variables to be used
source 'graphical.data'

mkdir output

if [ $GRAPHICAL = 'T' ] ;
then
    #Plot stuff
    echo ' >> Plot'

    mv src/Plotting/*.pyc ./bin

    python2 bin/probabilities.pyc $VECTORS
    python2 bin/eigenvalues.pyc $VECTORS
    python2 bin/dynamics.pyc $TOTALTIME $INITIALVEC
    python2 bin/exmap.pyc

    mkdir output/Plots
    mv ./*.png ./output/Plots
fi

mv ./*.data ./*.out ./output

#Clean everything
make clean
rm -rf ./bin


